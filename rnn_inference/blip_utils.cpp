  #include "blip_utils.h"
  #include <iostream>
  #include <algorithm>
  #include <memory>
  
  // Include ALL ROOT headers FIRST
  #include "TTree.h"
  #include "TMath.h"
  
  // THEN undefine the conflicting macros
  #ifdef ClassDef
  #undef ClassDef
  #endif
  #ifdef ClassImp
  #undef ClassImp
  #endif
  
  // LASTLY include PyTorch headers
  #include <torch/torch.h>
  #include <torch/script.h>

 // Global pointers to hit arrays (set via setHitArrays)
  namespace HitData {
    Float_t* hit_rms = nullptr;
    Float_t* hit_amp = nullptr;
    Int_t* hit_wire = nullptr;
    Int_t* hit_plane = nullptr;
    Float_t* hit_driftT = nullptr;
  }

  // Initialization function
  void setHitArrays(Float_t* rms, Float_t* amp, Int_t* wire, 
		    Int_t* plane, Float_t* driftT) {
    HitData::hit_rms = rms;
    HitData::hit_amp = amp;
    HitData::hit_wire = wire;
    HitData::hit_plane = plane;
    HitData::hit_driftT = driftT;
  }

  // RNNPredictor class definition
  class RNNPredictor {
  private:
      torch::jit::script::Module model;
      torch::Device device;
      bool model_loaded;
      std::string loaded_model_path;
   
  public:
    RNNPredictor() : device(torch::kCPU), model_loaded(false) {}
    
    bool loadModel(const std::string& model_path) {
      // Skip if already loaded with same path
      if (model_loaded && loaded_model_path == model_path) {
        std::cout << "Model already loaded, skipping..." << std::endl;
        return true;
      }
      
      try {
	model = torch::jit::load(model_path);
	model.to(device);
	model.eval();
	model_loaded = true;
	loaded_model_path = model_path;
	std::cout << "Model loaded: " << model_path << std::endl;
	return true;
      } catch (const std::exception& e) {
        std::cerr << "Failed to load model: " << e.what() << std::endl;
        return false;
      }
    }

    // Direction Prediction Function
    // Paramaters:
    //	blip_hit_sequence = a blip vector containing each hit in a blip, each hit is a vector of its 5 hit features
    //	nHits = number of hits this blip has
    std::vector<float> predict(const std::vector<std::vector<float>>& blip_hit_seq, int nHits) {
      if (!model_loaded) {
        std::cerr << "Model not loaded!" << std::endl;
        return {};
      }
      try {
        int nFeatures = 5;
        auto blip = torch::zeros({1, nHits, nFeatures});
        for (int i = 0; i < nHits; i++) {
          for (int j = 0; j < nFeatures; j++) {
            blip[0][i][j] = blip_hit_seq[i][j];
          }
        }
        
        auto length = torch::tensor({nHits}, torch::kLong);
        torch::NoGradGuard no_grad;
        std::vector<torch::jit::IValue> input;
        input.push_back(blip);
        input.push_back(length);
        
        at::Tensor output = model.forward(input).toTensor();
        auto accessor = output.accessor<float, 2>();
        std::vector<float> reco_dir = { accessor[0][0], accessor[0][1], accessor[0][2] };
        return reco_dir;
      } catch (const std::exception& e) {
        std::cerr << "Prediction failed: " << e.what() << std::endl;
        return {};
      }
    }

  };

  // Global static predictor instance
  static RNNPredictor& getPredictor() {
    static RNNPredictor predictor;
    return predictor;
  }

  // Initialize the model once (call this at the start of the macro)
  bool initializeModel(const char* model_path) {
    return getPredictor().loadModel(model_path);
  } 

  // blipDirection implementation - single blip (singleton) predictor
  std::vector<float> blipDirection(
    int blipID, 
    std::map<int,std::vector<int>>& map_blipID_hitID, 
    TFile* input_file,
    const char* model_path
  ) {
    // Get the singleton predictor and ensure model is loaded
    RNNPredictor& predictor = getPredictor();
  
    float reco_dir_x, reco_dir_y, reco_dir_z;
    reco_dir_x = reco_dir_y = reco_dir_z = -9;
  
    std::vector<float> blip_amp;
    std::array<std::vector<float>, 3> blip_wire;   // [plane] -> vector of wires
    std::vector<float> blip_plane, blip_time;
  
    std::vector<std::vector<float>> blip_seq;
    int Nhits = 0;
  
    // Process hits for this blip
    for(auto hi : map_blipID_hitID[blipID]) {
      Nhits++;
  
      // Convert to units of wire spacings
      float rms = HitData::hit_rms[hi] * (0.5 * 1.1) / 3.0;
      float time = HitData::hit_driftT[hi] * (0.5 * 1.1) / 3.0;
      float amp = HitData::hit_amp[hi];
      float wire = HitData::hit_wire[hi];
      float plane = HitData::hit_plane[hi];
  
      blip_amp.push_back(amp);
      blip_wire[(int)plane].push_back(wire);
      blip_plane.push_back(plane);
      blip_time.push_back(time);
  
      // Add hits (with features) to blip input
      std::vector<float> hit_features = { rms, amp, wire, plane, time };
      blip_seq.push_back(hit_features);
    }//end hit loop
  
    // Check if all blip hits are only on one plane
    std::sort(blip_plane.begin(), blip_plane.end());    // is this fucking something???
    int uniqueCount = std::unique(blip_plane.begin(), blip_plane.end()) - blip_plane.begin();
    if(uniqueCount == 1) {
      std::cout << "BAD BLIP! Exiting\n";
      return {-9,-9,-9};
    }

    for(int pl = 0; pl < 3; pl++) {
      if(blip_wire[pl].empty()) continue;
      float max_wire = *std::max_element(blip_wire[pl].begin(), blip_wire[pl].end());
      for(auto& wire_num : blip_wire[pl]) wire_num = max_wire - wire_num;
    }

    // Relative time within each blip (zero-shift scaling)
    float t_min = *std::min_element(blip_time.begin(), blip_time.end());
    for(auto& t : blip_time) { t = t - t_min; }

    std::array<int, 3> plane_idx = {0, 0, 0};
    for(int hi = 0; hi < Nhits; hi++) {
      int wire_plane = blip_plane[hi];
      blip_seq[hi][2] = blip_wire[wire_plane][plane_idx[wire_plane]];
      plane_idx[wire_plane]++;
      blip_seq[hi][4] = blip_time[hi];
    }

    // Sort by time (index 4) ascending
    //std::sort(blip_seq.begin(), blip_seq.end(), 
    //          [](const auto& a, const auto& b) { return a[4] < b[4]; });

    // Make prediction
    if(Nhits > 2) {
      auto prediction = predictor.predict(blip_seq, Nhits);

      if(!prediction.empty()) {
        reco_dir_x = prediction[0];
        reco_dir_y = prediction[1];
        reco_dir_z = prediction[2];
      }
    } else {
      std::cout << "\nWarning: blip has <= 2 hits, direction components will be -9s\n";
      return {-9,-9,-9};
    }
  
    return {reco_dir_x, reco_dir_y, reco_dir_z};

  }


