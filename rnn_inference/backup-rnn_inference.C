
  // ROOT macro for RNN inference

  #include <iostream>
  #include <cmath>
  #include <vector>
  #include <map>
  #include <unordered_set>
  #include <iterator>
  #include <algorithm>
  #include <memory>

  #include "TROOT.h"
  #include "TFile.h"
  #include "TTree.h"
  #include "TVectorT.h"
  #include "TVector3.h"
  #include "TMath.h"

  #ifdef ClassDef
  #undef ClassDef
  #endif
  #ifdef ClassImp
  #undef ClassImp
  #endif

  #include <torch/torch.h>
  #include <torch/script.h>

  const char* model_path = "model_traced.pt";

 // wrapper class for ROOT macro use; simple interface for complex PyTorch C++ API (libtorch)
  class RNNPredictor {
  private:
    torch::jit::script::Module model;
    torch::Device device;   // what device???
    bool model_loaded;
 
  public:
    RNNPredictor() : device(torch::kCPU), model_loaded(false) {}
 
    bool loadModel(const std::string& model_path) {
      try {   // missing this: model = RNN(input_dim=5, hidden_dim=128, output_dim=3)
        model = torch::jit::load(model_path); model.to(device); model.eval();
        model_loaded = true; //std::cout << "Model loaded: " << model_path << std::endl;
        return true;
      } catch (const std::exception& e) { std::cerr << "Failed to load model: " << e.what() << std::endl; return false; }
    }

    //============================== Predictor Function ==============================// 
    std::vector<float> predict(const std::vector<std::vector<float>>& blip_hit_seq, int nHits) {
      if (!model_loaded) { std::cerr << "Model not loaded!" << std::endl; return {}; }
      try {
	// Convert to tensor [1, seq_len, nFeatures]
	//int seq_len = blip_hit_seq.size();
	int nFeatures = 5;
	auto blip = torch::zeros({1, nHits, nFeatures});   // batch_size, seq_len, features
	for (int i = 0; i < nHits; i++) { for (int j = 0; j < nFeatures; j++) { blip[0][i][j] = blip_hit_seq[i][j]; } }
	
	// lengths tensor
	auto length = torch::tensor({nHits}, torch::kLong);
	
	// inference
	torch::NoGradGuard no_grad;
	std::vector<torch::jit::IValue> input;
	input.push_back(blip);
	input.push_back(length);

	at::Tensor output = model.forward(input).toTensor();
	
	// Convert to vector
	auto accessor = output.accessor<float, 2>();   // why these dimensions??
 	std::vector<float> reco_dir = { accessor[0][0], accessor[0][1], accessor[0][2] };
	//for (int i = 0; i < output.size(1); ++i) { reco_dir.push_back(accessor[0][i]); }
	return reco_dir;   // would it be more efficient to do a bunch of blips at once?
      } catch (const std::exception& e) { std::cerr << "Prediction failed: " << e.what() << std::endl; return {}; }
    }

  };

  // initialize input branches, CHANGE for ntuples (add "if" to work with diff file types)
  int nblips;   // = nullptr;
  int nhits;
  int nparticles;
  Int_t blip_id[100000];
  Int_t blip_edepid[100000];
  Int_t hit_blipid[100000];
  Float_t hit_rms[100000];
  Float_t hit_amp[100000];
  Int_t hit_wire[100000];
  Int_t hit_plane[100000];
  Float_t hit_driftT[100000];

  Float_t part_Px[100000];
  Float_t part_Py[100000];
  Float_t part_Pz[100000];
  Int_t part_g4id[100000];
  Int_t blip_g4id[100000];

  //initialize/reset output variables
  float reco_dir_x, reco_dir_y, reco_dir_z, angular_res, seq_len;

    //============================== Blip Processing Function ==============================// 
    std::vector<float> blipDirection( int blipID, std::map<int,std::vector<int>>& map_blipID_hitID, TFile* input_file ) { 

      //TFile* inFile = new TFile(input_file,"READ");
      //if (!*inFile || *inFile->IsZombie()) { std::cerr << "Cannot open input file: " << *input_file << std::endl; return; } 
      std::string fTreeName = "blipana/anatree";
      TTree* inTree = (TTree*)input_file->Get(fTreeName.c_str());
      if (!inTree) { std::cerr << "Cannot find tree "<< fTreeName << std::endl; input_file->Close(); return {-9,-9,-9}; } 

      inTree->SetBranchAddress("hit_rms",    &hit_rms);
      inTree->SetBranchAddress("hit_amp",    &hit_amp);
      inTree->SetBranchAddress("hit_wire",   &hit_wire);
      inTree->SetBranchAddress("hit_plane",  &hit_plane);
      inTree->SetBranchAddress("hit_driftT", &hit_driftT);   //"hit_peakT",       &hit_peakT

      //std::vector<int> blip_num;
      float rms;
      float amp;
      float wire;
      float plane;   // was int, will cause issues??
      float time;   // was double, will cause issues??
      std::vector<int> event;
      //bool bad_blip; std::vector<std::vector<float>> blip_input;

      // initialize predictor
      RNNPredictor predictor;
      if (!predictor.loadModel(model_path)) { std::cerr << "Cannot load model, exiting." << std::endl; return {-9,-9,-9}; }

      // initialize/reset output variables
      reco_dir_x = reco_dir_y = reco_dir_z = angular_res = seq_len = -9;

      std::vector<float> blip_amp;
      std::vector<float> blip_wire0; std::vector<float> blip_wire1; std::vector<float> blip_wire2;
      std::vector<float> blip_plane; std::vector<float> blip_time;   // double

      std::vector<std::vector<float>> blip_seq; int Nhits=0;
      // for all the hit IDs in the blip ID container pertaining to nblip j
      for(auto hi : map_blipID_hitID[blipID] ) {
	Nhits++;

	// convert to units of wire spacings
	rms = hit_rms[hi]*(0.5*1.1)/3.0; time = hit_driftT[hi]*(0.5*1.1)/3.0;
	amp = hit_amp[hi]; wire = hit_wire[hi]; plane = hit_plane[hi];
	blip_plane.push_back(plane);

	blip_amp.push_back(amp);   // try standardizing !!!
	if( hit_plane[hi] == 0 ) blip_wire0.push_back(hit_wire[hi]);
	if( hit_plane[hi] == 1 ) blip_wire1.push_back(hit_wire[hi]);
	if( hit_plane[hi] == 2 ) blip_wire2.push_back(hit_wire[hi]);
	blip_plane.push_back(plane); blip_time.push_back(time);

	// add hits (w/ features) to blip input, ORDER MUST BE PRESERVED
	std::vector<float> hit_features = { rms, amp, static_cast<float>(hit_wire[hi]), plane, time };
	blip_seq.push_back(hit_features);

      }//end hit loop

      // if all blip hits are (oddly) only on one plane
      std::sort(blip_plane.begin(), blip_plane.end());
      int uniqueCount = std::unique(blip_plane.begin(), blip_plane.end()) - blip_plane.begin(); 
      if(uniqueCount==1) { return {-9,-9,-9}; }   // bad_blip=true // std::cout<<"BAD BLIP! Exiting\n"; 

      // centralize wire# per plane per blip
      std::vector<std::vector<float>> all_wires = {blip_wire0, blip_wire1, blip_wire2};
      std::vector<float> wire_max(3);
      for( int pl=0; pl<3; ++pl ) {
        std::vector<float> blip_wire = all_wires[pl];
        if (!blip_wire.empty()) {
          auto max_it = std::max_element( blip_wire.begin(), blip_wire.end() );
          wire_max[pl] = *max_it;
        }
      }
      double time_mean = TMath::Mean(blip_time.begin(), blip_time.end());
      double time_SD = TMath::StdDev(blip_time.begin(), blip_time.end());

      float time_element; float wire_element; float plane_element;
      for( int i=0; i<Nhits; i++ ) {
        time_element = blip_seq[i][4]; plane_element = blip_seq[i][3]; wire_element = blip_seq[i][2];
        blip_seq[i][4]= (time_element - time_mean)/time_SD;
        blip_seq[i][2] = wire_max[plane_element] - wire_element;
      }

      // sort by 5th feature (time, index 4) descending
      std::sort(blip_seq.begin(), blip_seq.end(), [](const auto& a, const auto& b) { return a[4] < b[4]; });

      //for(auto& a : blip_amp) amp.push_back(a);

      // should be no need for padding since we're not doing batch inference
      seq_len = Nhits;
      // make preselection, preprocess then make prediction & set branches
      if( Nhits > 2 ) { // || bad_blip MAKE THIS A REQUIREMENT AND PRINTOUT

        auto prediction = predictor.predict(blip_seq, Nhits);

        if (!prediction.empty()) {
          TVector3 reco_dir(prediction[0], prediction[1], prediction[2]);   // should I make this a branch???
          reco_dir_x = reco_dir.X(); reco_dir_y = reco_dir.Y(); reco_dir_z = reco_dir.Z();    
        }
      } else { return {-9,-9,-9}; }   //std::cout<<"\nWarning: you are including blips with <=2 hits, their direction components will be -9s\n"; 

      return {reco_dir_x, reco_dir_y, reco_dir_z};
    }

  // Main macro function
  void rnn_inference( const char* input_file, const char* output_file = "output.root" ) {
      
    std::cout << "=== ROOT Macro RNN Inference ===" << std::endl;
    std::cout << "Input: " << input_file << std::endl;
    std::cout << "Output: " << output_file << std::endl;
    std::cout << "Model: " << model_path << std::endl;

    // Open input file and get tree
    TFile* file = new TFile(input_file,"READ");
    if (!file || file->IsZombie()) { std::cerr << "Cannot open input file: " << input_file << std::endl; return; }
    std::string fTreeName = "blipana/anatree";
    TTree* tree = (TTree*)file->Get(fTreeName.c_str());
    if (!tree) { std::cerr << "Cannot find tree "<< fTreeName << std::endl; file->Close(); return; }

    tree->SetBranchAddress("nblips",     &nblips);
    tree->SetBranchAddress("nhits",      &nhits);
    tree->SetBranchAddress("nparticles", &nparticles);
    tree->SetBranchAddress("blip_id",	 &blip_id);
    tree->SetBranchAddress("blip_edepid",&blip_edepid);
    tree->SetBranchAddress("hit_blipid", &hit_blipid);

    tree->SetBranchAddress("blip_g4id",  &blip_g4id);
    tree->SetBranchAddress("part_g4id",  &part_g4id);
    tree->SetBranchAddress("part_Px",    &part_Px);
    tree->SetBranchAddress("part_Py",    &part_Py);
    tree->SetBranchAddress("part_Pz",    &part_Pz);


    // Create output file and tree
    TFile* outFile = TFile::Open(output_file, "RECREATE");
    TTree* outTree = new TTree("rnn_results", "RNN Predictions");

    std::map<int, int> map_g4id_index;
    std::map<int,std::vector<int>> map_blipID_hitID;

    int   blip_index, blip_event;

    TVector3 true_P; int max_hits=0; int p;

    outTree->Branch("reco_dir_x", &reco_dir_x);
    outTree->Branch("reco_dir_y", &reco_dir_y);
    outTree->Branch("reco_dir_z", &reco_dir_z);
    outTree->Branch("blip_index", &blip_index);
    outTree->Branch("blip_event", &blip_event);
    outTree->Branch("angular_res",&angular_res);
    outTree->Branch("seq_len",	  &seq_len);
    
    // Process events
    Long64_t nEvents = tree->GetEntries();
    std::cout << "Processing " << nEvents << " events..." << std::endl;
    
    for (Long64_t iEvent = 0; iEvent < nEvents; iEvent++) {
      tree->GetEntry(iEvent); 
      if (iEvent % 100 == 0) { std::cout<<"Processing event "<<iEvent<<"/"<<nEvents<<std::endl;
			       std::cout<<"Number of blips: "<<nblips<<std::endl; }

      // set blipID map, can easily check >2 hits here IF USING VECTOR
      map_blipID_hitID.clear();
      for(int ih=0; ih<nhits; ih++){ if (hit_blipid[ih] >= 0) map_blipID_hitID[hit_blipid[ih]].push_back(ih); }

      map_g4id_index.clear();
      for(int i=0; i<nparticles; i++){   // particle loop
	map_g4id_index[part_g4id[i]] = i;
      }//end particle loop

      for(int j=0; j<nblips; j++) {   // blip loop

        //if( blip_edepid[j] <= 0 ) continue;   //|| part_isContained[i] != 1 || blip_g4id[j] != part_g4id[p]
	p = map_g4id_index[blip_g4id[j]];
	true_P.SetX(part_Px[p]); true_P.SetY(part_Py[p]); true_P.SetZ(part_Pz[p]);
	std::vector<float> recoDirComponents = blipDirection( blip_id[j], map_blipID_hitID, file );
	TVector3 recoDir(recoDirComponents[0], recoDirComponents[1], recoDirComponents[2]);

        angular_res = recoDir.Angle(true_P);
        angular_res *= TMath::RadToDeg();
        if (angular_res > 90) angular_res = 180 - angular_res; 

        //if( Nhits > max_hits ) max_hits = Nhits;
	blip_index = j; blip_event = iEvent; 

	outTree->Fill();

      }//end blip loop

    }//end event loop

    // Save results
    outTree->Write();
    outFile->Close();
    file->Close();
    
    std::cout << "✓ Processing complete!" << std::endl;
    std::cout << "Results saved to: " << output_file << std::endl;
  }

  #ifndef __CINT__
  int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <input.root> <output.root> [model.pt]" << std::endl;
        return 1;
    }
    
    const char* input_file = argv[1];
    const char* output_file = argv[2];
   
    rnn_inference(input_file, output_file);
    return 0;
  }
  #endif

