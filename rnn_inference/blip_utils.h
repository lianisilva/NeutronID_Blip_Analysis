  #ifndef BLIP_UTILS_H
  #define BLIP_UTILS_H
  #include <vector>
  #include <map>
  #include <string>
  #include "TFile.h"
   
  // Forward declaration
  class RNNPredictor;
 
  // Initialize model once at the beginning (call this at start of your macro)
  bool initializeModel(const char* model_path);

  void setHitArrays(Float_t* hit_rms, Float_t* hit_amp, Int_t* hit_wire, 
		    Int_t* hit_plane, Float_t* hit_driftT);

  // Main function - will use singleton predictor internally
  std::vector<float> blipDirection(
    int blipID, 
    std::map<int,std::vector<int>>& map_blipID_hitID, 
    TFile* input_file,
    const char* model_path = "model_traced.pt"
  );
  
  #endif
