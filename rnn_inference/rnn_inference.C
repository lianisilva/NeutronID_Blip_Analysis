  #include "blip_utils.h"
  #include <iostream>
  #include "TFile.h"
  #include "TTree.h"
  #include "TVector3.h"
  #include "TMath.h"
  
  const char* model_path = "model_traced.pt";
  
  // Global branch arrays
  int nblips;
  int nhits;
  int nparticles;
  //Int_t blip_id[100000];
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
  Int_t hit_g4id[100000];
  
  void rnn_inference(const char* input_file, const char* output_file = "output.root", const char* model_path = "model_traced.pt") {
    
    std::cout << "=== ROOT Macro RNN Inference ===" << std::endl;
    std::cout << "Input: " << input_file << std::endl;
    std::cout << "Output: " << output_file << std::endl;
  
    // Load the shared library
    if(gSystem->Load("./libbliputils.so") < 0) {
      std::cerr << "Failed to load libbliputils.so" << std::endl;
      return;
    }
  
    std::cout << "Initializing model..." << std::endl;
    if(!initializeModel(model_path)) {
      std::cerr << "Failed to initialize model, exiting." << std::endl;
      return;
    }
    std::cout << "Model initialized successfully!" << std::endl;
  
    // Open input file
    TFile* file = new TFile(input_file, "READ");
    if(!file || file->IsZombie()) {
      std::cerr << "Cannot open input file: " << input_file << std::endl;
      return;
    }
  
    std::string fTreeName = "blipana/anatree";
    TTree* tree = (TTree*)file->Get(fTreeName.c_str());
    if(!tree) {
      std::cerr << "Cannot find tree " << fTreeName << std::endl;
      file->Close();
      return;
    }
  
    // Set branch addresses
    tree->SetBranchAddress("nblips", &nblips);
    tree->SetBranchAddress("nhits", &nhits);
    tree->SetBranchAddress("nparticles", &nparticles);
    tree->SetBranchAddress("blip_g4id", &blip_g4id);
    tree->SetBranchAddress("blip_edepid", &blip_edepid);
    tree->SetBranchAddress("hit_blipid", &hit_blipid);
    tree->SetBranchAddress("hit_rms", &hit_rms);
    tree->SetBranchAddress("hit_amp", &hit_amp);
    tree->SetBranchAddress("hit_wire", &hit_wire);
    tree->SetBranchAddress("hit_plane", &hit_plane);
    tree->SetBranchAddress("hit_driftT", &hit_driftT);
    tree->SetBranchAddress("hit_g4id", &hit_g4id);
    tree->SetBranchAddress("blip_g4id", &blip_g4id);
    tree->SetBranchAddress("part_g4id", &part_g4id);
    tree->SetBranchAddress("part_Px", &part_Px);
    tree->SetBranchAddress("part_Py", &part_Py);
    tree->SetBranchAddress("part_Pz", &part_Pz);

    // Make sure to set the hit array pointers
    setHitArrays(hit_rms, hit_amp, hit_wire, hit_plane, hit_driftT);

    // Create output file and tree
    TFile* outFile = TFile::Open(output_file, "RECREATE");
    TTree* outTree = new TTree("rnn_results", "RNN Predictions");
  
    float reco_dir_x, reco_dir_y, reco_dir_z, angular_res, seq_len;
    int blip_index, blip_event;
  
    outTree->Branch("reco_dir_x", &reco_dir_x);
    outTree->Branch("reco_dir_y", &reco_dir_y);
    outTree->Branch("reco_dir_z", &reco_dir_z);
    outTree->Branch("blip_index", &blip_index);
    outTree->Branch("blip_event", &blip_event);
    outTree->Branch("angular_res", &angular_res);
    outTree->Branch("seq_len", &seq_len);
  
    // Process events
    Long64_t nEvents = tree->GetEntries();
    std::cout << "Processing " << nEvents << " events..." << std::endl;
  
    for(int iEvent = 0; iEvent < nEvents; iEvent++) {
      tree->GetEntry(iEvent);
      
      if(iEvent % 100 == 0) {
        std::cout << "Processing event " << iEvent << "/" << nEvents << std::endl;
        std::cout << "Number of blips: " << nblips << std::endl;
      }
  
      // Build blipID -> hitID map
      std::map<int, std::vector<int>> map_blipID_hitID;
      for(int ih = 0; ih < nhits; ih++) { if(hit_g4id[ih] >= 0) { map_blipID_hitID[hit_g4id[ih]].push_back(ih); } }
  
      // Build g4id -> particle index map
      std::map<int, int> map_g4id_index;
      for(int i = 0; i < nparticles; i++) { map_g4id_index[part_g4id[i]] = i; }
  
      // Process each blip
      for(int j = 0; j < nblips; j++) {
	if( map_g4id_index.count(blip_g4id[j]) == 0 ) continue;
        int p = map_g4id_index[blip_g4id[j]];
        TVector3 true_P(part_Px[p], part_Py[p], part_Pz[p]);
  
        // Call the blipDirection function from the library
        std::vector<float> recoDirComponents = blipDirection(
	  blip_g4id[j], 
	  map_blipID_hitID, 
	  file,
	  model_path
        );
  
        TVector3 recoDir(recoDirComponents[0], recoDirComponents[1], recoDirComponents[2]);
  
        reco_dir_x = recoDirComponents[0];
        reco_dir_y = recoDirComponents[1];
        reco_dir_z = recoDirComponents[2];
        seq_len = map_blipID_hitID[blip_g4id[j]].size();

	// Calculate the angular resolution
        angular_res = recoDir.Angle(true_P);
        angular_res *= TMath::RadToDeg();
        if(angular_res > 90) angular_res = 180 - angular_res;
  
        blip_index = j;
        blip_event = iEvent;

	// Fill output tree branches for every blip
        outTree->Fill();
      }
    }
  
    // Save results
    outTree->Write();
    outFile->Close();
    file->Close();
  
    std::cout << "Processing complete!" << std::endl;
    std::cout << "Results saved to: " << output_file << std::endl;
  }
