
  #include "TROOT.h"
  #include <cmath>
  #include <algorithm>
  #include <iomanip>   // for std::boolalpha
  #include <iostream>
  #include <vector>
  #include <string>
  #include "merge_plot2.C"

  // TMVA headers
  #include "TMVA/RReader.hxx"
  #include "TMVA/Reader.h"
  #include "ROOT/RVec.hxx"

  int maxEvents = -2000;

  TDirectory* blipE_plots;
  TDirectory* recoVtxDist_plots;
  TDirectory* cosAlignment_plots;
  TDirectory* alignmentAngle_plots;
  TDirectory* resolution_plots;
  TDirectory* res_vs_E_plots;
  TDirectory* Edensity_plots;
  TDirectory* byDeadWire_plots;
  TDirectory* proxTrkDist_plots;
  TDirectory* cosmicDistToVtx_plots;
  TDirectory* dx_plots;
  TDirectory* dw_plots;
  TDirectory* size_plots;
  TDirectory* inShower_plots;
  TDirectory* energy_res_plots;
  TDirectory* BDT_score_plots;

  int			   n_pfps;
  //int			   selected;
  //int			   n_tracks;
  //int			   n_showers;
  float			   reco_nu_vtx_sce_x = 0;
  float			   reco_nu_vtx_sce_y = 0;
  float			   reco_nu_vtx_sce_z = 0;
  int			   slice_orig_pass_id = 0;
  bool			   isVtxInFiducial = 0;
  float 		   topological_score = 0;
  int			   crtveto = 0;
  double		   _closestNuCosmicDist = 0;
  //float			   bdt_numuCC = 0;
  int			   nblips_saved;
  std::vector<float>	   * blip_dx = 0;
  std::vector<float>	   * blip_dw = 0;
  std::vector<float>	   * blip_energy = 0;
  std::vector<float>	   * blip_charge = 0;
  std::vector<int>  	   * blip_nplanes = 0;
  std::vector<float>	   * blip_proxtrkdist = 0;
  std::vector<int>  	   * blip_proxtrkid = 0;
  std::vector<bool> 	   * blip_touchtrk = 0;
  std::vector<float>	   * blip_touchtrkid = 0;
  std::vector<int>  	   * blip_true_g4id = 0;
  std::vector<int>  	   * blip_true_pdg = 0;
  std::vector<float>	   * blip_true_energy = 0;
  std::vector<float>	   * blip_x = 0;
  std::vector<float>	   * blip_y = 0;
  std::vector<float>	   * blip_z = 0;
  std::vector<int>	   * blip_pl0_nwires = 0;
  std::vector<int>	   * blip_pl1_nwires = 0;
  std::vector<int>	   * blip_pl2_nwires = 0;
  std::vector<int>	   * blip_pl0_bydeadwire = 0;
  std::vector<int>	   * blip_pl1_bydeadwire = 0;
  std::vector<int>	   * blip_pl2_bydeadwire = 0;
  std::vector<bool>        * blip_rnn_dir_isValid = 0;
  std::vector<float>       * blip_rnn_dir_x = 0;
  std::vector<float>       * blip_rnn_dir_y = 0;
  std::vector<float>       * blip_rnn_dir_z = 0;
  std::vector<float>       * blip_true_dir_x = 0;
  std::vector<float>       * blip_true_dir_y = 0;
  std::vector<float>       * blip_true_dir_z = 0;
  std::vector<float>	   * blip_true_ncategory = 0;

  std::vector<int>	   * pfng2semlabel = 0;

  std::vector<float>	   * trk_sce_start_x_v = 0;
  std::vector<float>	   * trk_sce_start_y_v = 0;
  std::vector<float>	   * trk_sce_start_z_v = 0;
  std::vector<float>	   * trk_sce_end_x_v = 0;
  std::vector<float>	   * trk_sce_end_y_v = 0;
  std::vector<float>	   * trk_sce_end_z_v = 0;
  std::vector<float>	   * trk_llr_pid_score_v = 0;
  //std::vector<float>	   * trk_len = 0;
  std::vector<float>       * trk_score_v = 0;
  //std::vector<float>       * trk_pfp_id_v = 0;

  std::vector<float>	   * shr_start_x_v;
  std::vector<float>	   * shr_start_y_v;
  std::vector<float>	   * shr_start_z_v;
  std::vector<float>	   * shr_px_v;
  std::vector<float>	   * shr_py_v;
  std::vector<float>	   * shr_pz_v;
  std::vector<float>	   * shr_dist_v;
  std::vector<float>       * shr_openangle_v;
  //float       shr_openangle;
  std::vector<float>       * shr_energy_y_v;

  std::vector<int>	   * all_mc_trkid = 0;
  std::vector<int>	   * all_mc_pdg = 0;
  std::vector<float>	   * all_mc_E = 0;
  std::vector<float>	   * all_mc_endx = 0;
  std::vector<float>	   * all_mc_endy = 0;
  std::vector<float>	   * all_mc_endz = 0;
  std::vector<int>	   * all_mc_mother = 0;
  std::vector<std::string> * all_mc_process = 0;

  //Int_t    reco_id[1000];   // wc tree vars are ints not vectors!

  //bool   nu_isNeutrino;
  float  reco_nu_vtx_x;
  float  reco_nu_vtx_y;
  float  reco_nu_vtx_z;

  const int numHists=10;
  TH1D* h_blipE[numHists];
  TH1D* h_recoVtxDist[numHists];
  TH1D* h_cosAlignment[numHists];
  TH1D* h_alignmentAngle[numHists];
  TH1D* h_resolution[numHists];
  TH2D* h_res_vs_E[numHists];
  TH1D* h_Edensity[numHists];
  TH1D* h_byDeadWire[numHists];
  TH1D* h_proxTrkDist[numHists];
  TH1D* h_cosmicDistToVtx[numHists];
  TH1D* h_dx[numHists];
  TH1D* h_dw[numHists];
  TH1D* h_size[numHists];
  TH1D* h_inShower[numHists];
  TH1D* h_energy_res[numHists];
  TH1D* h_BDT_score[numHists];
  TH1D* h_nBlips;
  TH1D* h_BlipsPerEvent;
  TH1D* h_maxWires;
  TH1D* h_nwires_pl0;
  TH1D* h_nwires_pl1;
  TH1D* h_nwires_pl2;
  TH2D* h_trueNeutronEvsProtonE;
  TH2D* h_ShrE_vs_OpAng;
  TH1D* h_trueNeutronE;
  TH1D* h_preFoundNeutronE;
  TH1D* h_foundNeutronE;
  TH1D* h_prTrksPerE;
  TH1D* h_prBlipsPerE;

  TGraph* gr_signal_E_ds;
  TGraph* gr_elec_E_ds;
  TGraph* gr_cosmic_E_ds;

  // Input file
  //std::string fFileName = "checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root";   // retupled BNB nu
  //std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus.root";   // BNB nu with more blip vars
  std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260224.root";   // BNB nu with more blip vars (full events)
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";   // BNB nu (no all_mc)

  std::string nuTreeName = "nuselection/NeutrinoSelectionFilter";
  //std::string wcTreeName = "wcpselection/T_PFeval";

  std::map<int, int> map_allTrkID_partIDX;
  std::map<int, int> map_allTrkID_momTrkID;

  //============================== FUNCTIONS ==============================//
  // Function to determine if a particle descended from another particle.
  // Allows option to break lineage at photons for contiguous parentage.
  bool IsAncestorOf(int particleID, int ancestorID, bool breakAtPhots = false){
    if( particleID == ancestorID  )       return true;
    if( particleID < ancestorID   )       return false;
    if( map_allTrkID_momTrkID.count(ancestorID) == 0 )  return false;

    int iterations = 0; const int MAX_ITERATIONS = 1000; // prevent infinite loops

    while( particleID > ancestorID && iterations < MAX_ITERATIONS ) {
      iterations++;
      if( map_allTrkID_momTrkID.count(particleID) == 0 ) { return false; }

      int mother = map_allTrkID_momTrkID[particleID];
      if( mother == ancestorID ) { return true; } if( mother == 0 ) { return false; }
      //if( !map_g4id_pdg.count(mother) ) { return false; }
      //if ( breakAtPhots == true && map_g4id_pdg[mother] == 22 )      { return false; }
      particleID = mother;
    }
    if( iterations >= MAX_ITERATIONS ) { std::cout << "\nWarning: Hit max iterations in IsAncestorOf - infinite loop?!"; } 
    return false;
  }

  bool greater_than_1(int value){ return value > 1; }

  bool IsInFV( TVector3 point ) {
    return point.X() >= 10.	&&  point.X() <= 246.35  &&
	   point.Y() >= -106.5  &&  point.Y() <= 106.5   &&
	   point.Z() >= 10.	&&  point.Z() <= 1026.80 ;
  }

  bool IsWithinShowerCone(const TVector3& vtx,     // Shower Vtx
			  const TVector3& dir,     // Shower direction (momentum)
			  const TVector3& blip,    // Blip coordinates
			  double openingAngleDeg)  // Total opening angle of exclusion cone 
  {
    TVector3 d = blip - vtx;
    double angle = d.Angle(dir);
    double opang = openingAngleDeg * (TMath::Pi() / 180.) * 0.5;
    return ( angle < opang );
  }

  TFile* f = new TFile("fit_cut.root", "READ"); 
  TGraph* cutoff_graph = (TGraph*)f->Get("cutoff_graph");
  TTree* param_tree = (TTree*)f->Get("param_tree");        
  float A, B, C, D;
  param_tree->SetBranchAddress("Param_A", &A);
  param_tree->SetBranchAddress("Param_B", &B);
  param_tree->SetBranchAddress("Param_C", &C);
  param_tree->SetBranchAddress("Param_D", &D);
  param_tree->GetEntry(0);  // Only one entry with the parameters
  f->Close();
  float Fit_Function_PID_Cut( float blip_recoE ){ 
    float x = blip_recoE;
    return A * TMath::Log( B * x ) + C*x + D ;
  }

  //============================== LOAD BDT ==============================//
  //TMVA::Experimental::RReader model("n_p_select_noShow_100k.xml");
  TMVA::Reader *reader = new TMVA::Reader( );

  // Provide the features exactly as used during training
  //std::vector<std::string> featureNames = {"blip_E", "dist2Vtx", "blip_X", "blip_Y", "blip_Z", "blip_dX", "blip_dW", "blip_size", "blip_density"};
  float blip_E, dist2Vtx, blip_X, blip_Y, blip_Z, blip_dX, blip_dW, blip_size, blip_density;//, blip_dir; 
  reader->AddVariable("blip_E", &blip_E);
  reader->AddVariable("dist2Vtx", &dist2Vtx);
  reader->AddVariable("blip_X", &blip_X);
  reader->AddVariable("blip_Y", &blip_Y);
  reader->AddVariable("blip_Z", &blip_Z);
  reader->AddVariable("blip_dX", &blip_dX);
  reader->AddVariable("blip_dW", &blip_dW);
  reader->AddVariable("blip_size", &blip_size);
  reader->AddVariable("blip_density", &blip_density);
  //reader->AddVariable("blip_dir", &blip_dir);

  reader->BookMVA("BDT", "models_BDT/n_p_select_withDir.xml");

  // Set BDT variables directly from a saved file instead of importing entire model //
  // This process is slower and more involved, only use if there's difficulty with TMVA import //
  //TFile* BDT_results = new TFile("blipPlus_quasiFull_test_noShow_100k.root");
  //TTree* results = (TTree*)BDT_results->Get("results");
  //Double_t iEvent; Double_t test_blip_idx; Float_t prediction;
  //results->SetBranchAddress("iEvent", &iEvent);
  //results->SetBranchAddress("test_blip_idx", &test_blip_idx);
  //results->SetBranchAddress("prediction", &prediction);
  //results->GetEntry(0);
  //std::map<pair<int, int>, float> map_iEvent_blipIDX_to_pred;
  //for( int i = 0; i < results->GetEntries(); i++ ) {
  //  results->GetEntry(i);
  //  map_iEvent_blipIDX_to_pred[std::make_pair((int)iEvent, (int)test_blip_idx)] = prediction;
  //}


  // ********************************************************
  void surprise_nu_analysis(){
    //TH1::AddDirectory(kFALSE);   // remove automatic histogram registration to the temporary gDirectory

    // Read in the TFile and grab the TTree
    TFile* file = new TFile(("files/"+fFileName).c_str(),"READ");
    TTree* nuTree = (TTree*)file->Get(nuTreeName.c_str());
    //TTree* wcTree = (TTree*)file->Get(wcTreeName.c_str());
    std::cout<<"Reading data file "<<fFileName<<"\n\n";


    // Set branch address for all variables you plan to use
    //nuTree->SetBranchAddress("selected",	    &selected);
    //nuTree->SetBranchAddress("n_tracks",	    &n_tracks);
    nuTree->SetBranchAddress("nblips_saved",	    &nblips_saved);
    nuTree->SetBranchAddress("blip_true_g4id",	    &blip_true_g4id);
    nuTree->SetBranchAddress("blip_x",		    &blip_x);
    nuTree->SetBranchAddress("blip_y",		    &blip_y);
    nuTree->SetBranchAddress("blip_z",		    &blip_z);
    nuTree->SetBranchAddress("blip_dx",             &blip_dx);
    nuTree->SetBranchAddress("blip_dw",             &blip_dw);
    nuTree->SetBranchAddress("blip_energy",	    &blip_energy);
    //nuTree->SetBranchAddress("blip_charge",         &blip_charge);
    //nuTree->SetBranchAddress("blip_nplanes",        &blip_nplanes);
    nuTree->SetBranchAddress("blip_proxtrkdist",    &blip_proxtrkdist);
    nuTree->SetBranchAddress("blip_touchtrk",       &blip_touchtrk);
    nuTree->SetBranchAddress("blip_true_pdg",       &blip_true_pdg);
    nuTree->SetBranchAddress("blip_true_energy",    &blip_true_energy);

    //nuTree->SetBranchAddress("blip_pl0_nwires",	    &blip_pl0_nwires);
    //nuTree->SetBranchAddress("blip_pl1_nwires",	    &blip_pl1_nwires);
    //nuTree->SetBranchAddress("blip_pl2_nwires",	    &blip_pl2_nwires);
    nuTree->SetBranchAddress("blip_pl0_bydeadwire", &blip_pl0_bydeadwire);
    nuTree->SetBranchAddress("blip_pl1_bydeadwire", &blip_pl1_bydeadwire);
    nuTree->SetBranchAddress("blip_pl2_bydeadwire", &blip_pl2_bydeadwire);

    nuTree->SetBranchAddress("blip_rnn_dir_isValid",&blip_rnn_dir_isValid);
    nuTree->SetBranchAddress("blip_rnn_dir_x",	    &blip_rnn_dir_x);
    nuTree->SetBranchAddress("blip_rnn_dir_y",	    &blip_rnn_dir_y);
    nuTree->SetBranchAddress("blip_rnn_dir_z",	    &blip_rnn_dir_z);
    nuTree->SetBranchAddress("blip_true_dir_x",	    &blip_true_dir_x);
    nuTree->SetBranchAddress("blip_true_dir_y",	    &blip_true_dir_y);
    nuTree->SetBranchAddress("blip_true_dir_z",	    &blip_true_dir_z);

    nuTree->SetBranchAddress("n_pfps",		    &n_pfps);
    nuTree->SetBranchAddress("pfng2semlabel",	    &pfng2semlabel);

    nuTree->SetBranchAddress("trk_sce_start_x_v",   &trk_sce_start_x_v);
    nuTree->SetBranchAddress("trk_sce_start_y_v",   &trk_sce_start_y_v);
    nuTree->SetBranchAddress("trk_sce_start_z_v",   &trk_sce_start_z_v);
    nuTree->SetBranchAddress("trk_sce_end_x_v",     &trk_sce_end_x_v);
    nuTree->SetBranchAddress("trk_sce_end_y_v",     &trk_sce_end_y_v);
    nuTree->SetBranchAddress("trk_sce_end_z_v",     &trk_sce_end_z_v);
    nuTree->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v);
    //nuTree->SetBranchAddress("trk_len",		    &trk_len);
    nuTree->SetBranchAddress("trk_score_v",	    &trk_score_v);
    //nuTree->SetBranchAddress("trk_pfp_id_v",	    &trk_pfp_id_v);

    //nuTree->SetBranchAddress("n_showers",           &n_showers);
    nuTree->SetBranchAddress("shr_start_x_v",	    &shr_start_x_v);
    nuTree->SetBranchAddress("shr_start_y_v",	    &shr_start_y_v);
    nuTree->SetBranchAddress("shr_start_z_v",	    &shr_start_z_v);
    nuTree->SetBranchAddress("shr_px_v",	    &shr_px_v);
    nuTree->SetBranchAddress("shr_py_v",	    &shr_py_v);
    nuTree->SetBranchAddress("shr_pz_v",	    &shr_pz_v);
    nuTree->SetBranchAddress("shr_dist_v",	    &shr_dist_v);
    nuTree->SetBranchAddress("shr_openangle_v",     &shr_openangle_v);
    //nuTree->SetBranchAddress("shr_openangle",       &shr_openangle);
    nuTree->SetBranchAddress("shr_energy_y_v",	    &shr_energy_y_v);

    nuTree->SetBranchAddress("all_mc_trkid",	    &all_mc_trkid);
    nuTree->SetBranchAddress("all_mc_pdg",	    &all_mc_pdg);
    nuTree->SetBranchAddress("all_mc_E",	    &all_mc_E);
    //nuTree->SetBranchAddress("all_mc_endx",         &all_mc_endx);
    //nuTree->SetBranchAddress("all_mc_endy",         &all_mc_endy);
    //nuTree->SetBranchAddress("all_mc_endz",         &all_mc_endz);
    nuTree->SetBranchAddress("all_mc_mother",       &all_mc_mother);
    nuTree->SetBranchAddress("all_mc_process",      &all_mc_process);

    nuTree->SetBranchAddress("reco_nu_vtx_sce_x",   &reco_nu_vtx_sce_x);
    nuTree->SetBranchAddress("reco_nu_vtx_sce_y",   &reco_nu_vtx_sce_y);
    nuTree->SetBranchAddress("reco_nu_vtx_sce_z",   &reco_nu_vtx_sce_z);
    nuTree->SetBranchAddress("slice_orig_pass_id",  &slice_orig_pass_id);
    nuTree->SetBranchAddress("isVtxInFiducial",     &isVtxInFiducial);
    nuTree->SetBranchAddress("topological_score",   &topological_score);
    nuTree->SetBranchAddress("crtveto",		    &crtveto);
    nuTree->SetBranchAddress("_closestNuCosmicDist",&_closestNuCosmicDist);
    //nuTree->SetBranchAddress("bdt_numuCC",	    &bdt_numuCC);
    //nuTree->SetBranchAddress("reco_nu_vtx_x",       &reco_nu_vtx_x);
    //nuTree->SetBranchAddress("reco_nu_vtx_y",       &reco_nu_vtx_y);
    //nuTree->SetBranchAddress("reco_nu_vtx_z",       &reco_nu_vtx_z);

    //wcTree->SetBranchAddress("reco_id",		    &reco_id);

    // make output file to before declaring plots
    TFile* outfile = new TFile("test_out.root", "recreate");
    blipE_plots           = outfile->mkdir("blipE_plots");
    recoVtxDist_plots     = outfile->mkdir("recoVtxDist_plots");
    cosAlignment_plots    = outfile->mkdir("cosAlignment_plots");
    alignmentAngle_plots  = outfile->mkdir("alignmentAngle_plots");
    resolution_plots      = outfile->mkdir("resolution_plots");
    res_vs_E_plots        = outfile->mkdir("res_vs_E_plots");
    Edensity_plots        = outfile->mkdir("Edensity_plots");
    byDeadWire_plots      = outfile->mkdir("byDeadWire_plots");
    proxTrkDist_plots     = outfile->mkdir("proxTrkDist_plots");
    cosmicDistToVtx_plots = outfile->mkdir("cosmicDistToVtx_plots");
    dx_plots              = outfile->mkdir("dx_plots");
    dw_plots              = outfile->mkdir("dw_plots");
    size_plots            = outfile->mkdir("size_plots");
    inShower_plots	  = outfile->mkdir("inShower_plots");
    energy_res_plots	  = outfile->mkdir("energy_res_plots");
    BDT_score_plots	  = outfile->mkdir("BDT_score_plots");

    // Make histograms
    blipE_plots->cd();
    h_blipE[0]		= new TH1D("blipE","Reconstructed Blip Energy;Energy [MeVee]",30,0,90);//20,0,100);
    for(int i=1; i<numHists; i++) { h_blipE[i] = (TH1D*)h_blipE[0]->Clone(("blipE_" + std::to_string(i)).c_str()); }
    recoVtxDist_plots->cd();
    h_recoVtxDist[0]	= new TH1D("recoVtxDist","Blip Distance to Reconstructed Neutrino Vertex;Distance [cm];Nblips",20,0,100);
    for(int i=1; i<numHists; i++) { h_recoVtxDist[i] = (TH1D*)h_recoVtxDist[0]->Clone(("recoVtxDist_" + std::to_string(i)).c_str()); }
    cosAlignment_plots->cd();
    h_cosAlignment[0]        = new TH1D("cosAlignment","|Cos(#Theta_{pRT})| Between Reconstructed Direction and Neutrino Vertex;|Cos(#Theta_{pRT})|",20,0,1);//40,-1,1);
    for(int i=1; i<numHists; i++) { h_cosAlignment[i] = (TH1D*)h_cosAlignment[0]->Clone(("cosAlignment_" + std::to_string(i)).c_str()); }
    alignmentAngle_plots->cd();
    h_alignmentAngle[0]      = new TH1D("alignmentAngle","Alignment Angle Between Reconstructed Direction and Neutrino Vertex;#Theta_{pRT} [degrees]",30,0,90);
    for(int i=1; i<numHists; i++) { h_alignmentAngle[i] = (TH1D*)h_alignmentAngle[0]->Clone(("alignmentAngle_" + std::to_string(i)).c_str()); }
    res_vs_E_plots->cd();
    h_res_vs_E[0]          = new TH2D("res_vs_E","Angular Resolution vs. True Blip Energy;True Energy [MeV];#Theta_{pRT} [degrees]",20,0,100,18,0,90);
    for(int i=1; i<numHists; i++) { h_res_vs_E[i] = (TH2D*)h_res_vs_E[0]->Clone(("res_vs_E_" + std::to_string(i)).c_str()); }
    resolution_plots->cd();
    h_resolution[0]          = new TH1D("resolution","Blip Angular Resolution;#Theta_{pRT} [degrees]",30,0,90);//46,-1,91);
    for(int i=1; i<numHists; i++) { h_resolution[i] = (TH1D*)h_resolution[0]->Clone(("resolution_" + std::to_string(i)).c_str()); }
    Edensity_plots->cd();
    h_Edensity[0]            = new TH1D("Edensity","Blip Energy Deposition Density;E/ds [MeVee/cm]",30,0,15);
    for(int i=1; i<numHists; i++) { h_Edensity[i] = (TH1D*)h_Edensity[0]->Clone(("Edensity_" + std::to_string(i)).c_str()); }
    byDeadWire_plots->cd();
    h_byDeadWire[0]          = new TH1D("byDeadWire","Blip By Dead Wire",2,0,2);
    for(int i=1; i<numHists; i++) { h_byDeadWire[i] = (TH1D*)h_byDeadWire[0]->Clone(("byDeadWire" + std::to_string(i)).c_str()); }
    proxTrkDist_plots->cd();
    h_proxTrkDist[0]         = new TH1D("proxTrkDist","Blip Distance to Nearest Track;Reconstructed Distance [cm]",21,0,105);
    for(int i=1; i<numHists; i++) { h_proxTrkDist[i] = (TH1D*)h_proxTrkDist[0]->Clone(("proxTrkDist_" + std::to_string(i)).c_str()); }
    cosmicDistToVtx_plots->cd();
    h_cosmicDistToVtx[0]      = new TH1D("cosmicDistToVtx","Reco. Neutrino Vertex Distance to Nearest Cosmic Track;Reconstructed Distance [cm]",41,0,510);
    for(int i=1; i<numHists; i++) { h_cosmicDistToVtx[i] = (TH1D*)h_cosmicDistToVtx[0]->Clone(("cosmicDistToVtx_" + std::to_string(i)).c_str()); }
    dx_plots->cd();
    h_dx[0]		= new TH1D("dx","Blip dX;dX [cm]",24,0,6);
    for(int i=1; i<numHists; i++) { h_dx[i] = (TH1D*)h_dx[0]->Clone(("dx_" + std::to_string(i)).c_str()); }
    dw_plots->cd();
    h_dw[0]		= new TH1D("dw","Blip dW;dW [cm]",20,0,5);
    for(int i=1; i<numHists; i++) { h_dw[i] = (TH1D*)h_dw[0]->Clone(("dw_" + std::to_string(i)).c_str()); }
    size_plots->cd();
    h_size[0]		= new TH1D("size","Blip Size;3D Projected Size [cm]",28,0,7);
    for(int i=1; i<numHists; i++) { h_size[i] = (TH1D*)h_size[0]->Clone(("size_" + std::to_string(i)).c_str()); }
    inShower_plots->cd();
    h_inShower[0]           = new TH1D("inShower","Blip Inside Shower;In or Out of Shower;nBlips",2,0,2);
    for(int i=1; i<numHists; i++) { h_inShower[i] = (TH1D*)h_inShower[0]->Clone(("inShower_" + std::to_string(i)).c_str()); }
    energy_res_plots->cd();
    h_energy_res[0]	    = new TH1D("energy_res","Blip Energy Resolution;(recoE-trueE)/trueE;nBlips",20,0,100);
    for(int i=1; i<numHists; i++) { h_energy_res[i] = (TH1D*)h_energy_res[0]->Clone(("energy_res_" + std::to_string(i)).c_str()); }
    BDT_score_plots->cd();
    h_BDT_score[0]	    = new TH1D("BDT_score","Blip BDT Prediction;BDT Score (1 is signal-like);nBlips",70,-12,2);
    for(int i=1; i<numHists; i++) { h_BDT_score[i] = (TH1D*)h_BDT_score[0]->Clone(("BDT_score_" + std::to_string(i)).c_str()); }

    //h_mostWiresPlane[]
    //  std::vector<int>         * blip_pl0_nwires = 0;
    //  std::vector<int>         * blip_pl1_nwires = 0;
    //  std::vector<int>         * blip_pl2_nwires = 0;

    outfile->cd();
    h_BlipsPerEvent	= new TH1D("BlipsPerEvent","Number of Blips per Event;# Blips",250,0,500);
    h_nwires_pl0	= new TH1D("nwires_pl0","Number of Wires per Blip on Plane 0;# Wires;# Blips",18,0,18);
    h_nwires_pl1	= new TH1D("nwires_pl1","Number of Wires per Blip on Plane 1;# Wires;# Blips",18,0,18);
    h_nwires_pl2	= new TH1D("nwires_pl2","Number of Wires per Blip on Plane 2;# Wires;# Blips",18,0,18);
    h_maxWires		= new TH1D("maxWires","Max Wires per Blip on one Plane;# Wires",18,0,18);
    h_trueNeutronEvsProtonE = new TH2D("trueNeutronEvsProtonE","True Proton Energy vs True Neutron Energy;Primary Neutron Kinetic Energy [MeV];Proton Kinetic Energy [MeV]",60,0,300,60,0,300);
    h_ShrE_vs_OpAng     = new TH2D("ShrE_vs_OpAng","Shower Energy vs. Opening Angle;Opening Angle [radians];Shower Energy on Collection Plane",157,0,1.57,200,0,1000);
    h_trueNeutronE      = new TH1D("trueNeutronE","True Neutron Energy;MeVee",600,0,600); h_trueNeutronE->Sumw2();
    h_preFoundNeutronE  = new TH1D("preFoundNeutronE","Preselection Identified Neutron Energy;MeVee",600,0,600); h_preFoundNeutronE->Sumw2();
    h_foundNeutronE     = new TH1D("foundNeutronE","Identified Neutron Energy;MeVee",600,0,600); h_foundNeutronE->Sumw2();
    h_prTrksPerE	= new TH1D("prTrksPerE", "Number of Reconstructed Objects vs True Energy of Protons From Primary Neutrons;Proton Kinetic Energy [MeV];Number of Reconstructed Objects",100,0,200);
    h_prBlipsPerE	= new TH1D("prBlipsPerE","Number of Reconstructed Objects vs True Energy of Protons From Primary Neutrons;Proton Kinetic Energy [MeV];Number of Reconstructed Objects",100,0,200);

    gr_signal_E_ds = new TGraph(); gr_signal_E_ds->SetNameTitle("signal_E_ds", "Blip Energy Deposition Density;Reconstructed Blip Energy (E_{blip}) [MeVee];E_{blip}/ds [MeVee/cm]");
    gr_elec_E_ds = new TGraph(); gr_elec_E_ds->SetNameTitle("gr_elec_E_ds", "Blip Energy Deposition Density;Reconstructed Blip Energy (E_{blip}) [MeVee];E_{blip}/ds [MeVee/cm]");
    gr_cosmic_E_ds = new TGraph(); gr_cosmic_E_ds->SetNameTitle("cosmic_E_ds", "Blip Energy Deposition Density;Reconstructed Blip Energy (E_{blip}) [MeVee];E_{blip}/ds [MeVee/cm]");

    float minX, minY, minZ = 0;
    float maxX, maxY, maxZ = 0;
 
    int numEvents = nuTree->GetEntries();
    std::cout<<"There are "<<numEvents<<" events in the nu tree.\n";
    //std::cout<<"There are "<<wcTree->GetEntries()<<" events in the WC tree.\n";  

    int totalBlips=0;
    float totalBlipE=0;
    int totalSelectedBlips=0;   // can't be float due to floating-point precision limit
    float matchedBlips=0;
    float neutronIDprBlips=0;
    float signal_blips=0;
    int selectedEvents=0;
    int nonSelected=0;
    int numEvents_numuCC=0;
    std::vector<int> nShrs_per_numuCC_event;

    int total_nu_Ns=0;
    std::map<int, std::vector<int>> map_n_p;

    int   category;

    //============================== EVENT LOOP ==============================//
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      //if( iEvent > 10000 ) continue;
      if( iEvent < 100000  && iEvent % 2 == 0 ) continue; selectedEvents++;
      nuTree->GetEntry(iEvent);

      int all_mc_length = static_cast<int>(all_mc_trkid->size());

      if(iEvent < 10 || iEvent % 10000 == 0) {
	std::cout<<"\n=======================================\n";
	std::cout<<"Reading event "<<iEvent<<" out of "<<numEvents<<"\n";
	std::cout<<"NBlips: "<<nblips_saved<<"   NParticles (all_mc): "<<all_mc_length<<"\n";
      }

      map_allTrkID_partIDX.clear();
      map_allTrkID_momTrkID.clear();
      map_n_p.clear();

      TVector3 reco_nu_vtx( reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z );
      float cosmicProx = _closestNuCosmicDist;

      //pfp_slice_idx
      //n_tracks
      //pfng2bkgavrg   //Average cosmic score of hits in the PFParticle

      //============================== NUMUCC INCLUSIVE SELECTION ==============================//
      // Preselection
      if( slice_orig_pass_id != 1 ||   // replaces nslice != 1 ||
	     isVtxInFiducial != 1 ||
	topological_score <= 0.06 ||
		     crtveto == 1 ||
	_closestNuCosmicDist <= 5 ) { continue; }
	       //bdt_numuCC > 0.2 || bdt_numuCC > 10 ) { continue; }

      //============================== PFP LOOP ==============================//
      int n_shrs = 0;
      bool numuCC = false;
      std::vector<std::tuple<TVector3, TVector3, float>> shr_vtx_P_E; //std::vector<std::pair<TVector3, TVector3>> ShrVtx_ShrP;
      for(size_t i_pfp = 0; i_pfp < n_pfps; i_pfp++){

	// Muon Track   (could we assume this is always the leading track?)
	if( trk_score_v->at(i_pfp) > 0.45 ) { 
	  TVector3 trkStart(trk_sce_start_x_v->at(i_pfp), trk_sce_start_y_v->at(i_pfp), trk_sce_start_z_v->at(i_pfp));
	  TVector3 trkEnd(trk_sce_end_x_v->at(i_pfp), trk_sce_end_y_v->at(i_pfp), trk_sce_end_z_v->at(i_pfp));
	  float trk_length = (trkEnd - trkStart).Mag();

	  if( IsInFV(trkStart) && IsInFV(trkEnd) &&
	      std::min( (trkStart-reco_nu_vtx).Mag(), (trkEnd-reco_nu_vtx).Mag() ) < 4 &&
	      trk_llr_pid_score_v->at(i_pfp) > 0.2 &&
	      trk_length > 10 && 
	      trk_score_v->at(i_pfp) > 0.8)
	  { numuCC = true; }//trk_id.push_back(trk_pfp_id_v->at(i_pfp));

	// Showers
	} else { //if( pfng2semlabel->at(i_pfp) == 2 ) {
	  TVector3 shrStart(shr_start_x_v->at(i_pfp), shr_start_y_v->at(i_pfp), shr_start_z_v->at(i_pfp));
	  TVector3 shr_p(shr_px_v->at(i_pfp), shr_py_v->at(i_pfp), shr_pz_v->at(i_pfp));   // already normalized
	  TVector3 shr_distVec = shr_p * shr_dist_v->at(i_pfp);   // ( |Px|*d, |Py|*d, |Pz|*d )
	  TVector3 shrEnd = shrStart + shr_distVec; float shrE = shr_energy_y_v->at(i_pfp);
	  shr_vtx_P_E.emplace_back(shrStart, shr_p, shrE); //push_back(std::make_tuple(shrStart, shr_p, shrE));

	  n_shrs++;

	  float opAng = shr_openangle_v->at(i_pfp);
	  if( opAng > 1e-3 && opAng <= 1.57 ) h_ShrE_vs_OpAng->Fill( opAng, shrE );

	  //std::cout<<"Opening Angle V: "<<shr_openangle_v->at(i_pfp)<<
          //	       "\nShower (dx, dy, dz) = ("<<shr_distVec.x()<<", "<<shr_distVec.y()<<", "<<shr_distVec.z()<<
          //           "\nStart (x, y, z) = ("<<shrStart.x()<<", "<<shrStart.y()<<", "<<shrStart.z()<<
          //           "\nEnd (x, y, z) = ("<<shrEnd.x()<<", "<<shrEnd.y()<<", "<<shrEnd.z()<<")\n\n";
	}
      }
      if( numuCC == false ) continue;
      numEvents_numuCC++; nShrs_per_numuCC_event.push_back(n_shrs);

      //============================== PARTICLE LOOP ==============================//
      for(int i=0; i<all_mc_length; i++){

	int mc_trkID = all_mc_trkid->at(i); float mc_E = all_mc_E->at(i);
	int mc_pdg = all_mc_pdg->at(i); int mc_momTrkID = all_mc_mother->at(i);
	bool mc_isPrimary = mc_momTrkID == 0;

	map_allTrkID_partIDX[mc_trkID] = i;
	map_allTrkID_momTrkID[mc_trkID] = mc_momTrkID;

	// Printout this particle information to the screen 
	if( 1 && iEvent<10 ){   //&& part_isContained[i]
	  printf("%4i, TrackID: %-4i   PDG: %-12i E=%7.3f MeV,   moth=%5i, %18s \n",//   end(X,Y,Z)=(%1.2f, %1.2f, %1.2f)\n",   // dL=%6.1f cm, part_pathlen[i], other way to figure out dL with this file?
	    i, mc_trkID, mc_pdg, mc_E, mc_momTrkID, all_mc_process->at(i).c_str()//, all_mc_endx->at(i), all_mc_endy->at(i), all_mc_endz->at(i)
	  );
	}

        // -------------------------------------------------------------
        // Identify any neutrons and save their information
        if ( mc_pdg == 2112 ) {
          if( mc_isPrimary ) {   // primary neutrons
            total_nu_Ns++; h_trueNeutronE->Fill( (mc_E - 0.93957)*1000 );   // subtract neutron rest mass energy and convert from GeV to MeV
            map_n_p[mc_trkID].clear();
          }
        }

        if( !mc_isPrimary && map_n_p.count(mc_momTrkID) ) {   // from primary neutron mothers
          if( mc_pdg == 2212 ) {   // protons
	    map_n_p[mc_momTrkID].push_back(mc_trkID);
	    //if( std::find( trk_id.begin(), trk_id.end(),  mc_trkID ) != trk_id.end() ) h_prTrksPerE->Fill(mc_E);
            h_trueNeutronEvsProtonE->Fill( (all_mc_E->at(map_allTrkID_partIDX[mc_momTrkID]) - 0.93957)*1000, (all_mc_E->at(i) - 0.93827)*1000 );
	  }
        }

      }//end particle loop


      // cut booleans for ease of selection
      bool preselect = true;  bool nuProximity = true; bool vtxIso = false;
      bool energyCut = false; bool densityCut = false; bool BDTcut = true;
      //============================== BLIP LOOP ==============================//
      for(size_t i=0; i<nblips_saved; i++){

        totalBlips++; category=6;   // default to "other" if not assigned

	float blipX = blip_x->at(i), blipY = blip_y->at(i), blipZ = blip_z->at(i);
	if( minX > blipX ) minX=blipX; if( minY > blipY ) minY=blipY; if( minZ > blipZ ) minZ=blipZ;
	if( maxX < blipX ) maxX=blipX; if( maxY < blipY ) maxY=blipY; if( maxZ < blipZ ) maxZ=blipZ;
        TVector3 blip_loc( blipX, blipY, blipZ ); TVector3 reco_trace_vec( blip_loc - reco_nu_vtx );
	float dist = reco_trace_vec.Mag();	

	float blipE = blip_energy->at(i); float trueE = blip_true_energy->at(i);
	totalBlipE+=blipE; float E_res = abs(blipE - trueE)/trueE;

	TVector3 blip_true_dir( blip_true_dir_x->at(i), blip_true_dir_y->at(i), blip_true_dir_z->at(i) );

	float dx = blip_dx->at(i); float dw = blip_dw->at(i);
	float size = sqrt( dx*dx + dw*dw );   // direct multiplication is fastest ig

	float density = blipE/size; float dist2Trk = blip_proxtrkdist->at(i);
	float cutoff = cutoff_graph->Eval(blipE);

	int pdg = blip_true_pdg->at(i);	int g4id = blip_true_g4id->at(i);

	int allIDX = -9; bool isPrimary; 
	int mom_trkID = -9; bool mom_isPrimary = false;
	int mom_allIDX = -9; int mom_pdg = -9;
	if( map_allTrkID_partIDX.count(g4id) != 0 ) {
	  allIDX = map_allTrkID_partIDX[g4id];
	  mom_trkID = all_mc_mother->at( allIDX );
	  isPrimary = mom_trkID == 0;   // if mother == 0, particle is primary, could also use process="primary"
 
	  if( !isPrimary && map_allTrkID_partIDX.count(mom_trkID) != 0 ) {
	    mom_allIDX = map_allTrkID_partIDX[mom_trkID];
	    mom_isPrimary = all_mc_mother->at( mom_allIDX ) == 0; mom_pdg = all_mc_pdg->at( mom_allIDX );
	  }
	}

        if( !isPrimary && mom_pdg == 2112 && pdg == 2212 ) h_prBlipsPerE->Fill(trueE);

        //============================== CUTS ==============================//
	bool nearDeath = false; bool isProton = true; bool isInShower = false; bool pred_sig = true;
	if( preselect ) {
	  if( blip_touchtrk->at(i) || dist2Trk <= 1 ) continue;
	  nearDeath = ( blip_pl0_bydeadwire->at(i) == 1 || blip_pl1_bydeadwire->at(i) == 1 || blip_pl2_bydeadwire->at(i) == 1 );
	}
	if( nuProximity && dist >= 100 ) continue;
        if( pdg == 2212 && mom_pdg == 2112 && mom_isPrimary )   // signal blips
          { h_preFoundNeutronE->Fill( (all_mc_E->at(mom_allIDX) - 0.93957)*1000 ); }

	if( vtxIso && dist <= 5 ) continue;
	if( energyCut && blipE <= 3 ) continue;
	if( densityCut && density < Fit_Function_PID_Cut(blipE) ) isProton = false; //if( blipE <= 10 && density < cutoff ) isProton = false;
	for( const auto& [vtx, P, E] : shr_vtx_P_E ) { isInShower = IsWithinShowerCone(vtx, P, blip_loc, 90.0 ); }
	// To manually check if a blip is deemed to be inside any shower in this event:
	  //std::cout<<"\nShower vertex: ("<<vtx.X()<<", "<<vtx.Y()<<", "<<vtx.Z()<<")"
	  //<<"\nShower direction: ("<<P.X()<<", "<<P.Y()<<", "<<P.Z()<<")"; }
	//std::cout<<"Blip location: ("<<blip_loc.X()<<", "<<blip_loc.Y()<<", "<<blip_loc.Z()<<")";
	//std::cout<<"\nBlip is in shower: "<<std::boolalpha<<isInShower<<"\n\n";

	// BDT prediction; Compute output scores for binary classification and convert 
	blip_E=blipE; dist2Vtx=dist; blip_X=blipX; blip_Y=blipY; blip_Z=blipZ;
	blip_dX=dx; blip_dW=dw; blip_size=size; blip_density=density;
	//blip_dir_x=blip_rnn_dir_x->at(i); blip_dir_y=blip_rnn_dir_y->at(i); blip_dir_z=blip_rnn_dir_z->at(i);
	double tmva_score = reader->EvaluateMVA( "BDT" ); double base_score = 0.01855186;
	double XGBoost_score = 0.5 * TMath::Log( (1 + tmva_score)/(1 - tmva_score) ) - base_score;   // = atanh(tmva_score), sigmoid??
	if( BDTcut && XGBoost_score <= -0.40 ) pred_sig = false;   // BDT score CUT

        // RNN model prediction
        double res_angle = -1.0; TVector3 reco_dir; float align_angle, offset;
        if( blip_rnn_dir_isValid->at(i) ) {
          reco_dir.SetXYZ( blip_rnn_dir_x->at(i), blip_rnn_dir_y->at(i), blip_rnn_dir_z->at(i) );
          align_angle = reco_dir.Angle(reco_trace_vec); offset = abs(TMath::Cos(align_angle));   // pointing offset
          align_angle = align_angle*TMath::RadToDeg(); res_angle = reco_dir.Angle(blip_true_dir)*TMath::RadToDeg();
          if( align_angle > 90 ) align_angle = 180-align_angle;
          if( res_angle > 90 ) res_angle = 180-res_angle;

          //if( offset < 0.25 ) continue;   // hard direction cut

          if( !isInShower && !nearDeath && pred_sig && isProton ) {
            h_cosAlignment[0]->Fill(offset); h_alignmentAngle[0]->Fill(align_angle);
            if( g4id < 0 ) { h_cosAlignment[7]->Fill(offset); h_alignmentAngle[7]->Fill(align_angle); }
	    else if( abs(pdg) == 11 || abs(pdg) == 22 ) {
	      h_cosAlignment[5]->Fill(offset); h_alignmentAngle[5]->Fill(align_angle);
	      h_resolution[5]->Fill(res_angle); h_res_vs_E[5]->Fill(trueE, res_angle);
	    }
            else { h_resolution[0]->Fill(res_angle); h_res_vs_E[0]->Fill(trueE, res_angle); }
          }
        }

	// all selected blips
	if( !isInShower && isProton ) {
	  if( !nearDeath ) {
	    if( pred_sig ) {
	      totalSelectedBlips++;
	      h_recoVtxDist[0]->Fill(dist); h_blipE[0]->Fill(blipE);
	      h_Edensity[0]->Fill(density); h_energy_res[0]->Fill(E_res);
	      h_proxTrkDist[0]->Fill(dist2Trk); h_cosmicDistToVtx[0]->Fill(cosmicProx);
	      h_dx[0]->Fill(dx); h_dw[0]->Fill(dw);
	      h_size[0]->Fill(size);
	    } h_BDT_score[0]->Fill(XGBoost_score);
	  } h_byDeadWire[0]->Fill(nearDeath);
	} if( !nearDeath && isProton && pred_sig ) { h_inShower[0]->Fill(isInShower); }

	if( g4id < 0 ) {   // cosmic separation
          if( !isInShower ) {
	    if( !nearDeath ) {
	      if( pred_sig ) {
		if( blipE <= 10 ) { gr_cosmic_E_ds->SetPoint(gr_cosmic_E_ds->GetN(), blipE, density); }
		if( isProton ) {
		  h_recoVtxDist[7]->Fill(dist); h_blipE[7]->Fill(blipE);
		  h_Edensity[7]->Fill(density); h_energy_res[7]->Fill(E_res);
		  h_proxTrkDist[7]->Fill(dist2Trk); h_cosmicDistToVtx[7]->Fill(cosmicProx);
		  h_dx[7]->Fill(dx); h_dw[7]->Fill(dw); h_size[7]->Fill(size);
		}
	      } if( isProton ) h_BDT_score[7]->Fill(XGBoost_score);
	    } if( isProton ) h_byDeadWire[7]->Fill(nearDeath);
	  } if( !nearDeath && isProton && pred_sig ) h_inShower[7]->Fill(isInShower);
	  continue;
	}

	if( abs(pdg) != 11 && abs(pdg) != 22 ) { matchedBlips++;
	} else {   // photons, electrons, and positrons
          if( !isInShower ) {
	    if( !nearDeath ) {
	      if( pred_sig ) {
		if( blipE <= 10 ) gr_elec_E_ds->SetPoint(gr_elec_E_ds->GetN(), blipE, density);
		if( isProton ) {
		  h_recoVtxDist[5]->Fill(dist); h_blipE[5]->Fill(blipE);
		  h_Edensity[5]->Fill(density); h_energy_res[5]->Fill(E_res);
		  h_proxTrkDist[5]->Fill(dist2Trk); h_cosmicDistToVtx[5]->Fill(cosmicProx);
		  h_dx[5]->Fill(dx); h_dw[5]->Fill(dw); h_size[5]->Fill(size);
		}
	      } if( isProton ) h_BDT_score[5]->Fill(XGBoost_score);
	    } if( isProton ) h_byDeadWire[5]->Fill(nearDeath);
	  } if( !nearDeath && isProton && pred_sig ) { h_inShower[5]->Fill(isInShower); }
	  continue;
	}

	//if( blip_true_pdg->at(i) != all_mc_pdg->at(allIDX) ) std::cout<<"WOAH, HOLD UR HORSES! blip_true_pdg && all_mc_pdg don't match!!"<<endl;

	//h_nwires_pl0->Fill( blip_pl0_nwires->at(i) < 1 ? 0 : blip_pl0_nwires->at(i) );
	//h_nwires_pl1->Fill( blip_pl1_nwires->at(i) < 1 ? 0 : blip_pl1_nwires->at(i) ); 
	//h_nwires_pl2->Fill( blip_pl2_nwires->at(i) < 1 ? 0 : blip_pl2_nwires->at(i) ); 
	//int maxWires = std::max({blip_pl0_nwires->at(i), blip_pl1_nwires->at(i), blip_pl2_nwires->at(i)});
	//h_maxWires->Fill(maxWires);

        // PROTONS
	if( pdg == 2212 ) {
	  //std::cout<<"Proton mom PDG: "<<mom_pdg<<" Mom is primary: "<<mom_isPrimary<<endl;
	  if( isPrimary ) { category=4;   // primary protons
	  } else if( mom_pdg == 2112 ) {   // non-primary protons from neutron mom (should all be nIne)
	    if( mom_isPrimary ) {   // from *primary* neutron mom
	      if( blipE <= 10 && !nearDeath && !isInShower && pred_sig ) gr_signal_E_ds->SetPoint(gr_signal_E_ds->GetN(), blipE, density);
	      if( !nearDeath && !isInShower && isProton && pred_sig ) {
		signal_blips++; h_foundNeutronE->Fill( (all_mc_E->at(mom_allIDX) - 0.93957)*1000 );
		if( map_n_p[mom_trkID][0] == g4id ) neutronIDprBlips++;
	      }
	      if( map_n_p.at(mom_trkID).size() == 1 ) { category=1;   // 1p from primary neutron moms
	      } else { category=2; }   // Np from primary neutron moms
	    } else { category=3; }   // non-primary neutron mom
          }
	}

	//if( category == 4 ) std::cout<<"Primary proton energy resolution: "<<E_res<<std::endl;
	//if( category == 6 ) { std::cout<<"\nOH NAUR!, naur catagaury set ;( \nPDG = "<<p_pdg<<std::endl; }

	if( isProton ) {
	  if( !isInShower ) {
	    if( !nearDeath ) {
	      if( pred_sig ) {
		h_recoVtxDist[category]->Fill(dist); h_blipE[category]->Fill(blipE);
		h_Edensity[category]->Fill(density); h_energy_res[category]->Fill(E_res);
		h_proxTrkDist[category]->Fill(dist2Trk); h_cosmicDistToVtx[category]->Fill(cosmicProx);
		h_dx[category]->Fill(dx); h_dw[category]->Fill(dw); h_size[category]->Fill(size);
		if( blip_rnn_dir_isValid->at(i) ) {
		  h_cosAlignment[category]->Fill(offset); h_alignmentAngle[category]->Fill(align_angle);
		  h_resolution[category]->Fill(res_angle); h_res_vs_E[category]->Fill((all_mc_E->at(allIDX)-0.93827)*1000, res_angle);
		}
	      } h_BDT_score[category]->Fill(XGBoost_score);
	    } h_byDeadWire[category]->Fill(nearDeath);
	  } if( !nearDeath && pred_sig ) { h_inShower[category]->Fill(isInShower); }
	}

      }//end blip loop

      h_BlipsPerEvent->Fill(nblips_saved);

    }//end events loop

    float events_w_overOneShr = std::count_if(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end(), greater_than_1);
    float events_w_oneShr = std::count(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end(), 1);
    float events_w_noShr = std::count(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end(), 0);

    if( nblips_saved == 0) std::cout<<"No blips saved";

    int pre_bin_lim = h_preFoundNeutronE->FindFixBin(50);
    float pre_below50 = h_preFoundNeutronE->Integral(0,pre_bin_lim);
    int found_bin_lim = h_foundNeutronE->FindFixBin(50);
    float found_below50 = h_foundNeutronE->Integral(0,found_bin_lim);
    int neutronE_bin_lim = h_trueNeutronE->FindFixBin(50);
    float neutronE_below50 = h_trueNeutronE->Integral(0,neutronE_bin_lim);

    printf("\nBlip multiplicity: %8.3f +/- %.3f\n",
	   double(totalBlips)/double(selectedEvents), h_BlipsPerEvent->GetRMS(11) );

    std::cout<<"\nTotal Events: "<<numEvents<<" events\n"
             <<"Selected Events: "<<selectedEvents<<" events\n"
             <<"NumuCC Inclusive Selected Events: "<<numEvents_numuCC<<" events\n"
	     <<"Selected events without blips: "<<nonSelected<<" events\n"
             <<"Percent of selected events without reco nu vtx: "<<nonSelected*100/(selectedEvents+nonSelected)<<"%\n\n"

	     <<"Average number of showers per numuCC event: "<<TMath::Mean(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end())<<" showers\n"
	     <<"Percent of numuCC events with one shower: "<<100*events_w_oneShr/numEvents_numuCC<<"%\n"
	     <<"Percent of numuCC events with more than one shower: "<<100*events_w_overOneShr/numEvents_numuCC<<"%\n\n"
             <<"Percent of numuCC events without any showers: "<<100*events_w_noShr/numEvents_numuCC<<"%\n"

	     <<"Total blips: "<<totalBlips<<" blips\n"
	     <<"Total selected blips: "<<totalSelectedBlips<<" blips\n"
	     <<"Signal blips: "<<signal_blips<<" blips\n"
             <<"Percent of total blips with a true g4id ALSO not protons/electrons: "<<(matchedBlips*100)/totalBlips<<"%\n"

	     <<"Maximum blip (x, y, z): ("<<maxX<<", "<<maxY<<", "<<maxZ<<")\n"
             <<"Minimum blip (x, y, z): ("<<minX<<", "<<minY<<", "<<minZ<<")\n\n"

             <<"Mean dist to vtx: "<<h_recoVtxDist[0]->GetMean()<<"\n"
             <<"Mean reconstructed energy: "<<h_blipE[0]->GetMean()<<"\n\n";

    float eff = neutronIDprBlips/total_nu_Ns; float pur = signal_blips/totalSelectedBlips;
    printf( "Under 50 efficiency before selection: %5.3f%%, after selection: %5.3f%%\nEfficiency: %5.3f%% +/- %.2f\nPurity: %9.3f%% +/- %.2f\nE*P: %12.7f +/- %.2f\n\n",
	     pre_below50/neutronE_below50*100, found_below50/neutronE_below50*100,
	     eff*100, 1/sqrt(neutronIDprBlips),
	     pur*100, 1/sqrt(signal_blips),
	     eff*pur, 1/sqrt(neutronIDprBlips + signal_blips)
	  );
   	// could also do (max - min)/2 to estimate statistical uncertainty
  
    merge_plot2(h_prBlipsPerE, "Blips", h_prTrksPerE, "Tracks");

    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    h_ShrE_vs_OpAng->DrawCopy("hist");

    TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
    //gr_cosmic_E_ds->SetMarkerColor(kSpring-5);
    //gr_cosmic_E_ds->SetMarkerStyle(20);
    //gr_cosmic_E_ds->SetMarkerSize(0.4);
    //gr_cosmic_E_ds->Draw("AP");

    gr_elec_E_ds->SetMarkerColor(kViolet-5);
    gr_elec_E_ds->SetMarkerStyle(20);
    gr_elec_E_ds->SetMarkerSize(0.4);
    gr_elec_E_ds->Draw("AP");

    gr_signal_E_ds->SetMarkerColor(kOrange+6);
    gr_signal_E_ds->SetMarkerStyle(20);
    gr_signal_E_ds->SetMarkerSize(0.4);
    gr_signal_E_ds->Draw("P same");

    cutoff_graph->SetLineWidth(5);
    cutoff_graph->Draw("L same");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    //leg->AddEntry(gr_cosmic_E_ds, "Cosmics", "p");
    leg->AddEntry(gr_elec_E_ds, "Photons, e+, and e-", "p");
    leg->AddEntry(gr_signal_E_ds, "Signal Protons", "p");
    //leg->AddEntry(cutoff_graph, "90% of e- blips fall below this line", "L");
    leg->Draw();

    // Simple variable bin sizes above and below 250 MeV
    // 10 MeV bins below 250
    //for(double e = 0; e <= 250; e += 10) { binEdges.push_back(e); }
    // 50 MeV bins above 250
    //for(double e = 250; e <= 600; e += 50) { binEdges.push_back(e); } 

    // Define variable bin widths so that each bin contains 80 preselection only blips and 40 final selected blips
    // Need to refine this algorithm particulary so that the 80/40 numbers are statistically based, not arbitrary as they are now
    TCanvas *c5 = new TCanvas("c5", "Canvas", 800, 600);
    // define variable bin edges
    std::vector<double> binEdges;
    binEdges.push_back(h_trueNeutronE->GetBinLowEdge(1));
    double currentSum = 0;
    double upper_bound = 0;
    int zerosInnaRow = 0;
    for( int i = 1; i <= h_preFoundNeutronE->GetNbinsX(); i++ ) {
      if( zerosInnaRow > 30 ) {
        upper_bound = h_preFoundNeutronE->GetBinLowEdge(i-30);
        binEdges.push_back(upper_bound);
        break;
      }
      if( h_preFoundNeutronE->GetBinContent(i) == 0 ) { zerosInnaRow++; continue; }
      else zerosInnaRow = 0;
      currentSum += h_preFoundNeutronE->GetBinContent(i);
      if( currentSum > 80 ) {
        upper_bound = h_preFoundNeutronE->GetBinLowEdge(i+1); // Upper edge of current bin
        binEdges.push_back(upper_bound); currentSum = 0;
      }
    }
    h_preFoundNeutronE = (TH1D*)h_preFoundNeutronE->Rebin(binEdges.size()-1, "", &binEdges[0]);
    TH1D* h_pre_Eff_vs_E = (TH1D*)h_preFoundNeutronE->Clone("pre_Eff_vs_E");
    h_pre_Eff_vs_E->SetTitle("Selection Signal Efficiency vs. Energy;True Neutron Energy [MeV];Selection (n, Np) Efficiency");
    h_pre_Eff_vs_E->SetLineColor(kBlack); h_pre_Eff_vs_E->SetMarkerColor(kBlack);
    TH1D* h_preTrueNeutronE = (TH1D*)h_trueNeutronE->Rebin(binEdges.size()-1, "h_preTrueNeutronE", &binEdges[0]);
    h_pre_Eff_vs_E->Divide(h_preTrueNeutronE); //h_pre_Eff_vs_E->DrawCopy();
    TH1D* h_pre_drawn = (TH1D*)h_pre_Eff_vs_E->DrawCopy();
 
    binEdges.clear();
    binEdges.push_back(h_trueNeutronE->GetBinLowEdge(1));
    currentSum = 0; upper_bound = 0;zerosInnaRow = 0;
    for( int i = 1; i <= h_foundNeutronE->GetNbinsX(); i++ ) {
      if( zerosInnaRow > 30 ) {
        upper_bound = h_foundNeutronE->GetBinLowEdge(i-30);
        binEdges.push_back(upper_bound);
        break;
      }
      if( h_foundNeutronE->GetBinContent(i) == 0 ) { zerosInnaRow++; continue; }
      else zerosInnaRow = 0;
      currentSum += h_foundNeutronE->GetBinContent(i);
      if( currentSum > 40 ) {
        upper_bound = h_foundNeutronE->GetBinLowEdge(i+1); // Upper edge of current bin
        binEdges.push_back(upper_bound); currentSum = 0;
      }
    }
    h_foundNeutronE = (TH1D*)h_foundNeutronE->Rebin(binEdges.size()-1, "", &binEdges[0]);   
    TH1D* h_selection_Eff_vs_E = (TH1D*)h_foundNeutronE->Clone("selection_Eff_vs_E");
    h_selection_Eff_vs_E->SetTitle("Selection Signal Efficiency vs. Energy;True Neutron Energy [MeV];Selection (n, Np) Efficiency");
    h_selection_Eff_vs_E->SetLineColor(kRed); h_selection_Eff_vs_E->SetMarkerColor(kRed);
    TH1D* h_selTrueNeutronE = (TH1D*)h_trueNeutronE->Rebin(binEdges.size()-1, "h_selTrueNeutronE", &binEdges[0]);
    h_selection_Eff_vs_E->Divide(h_selTrueNeutronE); //h_selection_Eff_vs_E->DrawCopy("same");
    TH1D* h_sel_drawn = (TH1D*)h_selection_Eff_vs_E->DrawCopy("same");
    
    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_pre_drawn, "Preselection");
    leg1->AddEntry(h_sel_drawn, "After Cuts");
    leg1->Draw();

    TCanvas *c6 = new TCanvas("c6", "Canvas", 800, 600);
    TH1D* h_bkg_BDT_score = new TH1D("bkg_BDT_score","BDT Score for Background Blips;Score;Probability",100,0,1);
    TList* bkgList = new TList(); bkgList->Add(h_BDT_score[3]); bkgList->Add(h_BDT_score[4]);
    bkgList->Add(h_BDT_score[5]); bkgList->Add(h_BDT_score[6]); bkgList->Add(h_BDT_score[7]);
    h_bkg_BDT_score->Merge(bkgList);
    h_bkg_BDT_score->SetLineColor(kOrange); h_bkg_BDT_score->SetFillColorAlpha(kOrange, 0.5);
    h_bkg_BDT_score->Scale(1.0 / h_bkg_BDT_score->Integral("width"));
    h_bkg_BDT_score->DrawCopy("hist");

    TH1D* h_signal_BDT_score = new TH1D("signal_BDT_score","BDT Score for Signal and Background Blips;Score;Probability",100,0,1);
    TList* sigList = new TList(); sigList->Add(h_BDT_score[1]); sigList->Add(h_BDT_score[2]);
    h_signal_BDT_score->Merge(sigList); //h_signal_BDT_score
    h_signal_BDT_score->SetLineColor(kBlue); h_signal_BDT_score->SetFillColorAlpha(kBlue, 0.5);
    h_signal_BDT_score->Scale(1.0 / h_signal_BDT_score->Integral("width"));
    h_signal_BDT_score->DrawCopy("hist same");

    TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg2->AddEntry(h_signal_BDT_score, "Signal");
    leg2->AddEntry(h_bkg_BDT_score, "Background");
    leg2->Draw();

    TCanvas *c7 = new TCanvas("c7", "Canvas", 800, 600);
    h_trueNeutronEvsProtonE->DrawCopy();

    Float_t a;
    TTree otree("evd_tree", "TTree with a single float branch for total evd");
    TBranch* branch = otree.Branch("total_event_displays", &a, "total_event_displays/F");
    a = selectedEvents;
    otree.Fill();
  
    //===========================================
    // Write the output and close the file
    outfile->Write();
    outfile->Close();

  }
