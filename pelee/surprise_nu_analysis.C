
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
  int maxEvts = -2000;
  int maxRuns = -2000;

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
  TDirectory* cos_theta_plots;
  //TDirectory* trueCosTheta_vs_momTrueE_plots;

  TDirectory* neutronEnergy_plots;
  TDirectory* position_plots;
  TDirectory* theta_plots;
  TDirectory* delta_pT_plots;
  TDirectory* shower_plots;
  TDirectory* secondaries_plots;


  int			   n_pfps = 0;
  //int			   selected;
  //int			   n_tracks;
  //int			   n_showers;
  int			   nu_pdg = 0;
  int			   ccnc = 0;
  float			   reco_nu_vtx_sce_x = 0;
  float			   reco_nu_vtx_sce_y = 0;
  float			   reco_nu_vtx_sce_z = 0;
  int			   slice_orig_pass_id = 0;
  //bool			   isVtxInFiducial = 0;
  //float 		   topological_score = 0;
  float			   slice_orig_topo_score = 0;
  int			   crtveto = 0;
  double		   _closestNuCosmicDist = 0;
  //float			   bdt_numuCC = 0;   // doesn't exist, all values 9999
  int			   nblips_saved = 0;
  float			   pt = 0;
  int			   run = 0;
  int			   runTreeRun = 0;
  float			   pot = 0;
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
  std::vector<float>       * blip_trkid = 0;
  std::vector<int>	   * blip_true_primancPDG = 0; // primary blip ancestor PDG
  std::vector<int>	   * blip_true_primancG4ID = 0; // primary ancestor G4 TrackID
  std::vector<bool>	   * blip_rnn_dir_isValid = 0; // Is RNN output valid?
  std::vector<float>	   * blip_rnn_dir_x = 0;       // RNN-predicted blip direction (protons)
  std::vector<float>	   * blip_rnn_dir_y = 0;       // RNN-predicted blip direction (protons)
  std::vector<float>	   * blip_rnn_dir_z = 0;       // RNN-predicted blip direction (protons)
  std::vector<float>	   * blip_true_dir_x = 0;   // Initial momentum direction of truth-matched particle
  std::vector<float>	   * blip_true_dir_y = 0;   // Initial momentum direction of truth-matched particle
  std::vector<float>	   * blip_true_dir_z = 0;   // Initial momentum direction of truth-matched particle
  // the above three^ are the unit vector components (each ranges from -1 to 1)
  std::vector<int>	   * blip_true_ncategory = 0;   // Help categorize origin of blip
	//  -9 = no truth match (data/overlay)
	//  0  = truth-matched, but not falling in category
	//  1  = primary (n,1p)
	//  2  = primary (n,Np)
	//  3  = secondary (n,1p)
	//  4  = secondary (n,Np)
	//  5  = primary (n,gamma)
	//  6  = secondary (n,gamma)
	//  7  = ncapture gamma
	//  8  = mu capture gamma
  std::vector<int>	   * pfng2mipfrac = 0;
  //std::vector<int>	   * pfng2semlabel = 0;
  //std::vector<int>         * trksemlbl = 0;

  std::vector<float>	   * trk_sce_start_x_v = 0;
  std::vector<float>	   * trk_sce_start_y_v = 0;
  std::vector<float>	   * trk_sce_start_z_v = 0;
  std::vector<float>	   * trk_sce_end_x_v = 0;
  std::vector<float>	   * trk_sce_end_y_v = 0;
  std::vector<float>	   * trk_sce_end_z_v = 0;
  std::vector<float>	   * trk_llr_pid_score_v = 0;
  std::vector<float>       * trk_score_v = 0;
  std::vector<float>       * trk_pfp_id_v = 0;
  std::vector<float>	   * trk_mcs_muon_mom_v = 0;
  std::vector<float>	   * trk_energy_proton_v = 0;

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
    // Shower start distance from reco neutrino vertex (pre-calculated for
    // convenience)
    //set_object_input_branch_address( etree, "shr_dist_v",

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

  int                        n_pfps_data = 0;
  std::vector<float>       * data_trk_score_v = 0;
  std::vector<float>       * data_shr_start_x_v = 0;
  std::vector<float>       * data_shr_start_y_v = 0;
  std::vector<float>       * data_shr_start_z_v = 0;
  std::vector<float>       * data_shr_px_v = 0;
  std::vector<float>       * data_shr_py_v = 0;
  std::vector<float>       * data_shr_pz_v = 0;
  std::vector<float>       * blip_energy_data = 0;
  std::vector<float>       * blip_x_data = 0;
  std::vector<float>       * blip_y_data = 0;
  std::vector<float>       * blip_z_data = 0;
  int                        nblips_saved_data = 0;


  const int numHists=11;
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
  TH1D* h_cos_theta[numHists];
  //TH2D* h_trueCosTheta_vs_momTrueE[numHists];

  TH2D* h_trueNeutronEvsProtonE;
  TH1D* h_trueNeutronE;
  TH1D* h_preFoundNeutronE;
  TH1D* h_foundNeutronE;
  TH2D* h_trueNeutronEvsRecoAngle;

  TH2D* h_signal_ZY;
  TH2D* h_signal_ZX;
  TH1D* h_signal_X;
  TH1D* h_signal_Y;
  TH1D* h_signal_Z;
  TH1D* h_trueSignalX;
  TH1D* h_trueSignalY;
  TH1D* h_trueSignalZ;
  //TH1D* h_bdt_numuCC;

  TH1D* h_trueSig_theta;
  TH1D* h_trueSig_cosTheta;
  TH2D* h_trueNeutronEvsCosTheta;
  TH1D* h_ProtonCosTheta;
  TH1D* h_NeutronCosTheta;
  TH1D* h_PionCosTheta;
  TH1D* h_OtherCosTheta;

  TH1D* h_delta_pT[4];

  //TH2D* h_ShrE_vs_OpAng;
  TH1D* h_shrEnergy;
  TH1D* h_nBlipsInShr;
  TH1D* h_shrBlipE_data;
  TH1D* h_shrBlipE_MC;

  TH1D* h_secondPr_cosTheta;
  //TH1D* h_trueCosTheta_secondPr;

  TH1D* h_BlipsPerEvent;
  // these are also per event
  TH1D* h_nBlips_xp0n;
  TH1D* h_nBlips_xpxn;
  TH1D* h_nBlips_0p0n;
  TH1D* h_nBlips_0pxn;
  TH1D* h_maxWires;
  //TH1D* h_nwires_pl0;
  //TH1D* h_nwires_pl1;
  //TH1D* h_nwires_pl2;
  TH1D* h_prTrksPerE;
  TH1D* h_prBlipsPerE;
  TH1D* h_selectedsMomE;

  TGraph* gr_signal_E_ds;
  TGraph* gr_elec_E_ds;
  TGraph* gr_cosmic_E_ds;

  // Input file
  //std::string fFileName = "checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root";   // retupled BNB nu
  std::string fOverlayFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root";
  //std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260224.root";   // BNB nu with more blip vars (full events)
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";   // BNB nu (no all_mc)

  std::string fDataFileName = "MCC9p10_Run4b_v10_04_07_11_BNB_opendata_surprise_blipPlus_20260522.root";

  std::string nuTreeName = "nuselection/NeutrinoSelectionFilter";
  std::string subRunTreeName = "nuselection/SubRun";
  //std::string wcTreeName = "wcpselection/T_PFeval";

  std::unordered_map<int, int> map_allTrkID_partIDX;
  std::unordered_map<int, int> map_allTrkID_momTrkID;

  //============================== FUNCTIONS ==============================//
  // Function to determine if a particle descended from another particle.
  // Allows option to break lineage at photons for contiguous parentage.
  bool IsAncestorOf(int particleID, int ancestorID, bool breakAtPhots = false){
    if( particleID == ancestorID  )       return true;
    if( particleID < ancestorID   )       return false;
    if( map_allTrkID_momTrkID.count(ancestorID) == 0 )  return false;   // may silently fail for valid primaries that have trkID = 0

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

  bool IsInFV_Exact( const TVector3& point ) {
    return point.X() >= 0.      &&  point.X() <= 256.	&&
           point.Y() >= -116.25 &&  point.Y() <= 116.25	&&
           point.Z() >= 0.      &&  point.Z() <= 1036.8	;
  }   // from R. Acciarri et al 2017 JINST 12 P02017

  bool IsInFV_Buffer( const TVector3& point ) {
    return point.X() >= 10.	&&  point.X() <= 246.35	&&
           point.Y() >= -106.5	&&  point.Y() <= 106.5	&&
           point.Z() >= 10.	&&  point.Z() <= 1026.8	;
  }   // Same as Burk and Andy's analysis but may want to loosen since our analysis isn't track focused

  bool IsInFV_BufferAndDeadZone( const TVector3& point ) {
    return point.X() >= 5.	&&  point.X() <= 251.	&&
	   point.Y() >= -110.	&&  point.Y() <= 110	&&
	   point.Z() >= 20.	&&  point.Z() <= 986.	&&
	   (point.Z() < 675.	||  point.Z() > 775.)	;
  }   // 5cm buffer and dead wire removal used for numuCC nu vtx selection

  bool IsWithinShowerCone(const TVector3& vtx,	    // Shower Vtx
			  const TVector3& dir,	    // Shower direction (momentum)
			  const TVector3& blip,	    // Blip coordinates
			  double openingAngleDeg) { // Total opening angle of exclusion cone 
    if( vtx.X() < -100 ) return false;
    TVector3 d = blip - vtx;
    double angle = d.Angle(dir);
    double opang = openingAngleDeg * (TMath::Pi() / 180.) * 0.5;
    return ( angle < opang );
  }

  // I think the problem with this is that it needs to be inside the box, but IsInFV ?!?
  double distToNearestEdge( TVector3 p ){
    return std::min({ p.X() - 0.,      256.  - p.X(),
                      p.Y() + 116.25,  116.25 - p.Y(),
                      p.Z() - 0.,      1036.8 - p.Z() });
  }

  //double origin[3] = {0, 116.25, 0};
  //const TGeoBBox *MicroBox = new TGeoBBox( 256.0, 232.5, 1036.8, origin );   // from R. Acciarri et al 2017 JINST 12 P02017
  //double distToNearestEdge( TVector3 location ){
  //  const double loc[3] = {location.X(), location.Y(), location.Z()};
  //  return MicroBox->Safety( loc );//DistFromInside( location );
  //}

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

  // Counting number of blips in each energy bin to calculate efficiency per energy range
  std::tuple<float, float, int> effPurBelowEnergy( TH1D* trueSelectedEnergyHist, TH1D* totalEnergyHist, TH1D* allSelectedEnergyHist, int neutronEnergyLimit ) {
    int bin_lim = trueSelectedEnergyHist->FindFixBin(neutronEnergyLimit);
    float nTruSelectedBelowE = trueSelectedEnergyHist->Integral(0,bin_lim);
    float nTotalBelowE = totalEnergyHist->Integral(0,bin_lim);
    float nAllSelectedBelowE = allSelectedEnergyHist->Integral(0,bin_lim);
    float eff = nTruSelectedBelowE/nTotalBelowE*100;
    float pur = nTruSelectedBelowE/nAllSelectedBelowE*100;
    return {eff, pur, nTruSelectedBelowE};
  }

  std::tuple<float, float, int> effPurAboveEnergy( TH1D* trueSelectedEnergyHist, TH1D* totalEnergyHist, TH1D* allSelectedEnergyHist, int neutronEnergyMin ) {
    int bin_min = trueSelectedEnergyHist->FindFixBin(neutronEnergyMin);
    float nTruSelectedAboveE = trueSelectedEnergyHist->Integral( bin_min, trueSelectedEnergyHist->GetNbinsX() );
    float nTotalAboveE = totalEnergyHist->Integral( bin_min, totalEnergyHist->GetNbinsX() );
    float nAllSelectedAboveE = allSelectedEnergyHist->Integral( bin_min, totalEnergyHist->GetNbinsX() );
    float eff = nTruSelectedAboveE/nTotalAboveE*100;
    float pur = nTruSelectedAboveE/nAllSelectedAboveE*100;
    return {eff, pur, nTruSelectedAboveE};
  }

  //float effOrPurAboveBelowEnergy( TH1D* selectedEnergyHist, TH1D* totalEnergyHist, TH1D* trueSigEnergyHist, int neutronEnergyMin, int neutronEnergyLimit ) {
  //  int bin_min = selectedEnergyHist->FindFixBin(neutronEnergyMin);
  //  int bin_lim = selectedEnergyHist->FindFixBin(neutronEnergyLimit);
  //  float nSelectedAboveBelowE = selectedEnergyHist->Integral( bin_min, bin_lim );
  //  float nTotalAboveBelowE = totalEnergyHist->Integral( bin_min, bin_lim );
  //  float eff = nSelectedAboveBelowE/nTotalAboveBelowE*100;
  //  return eff;
  //}

  // Normalize each histogram to unit area (for A.U. y-axis), use Integral("width") if bin width > 1 / you want true prob. density
  auto Normalize = [](TH1D* h) { double integral = h->Integral("width"); if (integral > 0) h->Scale(1.0 / integral); };

  //============================== LOAD BDT ==============================//
  //TMVA::Experimental::RReader model("n_p_select_noShow_100k.xml");
  TMVA::Reader *reader = new TMVA::Reader( );

  // Provide the features exactly as used during training
  //std::vector<std::string> featureNames = {"blip_E", "dist2Vtx", "blip_X", "blip_Y", "blip_Z", "blip_dX", "blip_dW", "blip_size", "blip_density"};
  float blip_E, dist2Vtx, dist2Trk, dist2Edge, blip_size, blip_density, alignment_angle, Z_dist2Vtx;//, dist2ShrVtx, angle_blip_shr;
  reader->AddVariable("blip_E", &blip_E);
  reader->AddVariable("dist2Vtx", &dist2Vtx);
  reader->AddVariable("dist2Trk", &dist2Trk);
  reader->AddVariable("dist2Edge", &dist2Edge);
  //reader->AddVariable("blip_X", &blip_X);
  //reader->AddVariable("blip_Y", &blip_Y);
  //reader->AddVariable("blip_Z", &blip_Z);
  //reader->AddVariable("blip_dX", &blip_dX);
  //reader->AddVariable("blip_dW", &blip_dW);
  reader->AddVariable("blip_size", &blip_size);
  reader->AddVariable("blip_density", &blip_density);
  reader->AddVariable("alignment_angle", &alignment_angle);
  reader->AddVariable("Z_dist2Vtx", &Z_dist2Vtx);
  //reader->AddVariable("dist2ShrVtx", &dist2ShrVtx);
  //reader->AddVariable("angle_blip_shr", &angle_blip_shr);

  reader->BookMVA("BDT", "models_BDT/n_p_select_truNuNoShr.xml");   // MAKE SURE TO CHANGE WHICH EVENTS IF USING TRUE NUMUCC
  reader->BookMVA("highE_BDT", "models_BDT/n_p_select_highE.xml");
  reader->BookMVA("lowE_BDT", "models_BDT/n_p_select_lowE.xml");

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

    // Read in the MC and Data TFiles and grab the TTrees
    TFile* overlayFile = new TFile(("files/"+fOverlayFileName).c_str(),"READ");
    TTree* overlayNuTree = (TTree*)overlayFile->Get(nuTreeName.c_str());
    TTree* subRunTree = (TTree*)overlayFile->Get(subRunTreeName.c_str());
    //TTree* wcTree = (TTree*)file->Get(wcTreeName.c_str());
    std::cout<<"\nReading data file "<<fOverlayFileName<<"\n";
    TFile* dataFile = new TFile(("files/"+fDataFileName).c_str(),"READ");
    TTree* dataNuTree = (TTree*)dataFile->Get(nuTreeName.c_str());
    std::cout<<"Reading data file "<<fDataFileName<<"\n";


    // Set branch address for all variables you plan to use
    //overlayNuTree->SetBranchAddress("selected",	    &selected);
    //overlayNuTree->SetBranchAddress("n_tracks",	    &n_tracks);
    overlayNuTree->SetBranchAddress("nblips_saved",	    &nblips_saved);
    overlayNuTree->SetBranchAddress("blip_true_g4id",	    &blip_true_g4id);
    overlayNuTree->SetBranchAddress("blip_x",		    &blip_x);
    overlayNuTree->SetBranchAddress("blip_y",		    &blip_y);
    overlayNuTree->SetBranchAddress("blip_z",		    &blip_z);
    overlayNuTree->SetBranchAddress("blip_dx",             &blip_dx);
    overlayNuTree->SetBranchAddress("blip_dw",             &blip_dw);
    overlayNuTree->SetBranchAddress("blip_energy",	    &blip_energy);
    //overlayNuTree->SetBranchAddress("blip_charge",         &blip_charge);
    //overlayNuTree->SetBranchAddress("blip_nplanes",        &blip_nplanes);
    overlayNuTree->SetBranchAddress("blip_proxtrkdist",    &blip_proxtrkdist);
    overlayNuTree->SetBranchAddress("blip_touchtrk",       &blip_touchtrk);
    overlayNuTree->SetBranchAddress("blip_true_pdg",       &blip_true_pdg);
    overlayNuTree->SetBranchAddress("blip_true_energy",    &blip_true_energy);

    //overlayNuTree->SetBranchAddress("blip_pl0_nwires",	    &blip_pl0_nwires);
    //overlayNuTree->SetBranchAddress("blip_pl1_nwires",	    &blip_pl1_nwires);
    //overlayNuTree->SetBranchAddress("blip_pl2_nwires",	    &blip_pl2_nwires);
    overlayNuTree->SetBranchAddress("blip_pl0_bydeadwire", &blip_pl0_bydeadwire);
    overlayNuTree->SetBranchAddress("blip_pl1_bydeadwire", &blip_pl1_bydeadwire);
    overlayNuTree->SetBranchAddress("blip_pl2_bydeadwire", &blip_pl2_bydeadwire);
    overlayNuTree->SetBranchAddress("blip_trkid",	   &blip_trkid);
    overlayNuTree->SetBranchAddress("blip_true_primancPDG",&blip_true_primancPDG); 
    //overlayNuTree->SetBranchAddress("blip_true_primancG4ID",&blip_true_primancG4ID);
    overlayNuTree->SetBranchAddress("blip_true_ncategory", &blip_true_ncategory);
    overlayNuTree->SetBranchAddress("blip_rnn_dir_isValid",&blip_rnn_dir_isValid);
    overlayNuTree->SetBranchAddress("blip_rnn_dir_x",	    &blip_rnn_dir_x);
    overlayNuTree->SetBranchAddress("blip_rnn_dir_y",	    &blip_rnn_dir_y);
    overlayNuTree->SetBranchAddress("blip_rnn_dir_z",	    &blip_rnn_dir_z);
    overlayNuTree->SetBranchAddress("blip_true_dir_x",	    &blip_true_dir_x);
    overlayNuTree->SetBranchAddress("blip_true_dir_y",	    &blip_true_dir_y);
    overlayNuTree->SetBranchAddress("blip_true_dir_z",	    &blip_true_dir_z);

    overlayNuTree->SetBranchAddress("n_pfps",		    &n_pfps);
    //overlayNuTree->SetBranchAddress("pfng2mipfrac",	    &pfng2mipfrac);
    //overlayNuTree->SetBranchAddress("pfng2semlabel",	    &pfng2semlabel);
    //overlayNuTree->SetBranchAddress("trksemlbl",	    &trksemlbl);   // doesn't seem to exist...

    overlayNuTree->SetBranchAddress("trk_sce_start_x_v",   &trk_sce_start_x_v);
    overlayNuTree->SetBranchAddress("trk_sce_start_y_v",   &trk_sce_start_y_v);
    overlayNuTree->SetBranchAddress("trk_sce_start_z_v",   &trk_sce_start_z_v);
    overlayNuTree->SetBranchAddress("trk_sce_end_x_v",     &trk_sce_end_x_v);
    overlayNuTree->SetBranchAddress("trk_sce_end_y_v",     &trk_sce_end_y_v);
    overlayNuTree->SetBranchAddress("trk_sce_end_z_v",     &trk_sce_end_z_v);
    overlayNuTree->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v);
    overlayNuTree->SetBranchAddress("trk_score_v",	    &trk_score_v);
    overlayNuTree->SetBranchAddress("trk_pfp_id_v",	    &trk_pfp_id_v);

    //overlayNuTree->SetBranchAddress("n_showers",           &n_showers);
    overlayNuTree->SetBranchAddress("shr_start_x_v",	    &shr_start_x_v);
    overlayNuTree->SetBranchAddress("shr_start_y_v",	    &shr_start_y_v);
    overlayNuTree->SetBranchAddress("shr_start_z_v",	    &shr_start_z_v);
    overlayNuTree->SetBranchAddress("shr_px_v",	    &shr_px_v);
    overlayNuTree->SetBranchAddress("shr_py_v",	    &shr_py_v);
    overlayNuTree->SetBranchAddress("shr_pz_v",	    &shr_pz_v);
    //overlayNuTree->SetBranchAddress("shr_dist_v",	    &shr_dist_v);
    //overlayNuTree->SetBranchAddress("shr_openangle_v",     &shr_openangle_v);
    //overlayNuTree->SetBranchAddress("shr_openangle",       &shr_openangle);
    overlayNuTree->SetBranchAddress("shr_energy_y_v",	    &shr_energy_y_v);

    overlayNuTree->SetBranchAddress("all_mc_trkid",	    &all_mc_trkid);
    overlayNuTree->SetBranchAddress("all_mc_pdg",	    &all_mc_pdg);
    overlayNuTree->SetBranchAddress("all_mc_E",	    &all_mc_E);
    overlayNuTree->SetBranchAddress("all_mc_endx",         &all_mc_endx);
    overlayNuTree->SetBranchAddress("all_mc_endy",         &all_mc_endy);
    overlayNuTree->SetBranchAddress("all_mc_endz",         &all_mc_endz);
    overlayNuTree->SetBranchAddress("all_mc_mother",       &all_mc_mother);
    //overlayNuTree->SetBranchAddress("all_mc_process",      &all_mc_process);

    overlayNuTree->SetBranchAddress("nu_pdg",		    &nu_pdg);
    overlayNuTree->SetBranchAddress("ccnc",		    &ccnc);
    overlayNuTree->SetBranchAddress("reco_nu_vtx_sce_x",   &reco_nu_vtx_sce_x);
    overlayNuTree->SetBranchAddress("reco_nu_vtx_sce_y",   &reco_nu_vtx_sce_y);
    overlayNuTree->SetBranchAddress("reco_nu_vtx_sce_z",   &reco_nu_vtx_sce_z);
    overlayNuTree->SetBranchAddress("slice_orig_pass_id",  &slice_orig_pass_id);
    //overlayNuTree->SetBranchAddress("isVtxInFiducial",     &isVtxInFiducial);   // this is a true var
    //overlayNuTree->SetBranchAddress("topological_score",   &topological_score);
    overlayNuTree->SetBranchAddress("slice_orig_topo_score",&slice_orig_topo_score);
    overlayNuTree->SetBranchAddress("crtveto",		    &crtveto);
    overlayNuTree->SetBranchAddress("_closestNuCosmicDist", &_closestNuCosmicDist);
    overlayNuTree->SetBranchAddress("pt",		    &pt);
    //overlayNuTree->SetBranchAddress("reco_nu_vtx_x",       &reco_nu_vtx_x);
    //overlayNuTree->SetBranchAddress("reco_nu_vtx_y",       &reco_nu_vtx_y);
    //overlayNuTree->SetBranchAddress("reco_nu_vtx_z",       &reco_nu_vtx_z);
    overlayNuTree->SetBranchAddress("run",		    &run);

    dataNuTree->SetBranchAddress("n_pfps",                  &n_pfps_data);
    dataNuTree->SetBranchAddress("trk_score_v",                  &data_trk_score_v);
    dataNuTree->SetBranchAddress("shr_start_x_v",        &data_shr_start_x_v);
    dataNuTree->SetBranchAddress("shr_start_y_v",        &data_shr_start_y_v);
    dataNuTree->SetBranchAddress("shr_start_z_v",        &data_shr_start_z_v);
    dataNuTree->SetBranchAddress("shr_px_v",             &data_shr_px_v);
    dataNuTree->SetBranchAddress("shr_py_v",             &data_shr_py_v);
    dataNuTree->SetBranchAddress("shr_pz_v",             &data_shr_pz_v);
    dataNuTree->SetBranchAddress("blip_energy",          &blip_energy_data);
    dataNuTree->SetBranchAddress("blip_x",               &blip_x_data);
    dataNuTree->SetBranchAddress("blip_y",                  &blip_y_data);
    dataNuTree->SetBranchAddress("blip_z",                  &blip_z_data);
    dataNuTree->SetBranchAddress("nblips_saved",                  &nblips_saved_data);


    subRunTree->SetBranchAddress("run",		    &runTreeRun);
    //subRunTree->SetBranchAddress("subRun",          &subRun);
    subRunTree->SetBranchAddress("pot",		    &pot);

    // make output file to before declaring plots
    TFile* outfile = new TFile("test_out.root", "recreate");
    blipE_plots			= outfile->mkdir("blipE_plots");
    recoVtxDist_plots     	= outfile->mkdir("recoVtxDist_plots");
    cosAlignment_plots		= outfile->mkdir("cosAlignment_plots");
    alignmentAngle_plots	= outfile->mkdir("alignmentAngle_plots");
    resolution_plots		= outfile->mkdir("resolution_plots");
    res_vs_E_plots		= outfile->mkdir("res_vs_E_plots");
    Edensity_plots		= outfile->mkdir("Edensity_plots");
    byDeadWire_plots		= outfile->mkdir("byDeadWire_plots");
    proxTrkDist_plots		= outfile->mkdir("proxTrkDist_plots");
    cosmicDistToVtx_plots	= outfile->mkdir("cosmicDistToVtx_plots");
    dx_plots			= outfile->mkdir("dx_plots");
    dw_plots			= outfile->mkdir("dw_plots");
    size_plots			= outfile->mkdir("size_plots");
    inShower_plots		= outfile->mkdir("inShower_plots");
    energy_res_plots		= outfile->mkdir("energy_res_plots");
    BDT_score_plots		= outfile->mkdir("BDT_score_plots");
    cos_theta_plots		= outfile->mkdir("cos_theta_plots");
    //trueCosTheta_vs_momTrueE_plots = outfile->mkdir("trueCosTheta_vs_momTrueE");

    neutronEnergy_plots		= outfile->mkdir("neutronEnergy_plots");
    position_plots		= outfile->mkdir("position_plots");
    theta_plots			= outfile->mkdir("theta_plots");
    delta_pT_plots		= outfile->mkdir("delta_pT_plots");
    shower_plots		= outfile->mkdir("shower_plots");
    secondaries_plots		= outfile->mkdir("secondaries_plots");

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
    h_resolution[0]         = new TH1D("resolution","Blip Angular Resolution;#Theta_{pRT} [degrees]",30,0,90);//46,-1,91);
    for(int i=1; i<numHists; i++) { h_resolution[i] = (TH1D*)h_resolution[0]->Clone(("resolution_" + std::to_string(i)).c_str()); }
    Edensity_plots->cd();
    h_Edensity[0]           = new TH1D("Edensity","Blip Energy Deposition Density;E/ds [MeVee/cm]",20,0,10);
    for(int i=1; i<numHists; i++) { h_Edensity[i] = (TH1D*)h_Edensity[0]->Clone(("Edensity_" + std::to_string(i)).c_str()); }
    byDeadWire_plots->cd();
    h_byDeadWire[0]         = new TH1D("byDeadWire","Blip By Dead Wire",2,0,2);
    for(int i=1; i<numHists; i++) { h_byDeadWire[i] = (TH1D*)h_byDeadWire[0]->Clone(("byDeadWire" + std::to_string(i)).c_str()); }
    proxTrkDist_plots->cd();
    h_proxTrkDist[0]        = new TH1D("proxTrkDist","Blip Distance to Nearest Track;Reconstructed Distance [cm]",21,0,105);
    for(int i=1; i<numHists; i++) { h_proxTrkDist[i] = (TH1D*)h_proxTrkDist[0]->Clone(("proxTrkDist_" + std::to_string(i)).c_str()); }
    cosmicDistToVtx_plots->cd();
    h_cosmicDistToVtx[0]    = new TH1D("cosmicDistToVtx","Reco. Neutrino Vertex Distance to Nearest Cosmic Track;Reconstructed Distance [cm]",41,0,510);
    for(int i=1; i<numHists; i++) { h_cosmicDistToVtx[i] = (TH1D*)h_cosmicDistToVtx[0]->Clone(("cosmicDistToVtx_" + std::to_string(i)).c_str()); }
    dx_plots->cd();
    h_dx[0]		    = new TH1D("dx","Blip dX;dX [cm]",24,0,6);
    for(int i=1; i<numHists; i++) { h_dx[i] = (TH1D*)h_dx[0]->Clone(("dx_" + std::to_string(i)).c_str()); }
    dw_plots->cd();
    h_dw[0]		    = new TH1D("dw","Blip dW;dW [cm]",20,0,5);
    for(int i=1; i<numHists; i++) { h_dw[i] = (TH1D*)h_dw[0]->Clone(("dw_" + std::to_string(i)).c_str()); }
    size_plots->cd();
    h_size[0]		    = new TH1D("size","Blip Size;3D Projected Size [cm]",28,0,7);
    for(int i=1; i<numHists; i++) { h_size[i] = (TH1D*)h_size[0]->Clone(("size_" + std::to_string(i)).c_str()); }
    inShower_plots->cd();
    h_inShower[0]           = new TH1D("inShower","Blip Inside Shower;In or Out of Shower;nBlips",2,0,2);
    for(int i=1; i<numHists; i++) { h_inShower[i] = (TH1D*)h_inShower[0]->Clone(("inShower_" + std::to_string(i)).c_str()); }
    energy_res_plots->cd();
    h_energy_res[0]	    = new TH1D("energy_res","Blip Energy Resolution;(recoE-trueE)/trueE;nBlips",20,0,4);
    for(int i=1; i<numHists; i++) { h_energy_res[i] = (TH1D*)h_energy_res[0]->Clone(("energy_res_" + std::to_string(i)).c_str()); }
    BDT_score_plots->cd();
    h_BDT_score[0]	    = new TH1D("BDT_score","Blip BDT Prediction;BDT Score (1 is signal-like);nBlips",70,-12,2);
    for(int i=1; i<numHists; i++) { h_BDT_score[i] = (TH1D*)h_BDT_score[0]->Clone(("BDT_score_" + std::to_string(i)).c_str()); }
    cos_theta_plots->cd();
    h_cos_theta[0]	    = new TH1D("cos_theta","Reconstructed Blip Cos(Theta Angle);Cos(#Theta_{n});nBlips",20,-1,1);
    for(int i=1; i<numHists; i++) { h_cos_theta[i] = (TH1D*)h_cos_theta[0]->Clone(("cos_theta_" + std::to_string(i)).c_str()); }
    //trueCosTheta_vs_momTrueE_plots->cd();
    //h_trueCosTheta_vs_momTrueE[0]= new TH2D("trueCosTheta_vs_momTrueE","True Cos(Theta) vs True Energy of Mother Particle;Mother True Kinetic Energy [MeV];True Cos(#Theta)",18,20,200,10,-1,1);
    //for(int i=1; i<numHists; i++) { h_trueCosTheta_vs_momTrueE[i] = (TH2D*)h_trueCosTheta_vs_momTrueE[0]->Clone(("trueCosTheta_vs_momTrueE_" + std::to_string(i)).c_str()); }

    //h_mostWiresPlane[]
    //  std::vector<int>         * blip_pl0_nwires = 0;
    //  std::vector<int>         * blip_pl1_nwires = 0;
    //  std::vector<int>         * blip_pl2_nwires = 0;

    neutronEnergy_plots->cd();
    h_trueNeutronE		= new TH1D("trueNeutronE","True Neutron Energy;MeVee",600,0,600); h_trueNeutronE->Sumw2();
    h_preFoundNeutronE		= new TH1D("preFoundNeutronE","Preselection Identified Neutron True Kinetic Energy;MeV",600,0,600); h_preFoundNeutronE->Sumw2();
    h_foundNeutronE		= new TH1D("foundNeutronE","Identified Neutron True Kinetic Energy;MeV",600,0,600); h_foundNeutronE->Sumw2();
    h_trueNeutronEvsProtonE	= new TH2D("trueNeutronEvsProtonE","True Proton Energy vs True Neutron Energy;Primary Neutron Kinetic Energy [MeV];Proton Kinetic Energy [MeV]",60,0,300,60,0,300); 
    h_trueNeutronEvsRecoAngle	= new TH2D("trueNeutronEvsRecoAngle","True Neutron Energy vs RNN Reconstructed Angle;#Theta_{pRT} [degrees];Primary Neutron Kinetic Energy [MeV]",36,0,180,60,0,300);

    position_plots->cd();
    h_signal_ZY		= new TH2D("signal_ZY","Signal Blip Z vs Y Position;Y [cm];Z [cm]",60,-150,150,48,-100,1100);
    h_signal_ZX         = new TH2D("signal_ZX","Signal Blip Z vs X Position;X [cm];Z [cm]",60,-25,275,48,-100,1100);
    h_signal_X		= new TH1D("signal_X","Signal Blip Selection Efficiency vs X Position;X [cm];Blip Reco Efficiency",30,-25,275);
    h_signal_Y		= new TH1D("signal_Y","Signal Blip Selection Efficiency vs Y Position;Y [cm];Efficiency",30,-150,150);
    h_signal_Z		= new TH1D("signal_Z","Signal Blip Selection Efficiency vs Z Position;Z [cm];Efficiency",24,-100,1100);
    h_trueSignalX	= new TH1D("trueSignalX","True Signal Proton X Position;X [cm]",30,-25,275);
    h_trueSignalY	= new TH1D("trueSignalY","True Signal Proton Y Position;Y [cm]",30,-150,150);
    h_trueSignalZ	= new TH1D("trueSignalZ","True Signal Proton Z Position;Z [cm]",24,-100,1100);

    theta_plots->cd();
    h_trueNeutronEvsCosTheta	= new TH2D("trueNeutronEvsCosTheta","True Neutron Energy vs Cos Angle Between Neutrino Direction and Blip Position;Cos(#Theta_{sig_n});Primary Neutron Kinetic Energy [MeV]",40,-1,1,60,0,300);
    h_trueSig_theta	= new TH1D("trueSig_theta","Angle Between Neutrino Direction (Z-axis) and Vertex to Blip Position (True Signal);#Theta_{sig_n} [degrees]",50,-10,190);
    h_trueSig_cosTheta	= new TH1D("trueSig_cosTheta","Cos(Angle) Between Neutrino Direction (Z-axis) and Vertex to Blip Position (True Signal);Cos(#Theta_{sig_n})",20,-1,1);
    h_ProtonCosTheta	= new TH1D("ProtonCosTheta","Cos(Angle) Between Neutrino Direction (Z-axis) and Vertex to Blip Position (Proton Primary Ancestor);Cos(#Theta_{primancP})",20,-1,1);
    h_NeutronCosTheta	= new TH1D("NeutronCosTheta","Cos(Angle) Between Neutrino Direction (Z-axis) and Vertex to Blip Position (Neutron Primary Ancestor);Cos(#Theta_{primancN})",20,-1,1);
    h_PionCosTheta	= new TH1D("PionCosTheta","Cos(Angle) Between Neutrino Direction (Z-axis) and Vertex to Blip Position (Pi+/- Primary Ancestor);Cos(#Theta_{primancPi})",20,-1,1);
    h_OtherCosTheta	= new TH1D("OtherCosTheta","Cos(Angle) Between Neutrino Direction (Z-axis) and Vertex to Blip Position (Other Primary Ancestor);Cos(#Theta_{primancOther})",20,-1,1);

    delta_pT_plots->cd();
    h_delta_pT[0]	= new TH1D("delta_pT","Transverse Missing Momentum;#delta p_{T};nEvents",10,0,2.5);
    for(int i=1; i<4; i++) { h_delta_pT[i] = (TH1D*)h_delta_pT[0]->Clone(("delta_pT" + std::to_string(i)).c_str()); }

    shower_plots->cd();
    //h_ShrE_vs_OpAng     = new TH2D("ShrE_vs_OpAng","Shower Energy vs. Opening Angle;Opening Angle [radians];Shower Energy on Collection Plane",157,0,1.57,200,0,1000);
    h_shrEnergy		= new TH1D("shrEnergy","Reconstructed Shower Energy;Reco Energy [MeVee]",100,0,500);
    h_nBlipsInShr	= new TH1D("nBlipsInShr","Number of Blips Inside Showers per Event;nBlips",50,0,50);
    h_shrBlipE_data	= new TH1D("shrBlipE_data","Energy Spectrum of Open Data Blips Inside Showers;Reco. Blip Energy [MeVee];nBlips",110,0,55);
    h_shrBlipE_MC	= new TH1D("shrBlipE_MC","Energy Spectrum of MC Blips Inside Showers;Reco. Blip Energy [MeVee];nBlips",110,0,55);

    secondaries_plots->cd();
    h_secondPr_cosTheta	    = new TH1D("secondPr_cosTheta","Cos(Angle) Between Neutrino Direction (Z-axis) and Vertex to Blip Position (True Secondary Protons);Cos(#Theta_{secondaryP})",20,-1,1);
    //h_trueCosTheta_secondPr = new TH1D("trueCosTheta_secondPr","True Cos(Theta) (angle between true blip direction and neutrino direciton) of Secondary Protons;Cos(#Theta_{trueSecondP})",20,-1,1);

    outfile->cd();
    h_BlipsPerEvent	= new TH1D("BlipsPerEvent","Number of Blips per Event;# Blips",250,0,500);
    h_nBlips_xp0n	= new TH1D("nBlips_xp0n","Blip Multiplicity;Blip Multiplicity;A.U.",25,0,250);
    h_nBlips_xpxn	= new TH1D("nBlips_xpxn","Blip Multiplicity;Blip Multiplicity;A.U.",25,0,250);
    h_nBlips_0p0n	= new TH1D("nBlips_0p0n","Blip Multiplicity;Blip Multiplicity;A.U.",25,0,250);
    h_nBlips_0pxn	= new TH1D("nBlips_0pxn","Blip Multiplicity;Blip Multiplicity;A.U.",25,0,250);
    //h_nwires_pl0	= new TH1D("nwires_pl0","Number of Wires per Blip on Plane 0;# Wires;# Blips",18,0,18);
    //h_nwires_pl1	= new TH1D("nwires_pl1","Number of Wires per Blip on Plane 1;# Wires;# Blips",18,0,18);
    //h_nwires_pl2	= new TH1D("nwires_pl2","Number of Wires per Blip on Plane 2;# Wires;# Blips",18,0,18);
    h_maxWires		= new TH1D("maxWires","Max Wires per Blip on one Plane;# Wires",18,0,18);
    h_prTrksPerE	= new TH1D("prTrksPerE","Number of Reconstructed Objects vs True Energy of Protons From Primary Neutrons;Proton Kinetic Energy [MeV];Number of Reconstructed Objects",100,0,200);
    h_prBlipsPerE	= new TH1D("prBlipsPerE","Number of Reconstructed Objects vs True Energy of Protons From Primary Neutrons;Proton Kinetic Energy [MeV];Number of Reconstructed Objects",100,0,200);
    //h_bdt_numuCC	= new TH1D("bdt_numuCC","NDT numuCC;score?",1000,9900,10050);
    h_selectedsMomE     = new TH1D("selectedsMomE","True Energy of the Mothers of All Selected Blips",600,0,600);

    gr_signal_E_ds = new TGraph(); gr_signal_E_ds->SetNameTitle("signal_E_ds", "Blip Energy Deposition Density;Reconstructed Blip Energy (E_{blip}) [MeVee];E_{blip}/ds [MeVee/cm]");
    gr_elec_E_ds = new TGraph(); gr_elec_E_ds->SetNameTitle("gr_elec_E_ds", "Blip Energy Deposition Density;Reconstructed Blip Energy (E_{blip}) [MeVee];E_{blip}/ds [MeVee/cm]");
    gr_cosmic_E_ds = new TGraph(); gr_cosmic_E_ds->SetNameTitle("cosmic_E_ds", "Blip Energy Deposition Density;Reconstructed Blip Energy (E_{blip}) [MeVee];E_{blip}/ds [MeVee/cm]");

    float totalPOT=0;
    if(0) {
      //============================== RUN LOOP ==============================//
      int numRuns = subRunTree->GetEntries();
      std::cout<<"There are "<<numRuns<<" events in the overlay sub run tree.\n";
      std::unordered_map<int, float> map_run_pot;
      if( maxRuns < 0 ) maxRuns = numRuns;
      for(int iRun = 0; iRun < maxRuns; iRun++ ){
	subRunTree->GetEntry(iRun);
	totalPOT += pot;
	map_run_pot[runTreeRun] = pot;
      }//end run loop
    }

    if(0) {
      //============================== DATA EVENT LOOP ==============================//
      int numDataEvts = dataNuTree->GetEntries();
      std::cout<<"There are "<<numDataEvts<<" events in the data nu tree.\n";
      if( maxEvts < 0 ) maxEvts = numDataEvts;
      for(int iEvt = 0; iEvt < maxEvts; iEvt++ ){     
	dataNuTree->GetEntry(iEvt);
	
	std::vector<std::pair<TVector3, TVector3>> data_shr_vtx_P;
	for(size_t i_pfp = 0; i_pfp < n_pfps_data; i_pfp++){
	  if( data_trk_score_v->at(i_pfp) < 0.45 ) {
	    TVector3 data_shrStart(data_shr_start_x_v->at(i_pfp), data_shr_start_y_v->at(i_pfp), data_shr_start_z_v->at(i_pfp));
	    TVector3 data_shr_p(data_shr_px_v->at(i_pfp), data_shr_py_v->at(i_pfp), data_shr_pz_v->at(i_pfp));
	    data_shr_vtx_P.push_back( std::make_pair(data_shrStart, data_shr_p) );
	  }
	}
	
	for(size_t i=0; i<nblips_saved_data; i++){
	  float blipX_data = blip_x_data->at(i); float blipY_data = blip_y_data->at(i); float blipZ_data = blip_z_data->at(i);
	  float blipE_data = blip_energy_data->at(i); TVector3 blip_loc_data( blipX_data, blipY_data, blipZ_data );
	  for( const auto& [shr_vtx, shr_P] : data_shr_vtx_P ) {
	    if( IsWithinShowerCone(shr_vtx, shr_P, blip_loc_data, 90.0) ) { h_shrBlipE_data->Fill(blipE_data); }
	  }
	}
      }//end data event loop
    }

    float minX, minY, minZ = 0;
    float maxX, maxY, maxZ = 0;
 
    int numEvents = overlayNuTree->GetEntries();
    std::cout<<"There are "<<numEvents<<" events in the overlay nu tree.\n"; 
    //std::cout<<"There are "<<wcTree->GetEntries()<<" events in the WC tree.\n";  

    int totalBlips=0;
    float totalBlipE=0;
    int totalSelectedBlips=0;   // can't be float due to floating-point precision limit
    float matchedBlips=0;
    float neutronIDprBlips=0;
    //float oldMethodNeutronEff=0;
    float signal_blips=0;
    int selectedEvents=0;
    int nonSelected=0;
    int numEvents_numuCC=0;
    int numTrueEvents_numuCC=0;
    float invalid_RNN_blips=0;
    float numBlipsAlsoTrks=0;
    float numSigBlipsAlsoTrks=0;
    std::vector<int> nShrs_per_numuCC_event;

    int total_nu_Ns=0;
    float nTrueSignalEvents=0; float nTotalTaggedEvents=0; float nTotalTrueTaggedEvents;
    float trueSignalPOT=0; float totalTaggedPOT=0;
    std::unordered_map<int, std::vector<int>> map_n_p;
    int category;

    //============================== MC EVENT LOOP ==============================//
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      //if( iEvent > 10000 ) break;
      //if( iEvent < 100000  && iEvent % 2 == 0 ) continue;
      selectedEvents++;
      overlayNuTree->GetEntry(iEvent);

      int all_mc_length = static_cast<int>(all_mc_trkid->size());
      if(iEvent < 10 || iEvent % 10000 == 0) {
	std::cout<<"\n=======================================\n";
	std::cout<<"Reading event "<<iEvent<<" out of "<<numEvents<<"\n";
	std::cout<<"NBlips: "<<nblips_saved<<"   NParticles (all_mc): "<<all_mc_length<<"\n";
      }

      std::set<int> tagged_neutron_trkIDs;
      bool thisTrueEvtWasTagged = false; bool thisTotalEvtWasTagged = false;
 
      TVector3 reco_nu_vtx( reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z );
      float cosmicProx = _closestNuCosmicDist;

      map_allTrkID_partIDX.clear();
      map_allTrkID_momTrkID.clear();
      map_n_p.clear();

      //if( IsInFV_Buffer(reco_nu_vtx) == false ) {
      //std::cout<<"Neutrino vtx location: ("<<reco_nu_vtx.X()<<", "<<reco_nu_vtx.Y()<<", "<<reco_nu_vtx.Z()<<
      //")\nIs in FV? "<<std::boolalpha<<IsInFV_Buffer(reco_nu_vtx)<<"\n\n"; }

      //std::cout<<"\nccnc: "<<ccnc<<"\nnu_pdg: "<<nu_pdg<<"\n";
      if( ccnc == 0 && nu_pdg == 14 ) numTrueEvents_numuCC++;   // truth CC and numu pdg
      //============================== NUMUCC INCLUSIVE SELECTION ==============================//
      // Preselection
		if( 1 && ( slice_orig_pass_id != 1 ||   // replaces nslice != 1
      IsInFV_BufferAndDeadZone( reco_nu_vtx ) != 1 ||
		     slice_orig_topo_score <= 0.06 ||   //replaces topological_score
				      crtveto == 1 ||
			 _closestNuCosmicDist <= 5 ) ) { continue; }   // not SCE in ubana code

      //============================== PFP LOOP ==============================//
      int n_shrs = 0; std::vector<float> trk_id;
      bool numuCC = false;   // mip = muon, hip = proton
      std::vector<std::tuple<TVector3, TVector3, float>> shr_vtx_P_E; //std::vector<std::pair<TVector3, TVector3>> ShrVtx_ShrP;
      for(size_t i_pfp = 0; i_pfp < n_pfps; i_pfp++){

	// Muon Tracks   (could we assume this is always the leading track?)
	if( trk_score_v->at(i_pfp) > 0.8 ) {   // 0.45 if 0.8 checked elsewhere
	  TVector3 trkStart(trk_sce_start_x_v->at(i_pfp), trk_sce_start_y_v->at(i_pfp), trk_sce_start_z_v->at(i_pfp));
	  TVector3 trkEnd(trk_sce_end_x_v->at(i_pfp), trk_sce_end_y_v->at(i_pfp), trk_sce_end_z_v->at(i_pfp));
	  float trk_length = (trkEnd - trkStart).Mag();

	  if( IsInFV_Buffer(trkStart) && IsInFV_Buffer(trkEnd) &&
	      std::min( (trkStart-reco_nu_vtx).Mag(), (trkEnd-reco_nu_vtx).Mag() ) < 4 &&
	      trk_llr_pid_score_v->at(i_pfp) > 0.2 &&   //pfng2semlabel->at(i_pfp) == 0 &&   //makes efficiency worse
	      trk_length > 10 )
	    {
	    //if( numuCC == true ) std::cout<<"Jeez, two big 'ol muons?!"<<std::endl;
	    numuCC = true;
	    }
	  trk_id.push_back(trk_pfp_id_v->at(i_pfp));
	}

	//if( pfng2semlabel->at(i_pfp) == 0 ) {
	//  muon_mom_mag = trk_mcs_muon_mom_v->at(i_pfp);
	//  muon_mom_vec = (trkEnd - trkStart)*muon_mom_mag;   // is this correct? trkStart/End are unit vecs?
	//}
	//if( pfng2semlabel->at(i_pfp) == 1 ) {
	//  proton_momentum = sqrt( E*E - TMath::Sq(938.3) ) / c ;   // proton rest energy (m0*c^2) = ~938.3 MeV
	//}

	// Showers WHAT HAPPENS WHEN shower information is not available?!?Q!@
	if( trk_score_v->at(i_pfp) < 0.45 ) {
	  TVector3 shrStart(shr_start_x_v->at(i_pfp), shr_start_y_v->at(i_pfp), shr_start_z_v->at(i_pfp));
	  TVector3 shr_p(shr_px_v->at(i_pfp), shr_py_v->at(i_pfp), shr_pz_v->at(i_pfp));   // already normalized
	  //TVector3 shr_distVec = shr_p * shr_dist_v->at(i_pfp);   // ( |Px|*d, |Py|*d, |Pz|*d )
	  //TVector3 shrEnd = shrStart + shr_distVec;
	  float shrE = shr_energy_y_v->at(i_pfp);   // energies don't have 1/0.83 correction factor applied
	  h_shrEnergy->Fill(shrE);
	  shr_vtx_P_E.emplace_back(shrStart, shr_p, shrE); //push_back(std::make_tuple(shrStart, shr_p, shrE));

	  n_shrs++;

	  //float opAng = shr_openangle_v->at(i_pfp);
	  //if( opAng > 1e-3 && opAng <= 1.57 ) h_ShrE_vs_OpAng->Fill( opAng, shrE );
	  //std::cout<<"Opening Angle V: "<<shr_openangle_v->at(i_pfp)<<
          //	       "\nShower (dx, dy, dz) = ("<<shr_distVec.x()<<", "<<shr_distVec.y()<<", "<<shr_distVec.z()<<
          //           "\nStart (x, y, z) = ("<<shrStart.x()<<", "<<shrStart.y()<<", "<<shrStart.z()<<
          //           "\nEnd (x, y, z) = ("<<shrEnd.x()<<", "<<shrEnd.y()<<", "<<shrEnd.z()<<")\n\n";
	}
      }
      if( 1 && numuCC == false ) continue;
      numEvents_numuCC++; nShrs_per_numuCC_event.push_back(n_shrs);
      int num_nu_Ns = 0; bool hasFSprot=false; int blipsWithin100=0;
      //============================== PARTICLE LOOP ==============================//
      for(int i=0; i<all_mc_length; i++){

	int mc_trkID = all_mc_trkid->at(i); float mc_E = all_mc_E->at(i);
	int mc_pdg = all_mc_pdg->at(i); int mc_momTrkID = all_mc_mother->at(i);
	bool mc_isPrimary = mc_momTrkID == 0;

	map_allTrkID_partIDX[mc_trkID] = i;
	map_allTrkID_momTrkID[mc_trkID] = mc_momTrkID;

	// Printout this particle information to the screen 
	if( 1 && iEvent<10 ){   //&& part_isContained[i]
	  printf("%4i, TrackID: %-4i   PDG: %-12i E=%7.3f MeV,   moth=%5i\n",//, %18s \n",//   end(X,Y,Z)=(%1.2f, %1.2f, %1.2f)\n",   // dL=%6.1f cm, part_pathlen[i], other way to figure out dL with this file?
	    i, mc_trkID, mc_pdg, mc_E, mc_momTrkID//, all_mc_process->at(i).c_str()//, all_mc_endx->at(i), all_mc_endy->at(i), all_mc_endz->at(i)
	  );
	}

        // -------------------------------------------------------------
        // Identify any neutrons and save their information
        if ( mc_pdg == 2112 ) {
          if( mc_isPrimary ) {   // primary neutrons
            total_nu_Ns++; num_nu_Ns++;
	    h_trueNeutronE->Fill( (mc_E - 0.93957)*1000 );   // subtract neutron rest mass energy and convert from GeV to MeV
            map_n_p[mc_trkID].clear();
          }
        }

	auto it = map_n_p.find(mc_momTrkID);
	if( mc_pdg == 2212 ) {   // protons
	  if( mc_isPrimary ) hasFSprot = true;
	  else if( it != map_n_p.end() ) {   // from primary neutron mothers
	    it->second.push_back(mc_trkID); h_trueSignalX->Fill(all_mc_endx->at(i));
	    if( std::find( trk_id.begin(), trk_id.end(),  mc_trkID ) != trk_id.end() ) h_prTrksPerE->Fill( (mc_E- 0.93827)*1000 );
	    //if( (all_mc_E->at(map_allTrkID_partIDX[mc_momTrkID]) - 0.93957)*1000 < 5 ) std::cout<<"LESS THAN 5 MeV"<<std::endl; else std::cout<<"more than 5 MeV"<<std::endl;
	    //std::cout<<"FS neutron energy: "<<(all_mc_E->at(map_allTrkID_partIDX[mc_momTrkID]) - 0.93957)*1000<<"\nDaughter proton energy: "<<(all_mc_E->at(i) - 0.93827)*1000<<"\n";
	    //std::cout<<"Process: "<<all_mc_process->at(i).c_str()<<"\n\n";
            h_trueNeutronEvsProtonE->Fill( (all_mc_E->at(map_allTrkID_partIDX[mc_momTrkID]) - 0.93957)*1000, ( mc_E - 0.93827)*1000 );
	  }
        }

      }//end particle loop

      if( num_nu_Ns == 0 ) h_delta_pT[0]->Fill(pt);
      else if( num_nu_Ns == 1 ) h_delta_pT[1]->Fill(pt);
      else if( num_nu_Ns == 2 ) h_delta_pT[2]->Fill(pt);
      else h_delta_pT[3]->Fill(pt);

      bool evtIs_xp0n=0, evtIs_0p0n=0, evtIs_xpxn=0, evtIs_0pxn=0;
      if( num_nu_Ns == 0 ) {
	if( hasFSprot ) evtIs_xp0n=1;
	  else evtIs_0p0n=1;
      } else {
	if( hasFSprot ) evtIs_xpxn=1;
	  else evtIs_0pxn=1;
      }

      // cut booleans for ease of selection
      bool preselect = 1; bool RNNvalidCut = 1; bool nuProximity = 1;
      bool vtxIso = false; bool energyCut = 0; bool densityCut = 0;
      bool shrCut = 1; bool BDTcut = 1;
      int nBlipsInShr=0;
      //============================== BLIP LOOP ==============================//
      for(size_t i=0; i<nblips_saved; i++){

        totalBlips++; category=blip_true_ncategory->at(i);
	if( category == -9 ) category = 9;

	float blipX = blip_x->at(i), blipY = blip_y->at(i), blipZ = blip_z->at(i);
	if( minX > blipX ) minX=blipX; if( minY > blipY ) minY=blipY; if( minZ > blipZ ) minZ=blipZ;
	if( maxX < blipX ) maxX=blipX; if( maxY < blipY ) maxY=blipY; if( maxZ < blipZ ) maxZ=blipZ;
        TVector3 blip_loc( blipX, blipY, blipZ ); TVector3 reco_trace_vec( blip_loc - reco_nu_vtx );
	float dist = reco_trace_vec.Mag(); if( dist >= 100 ) blipsWithin100++; float Zdist = abs(blipZ - reco_nu_vtx.Z());
	float reco_theta = reco_trace_vec.Theta(); float cos_theta = TMath::Cos(reco_theta);   //.Angle(TVector3(0,0,1)); = .Theta() 

	float blipE = blip_energy->at(i); float trueE = blip_true_energy->at(i);
	totalBlipE+=blipE; float E_res = abs(blipE - trueE)/trueE;

	TVector3 blip_true_dir( blip_true_dir_x->at(i), blip_true_dir_y->at(i), blip_true_dir_z->at(i) );
	//float cos_blipTrueDir_nuDir = TMath::Cos(blip_true_dir.Theta());   // blip true direction dot product with Z-axis

	float dx = blip_dx->at(i); float dw = blip_dw->at(i);
	float size = sqrt( dx*dx + dw*dw );   // direct multiplication is fastest ig

	float density = blipE/size; float dist_to_trk = blip_proxtrkdist->at(i);
	float cutoff = cutoff_graph->Eval(blipE);

	int pdg = blip_true_pdg->at(i);	int g4id = blip_true_g4id->at(i);
	int primanc_pdg = blip_true_primancPDG->at(i);

	// This could be dangerous when accessing vectors, may want to switch to better logic
	// Can't just do "if( g4id < 0 ) {" to set these since photons,e+,e- also have no mc_truth info. saved
	int allIDX = -9; bool isPrimary = false; 
	int mom_trkID = -9; bool mom_isPrimary = false;
	int mom_allIDX = -9; int mom_pdg = -9;
	float mom_E = -9;
	auto it = map_allTrkID_partIDX.find(g4id);
	if( it != map_allTrkID_partIDX.end() ) {
	  allIDX = it->second;
	  mom_trkID = all_mc_mother->at( allIDX );
	  isPrimary = mom_trkID == 0;   // if mother == 0, particle is primary, could also use process="primary"
 
	  if( !isPrimary && map_allTrkID_partIDX.count(mom_trkID) != 0 ) {
	    mom_allIDX = map_allTrkID_partIDX[mom_trkID]; mom_E = (all_mc_E->at(mom_allIDX) - 0.93957)*1000;
	    mom_isPrimary = all_mc_mother->at( mom_allIDX ) == 0; mom_pdg = all_mc_pdg->at( mom_allIDX );
	  }
	}

	if( mom_isPrimary && pdg == 2212 && mom_pdg == 2112 ) {  // signal blips (protons from primary neutron moms)
	  h_signal_ZY->Fill(blipY, blipZ); h_signal_ZX->Fill(blipX, blipZ); h_trueSignalY->Fill(blipY); h_trueSignalZ->Fill(blipZ);

	  h_trueSig_theta->Fill(reco_theta*TMath::RadToDeg()); h_trueSig_cosTheta->Fill( cos_theta );
	  h_trueNeutronEvsCosTheta->Fill( cos_theta, mom_E );

	  if( !thisTrueEvtWasTagged ) { nTrueSignalEvents++; thisTrueEvtWasTagged=true; } //trueSignalPOT += map_run_pot[run]; }

	  h_prBlipsPerE->Fill(trueE);

	  if( blip_trkid->at(i) > -1 ) numSigBlipsAlsoTrks++;
	}

	if( blip_trkid->at(i) > -1 ) numBlipsAlsoTrks++;

	// fill true secondary proton plots
	if( category == 3 || category == 4 ) {   // can do this same category "if" for signal
	  h_secondPr_cosTheta->Fill(cos_theta);  //h_trueCosTheta_secondPr->Fill( trueCosTheta );
	  if( 0 && cos_theta > 0.5 ) {
	    std::cout<<"\nReco Cos(theta) = "<<cos_theta<<"\nTrue Cos(theta) of Neutron Mom = "//<<trueCosN//<<"\nBlip PDG: "<<pdg<<
	    <<"\nMother PDG: "<<mom_pdg<<"\nMother energy = "<<mom_E;//<<"\nMother is Primary? "<<std::boolalpha<<mom_isPrimary;
	    std::cout<<"\nPrimary ancestor PDG: "<< [primanc_pdg]() -> std::string {
	      if ( primanc_pdg == 2212 ) return "Proton"; if ( primanc_pdg == 2112 ) return "Neutron"; if ( primanc_pdg == 211 )  return "Pi+";
	      return std::to_string(primanc_pdg);
	    }() << "\n";
	  }
	}

        //============================== CUTS ==============================//
	bool nearDeath = false; bool isProton = true; bool isInShower = false; bool pred_sig = true;
	if( preselect ) {
	  if( blip_touchtrk->at(i) || dist_to_trk <= 1 || !IsInFV_Exact( blip_loc ) ) continue;
	  nearDeath = ( blip_pl0_bydeadwire->at(i) == 1 || blip_pl1_bydeadwire->at(i) == 1 || blip_pl2_bydeadwire->at(i) == 1 );
	}
	if( nuProximity && dist >= 100 ) continue;

	// Plot the energy of the mother of signal blips that passed the preselection & nu proximity cut
	if( pdg == 2212 && mom_pdg == 2112 && mom_isPrimary ){ h_preFoundNeutronE->Fill( mom_E ); }

	if( vtxIso && dist <= 5 ) continue;
	if( energyCut && blipE <= 3 ) continue;
	if( densityCut && density < Fit_Function_PID_Cut(blipE) ) isProton = false; //if( blipE <= 10 && density < cutoff ) isProton = false;

	//double minDist2Shr = 100000; TVector3 blip_to_shrVtx(0,0,0); double blip_shr_angle = 0;
	for( const auto& [shr_vtx, shr_P, E] : shr_vtx_P_E ) {
	  //if ( IsWithinShowerCone(shr_vtx, shr_P, blip_loc, 90.0) && (blip_loc - shr_vtx).Mag() < minDist2Shr ) {
	  //  blip_to_shrVtx = blip_loc - shr_vtx;
	  //  blip_shr_angle = blip_to_shrVtx.Angle(shr_P);
	  //  minDist2Shr = blip_to_shrVtx.Mag();
          //}
	  if( IsWithinShowerCone(shr_vtx, shr_P, blip_loc, 90.0) ) {
	    nBlipsInShr++; if( shrCut ) isInShower = true; h_shrBlipE_MC->Fill(blipE);
	  }
	}
	// To manually check if a blip is deemed to be inside any shower in this event:
	//std::cout<<"\nShower vertex: ("<<vtx.X()<<", "<<vtx.Y()<<", "<<vtx.Z()<<")"
	//<<"\nShower direction: ("<<P.X()<<", "<<P.Y()<<", "<<P.Z()<<")"; }
	//std::cout<<"Blip location: ("<<blip_loc.X()<<", "<<blip_loc.Y()<<", "<<blip_loc.Z()<<")";
	//std::cout<<"\nBlip is in shower: "<<std::boolalpha<<isInShower<<"\n\n";

        // RNN model prediction
        double res_angle = -1.0; TVector3 reco_dir; float align_angle, offset;
	bool RNN_isValid = blip_rnn_dir_isValid->at(i);
        if( RNN_isValid ) {
          reco_dir.SetXYZ( blip_rnn_dir_x->at(i), blip_rnn_dir_y->at(i), blip_rnn_dir_z->at(i) );
          align_angle = reco_dir.Angle(reco_trace_vec); offset = abs(TMath::Cos(align_angle));   // pointing offset
	  res_angle = reco_dir.Angle(blip_true_dir)*TMath::RadToDeg();
          if( res_angle > 90 ) res_angle = 180-res_angle;
        }

	// BDT prediction; Compute output scores for binary classification and convert 
	blip_E=blipE; dist2Vtx=dist; dist2Trk=dist_to_trk; dist2Edge=distToNearestEdge( blip_loc );
	blip_size=size; blip_density=density; alignment_angle=align_angle; Z_dist2Vtx=Zdist;
	//dist2ShrVtx=blip_to_shrVtx.Mag(); angle_blip_shr=blip_shr_angle;
	//blip_X=blipX; blip_Y=blipY; blip_Z=blipZ; blip_dX=dx; blip_dW=dw;
	//blip_dir_x=blip_rnn_dir_x->at(i); blip_dir_y=blip_rnn_dir_y->at(i); blip_dir_z=blip_rnn_dir_z->at(i);
	double tmva_score; double base_score; double XGBoost_score;
        tmva_score = reader->EvaluateMVA( "BDT" ); // I think this is where : Rebuilding Dataset Default
	base_score = 0.009503577;
	XGBoost_score = 0.5 * TMath::Log( (1 + tmva_score)/(1 - tmva_score) ) - base_score;   // = atanh(tmva_score), sigmoid??
        if( BDTcut && XGBoost_score < -1.1 ) pred_sig = false;
	if( 0 ) {   // turn on if using 2 different BDTs
	  if( blipE > 3 ) {
	    tmva_score = reader->EvaluateMVA( "highE_BDT" );    base_score = 0.02169922;
	    XGBoost_score = 0.5 * TMath::Log( (1 + tmva_score)/(1 - tmva_score) ) - base_score;
	    if( BDTcut && XGBoost_score < -1.2 ) pred_sig = false;
	  } else {
	    tmva_score = reader->EvaluateMVA( "lowE_BDT" );
	    base_score = 0.005817941;
	    XGBoost_score = 0.5 * TMath::Log( (1 + tmva_score)/(1 - tmva_score) ) - base_score;
	    if( BDTcut && XGBoost_score < -2.6 ) pred_sig = false;
	  }
	}


	if( !RNN_isValid ) {
	  if( !isInShower && !nearDeath && pred_sig ) invalid_RNN_blips++;
	  if( RNNvalidCut ) continue;
	}
	// convert alignment angle for plotting purposes (BDT uses raw radian form)
	if( blip_rnn_dir_isValid->at(i) ) {
	  align_angle = align_angle*TMath::RadToDeg();
	  if( align_angle > 90 ) align_angle = 180-align_angle;
	}

	// cosmics
        if( g4id < 0 && !isInShower && !nearDeath && pred_sig && blipE <= 10 ) { gr_cosmic_E_ds->SetPoint(gr_cosmic_E_ds->GetN(), blipE, density); }

	if( abs(pdg) != 11 && abs(pdg) != 22 ) { matchedBlips++; }
	else if( !isInShower && !nearDeath && pred_sig && blipE <= 10 ) {   // photons, electrons, and positrons
	  gr_elec_E_ds->SetPoint(gr_elec_E_ds->GetN(), blipE, density);
	}

	//if( blip_true_pdg->at(i) != all_mc_pdg->at(allIDX) ) std::cout<<"WOAH, HOLD UR HORSES! blip_true_pdg && all_mc_pdg don't match!!"<<endl;

	//h_nwires_pl0->Fill( blip_pl0_nwires->at(i) < 1 ? 0 : blip_pl0_nwires->at(i) );
	//h_nwires_pl1->Fill( blip_pl1_nwires->at(i) < 1 ? 0 : blip_pl1_nwires->at(i) ); 
	//h_nwires_pl2->Fill( blip_pl2_nwires->at(i) < 1 ? 0 : blip_pl2_nwires->at(i) ); 
	//int maxWires = std::max({blip_pl0_nwires->at(i), blip_pl1_nwires->at(i), blip_pl2_nwires->at(i)});
	//h_maxWires->Fill(maxWires);

        // PROTONS
	if( pdg == 2212 && !isPrimary && mom_pdg == 2112 && mom_isPrimary ) {   // non-primary protons from primary neutron moms (should all be nIne)
	  if( blipE <= 10 && !nearDeath && !isInShower && pred_sig ) gr_signal_E_ds->SetPoint(gr_signal_E_ds->GetN(), blipE, density);
	  if( !nearDeath && !isInShower && isProton && pred_sig ) {   // SIGNAL PROTONS
	    signal_blips++; h_signal_X->Fill(all_mc_endx->at(allIDX)); h_signal_Y->Fill(blipY); h_signal_Z->Fill(blipZ);
	    h_foundNeutronE->Fill( mom_E );
	    //auto it = std::find(tagged_neutron_trkIDs.begin(), tagged_neutron_trkIDs.end(), mom_trkID);
            //if( it == tagged_neutron_trkIDs.end() ) { neutronIDprBlips++; tagged_neutron_trkIDs.push_back(mom_trkID); }
	    if( tagged_neutron_trkIDs.insert(mom_trkID).second ) { neutronIDprBlips++; }   // true if newly inserted, same as above
	    //if( !map_n_p[mom_trkID].empty() && map_n_p[mom_trkID][0] == g4id ) oldMethodNeutronEff++;
	  }
	}

	//if( category == 0 ) { std::cout<<"\nOH NAUR!, naur catagaury set ;( \nPDG = "<<p_pdg<<std::endl; }

	float blip_dir = reco_dir.Angle(reco_trace_vec);
	if( mom_pdg == 2112 && mom_isPrimary && RNN_isValid ) {
	  h_trueNeutronEvsRecoAngle->Fill( blip_dir*TMath::RadToDeg(), mom_E );
	}

	if( isProton ) {
	  if( !isInShower ) {
	    if( !nearDeath ) {
	      if( pred_sig ) {
		if( category == 3 || category == 4 ) {
		  if( primanc_pdg == 2212 ) h_ProtonCosTheta->Fill(cos_theta);
		  else if( primanc_pdg == 2112 ) h_NeutronCosTheta->Fill(cos_theta);
		  else if( abs(primanc_pdg) == 211 ) h_PionCosTheta->Fill(cos_theta); else { h_OtherCosTheta->Fill(cos_theta); std::cout<<"\nOther PDG: "<<primanc_pdg; }
		}
		if( !thisTotalEvtWasTagged ) {
		  nTotalTaggedEvents++; if( thisTrueEvtWasTagged ) nTotalTrueTaggedEvents++;
		  thisTotalEvtWasTagged=true; //totalTaggedPOT += map_run_pot[run];
		}
		totalSelectedBlips++; if( mom_E >= 0 ) h_selectedsMomE->Fill( mom_E );
		// Fill "all blips" histogram
                h_recoVtxDist[10]->Fill(dist); h_blipE[10]->Fill(blipE);
                h_Edensity[10]->Fill(density); h_energy_res[10]->Fill(E_res);
		//std::cout<<"Blip true energy: "<<trueE<<"\nBlip reco energy: "<<blipE<<"\n\n";
                h_proxTrkDist[10]->Fill(dist_to_trk); h_cosmicDistToVtx[10]->Fill(cosmicProx);
                h_dx[10]->Fill(dx); h_dw[10]->Fill(dw); h_size[10]->Fill(size);
		h_cos_theta[10]->Fill(cos_theta); h_cos_theta[category]->Fill(cos_theta);
		h_cosAlignment[10]->Fill(offset); h_alignmentAngle[10]->Fill(align_angle);
		// Fill specific categories
		h_recoVtxDist[category]->Fill(dist); h_blipE[category]->Fill(blipE);
		h_Edensity[category]->Fill(density); h_energy_res[category]->Fill(E_res);
		h_proxTrkDist[category]->Fill(dist_to_trk); h_cosmicDistToVtx[category]->Fill(cosmicProx);
		h_dx[category]->Fill(dx); h_dw[category]->Fill(dw); h_size[category]->Fill(size);
		h_cosAlignment[category]->Fill(offset); h_alignmentAngle[category]->Fill(align_angle);
		if( allIDX >= 0 && g4id >= 0 ) {
                  h_resolution[10]->Fill(res_angle); //h_trueCosTheta_vs_momTrueE[10]->Fill(mom_E, trueCosTheta);
                  h_res_vs_E[10]->Fill((all_mc_E->at(allIDX)-0.93827)*1000, res_angle);
		  h_resolution[category]->Fill(res_angle); //h_trueCosTheta_vs_momTrueE[category]->Fill(mom_E, trueCosTheta);
		  h_res_vs_E[category]->Fill((all_mc_E->at(allIDX)-0.93827)*1000, res_angle);
		}
	      }
	      h_BDT_score[10]->Fill(XGBoost_score);
	      h_BDT_score[category]->Fill(XGBoost_score);
	    }
	    if( pred_sig ) {
              h_byDeadWire[10]->Fill(nearDeath);
	      h_byDeadWire[category]->Fill(nearDeath);
	    }
	  }
	  if( !nearDeath && pred_sig ) {
            h_inShower[10]->Fill(isInShower);
	    h_inShower[category]->Fill(isInShower);
	  }
	}

      }//end blip loop

      if( evtIs_xp0n ) h_nBlips_xp0n->Fill(totalSelectedBlips);
	else if( evtIs_0p0n ) h_nBlips_0p0n->Fill(totalSelectedBlips);
	  else if( evtIs_xpxn ) h_nBlips_xpxn->Fill(totalSelectedBlips);
	    else h_nBlips_0pxn->Fill(totalSelectedBlips); //std::cout<<evtIs_0pxn<<std::endl; }

      //std::cout<<"This event is: ";
      //std::map<std::string, bool> variables = {
      //    {"xp0n", evtIs_xp0n},
      //    {"0p0n", evtIs_0p0n},
      //    {"xpxn", evtIs_xpxn},
      //    {"0pxn", evtIs_0pxn}
      //};
      //for (const auto& pair : variables) { if (pair.second == true) std::cout << pair.first << std::endl; }

      h_BlipsPerEvent->Fill(nblips_saved); if( nBlipsInShr > 0 ) h_nBlipsInShr->Fill(nBlipsInShr);

    }//end events loop

    float events_w_overOneShr = std::count_if(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end(), greater_than_1);
    float events_w_oneShr = std::count(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end(), 1);
    float events_w_noShr = std::count(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end(), 0);

    if( nblips_saved == 0) std::cout<<"No blips saved"; int numTrueSigBlips = h_trueSignalY->GetEntries();

    printf("\nBlip multiplicity: %8.3f +/- %.3f\n",
	   double(totalBlips)/double(selectedEvents), h_BlipsPerEvent->GetRMS(11) );

    std::cout<<"\nTotal Events: "<<numEvents<<" events\n"
             //<<"After removing training events... Total Events = "<<selectedEvents<<" events\n"
             <<"NumuCC inclusive selected events: "<<numEvents_numuCC<<" events\n"
	     <<"True NumuCC inclusive events: "<<numTrueEvents_numuCC<<" events\n"
	     <<"TRUE signal events: "<<nTrueSignalEvents<<" events\n"//, Sum POT = "<<trueSignalPOT
	     //<<", numEvents/POT = "<<nTrueSignalEvents/trueSignalPOT<<"\n"
	     <<"Total TAGGED events: "<<nTotalTaggedEvents<<" events\n"//, Sum POT = "<<totalTaggedPOT
	     //<<", numEvents/POT = "<<nTotalTaggedEvents/totalTaggedPOT<<"\n"
	     <<"Total TRUE TAGGED events: "<<nTotalTrueTaggedEvents<<" events\n"
	     <<"Total sample POT: "<<totalPOT<<"\n\n"
	     //<<"Selected events without blips: "<<nonSelected<<" events\n"
             //<<"Percent of selected events without reco nu vtx: "<<nonSelected*100/(selectedEvents+nonSelected)<<"%\n\n"

	     <<"Average number of showers per numuCC event: "<<TMath::Mean(nShrs_per_numuCC_event.begin(), nShrs_per_numuCC_event.end())<<" showers\n"
	     //<<"Percent of numuCC events with one shower: "<<100*events_w_oneShr/numEvents_numuCC<<"%\n"
	     //<<"Percent of numuCC events with more than one shower: "<<100*events_w_overOneShr/numEvents_numuCC<<"%\n"
             //<<"Percent of numuCC events without any showers: "<<100*events_w_noShr/numEvents_numuCC<<"%\n\n"

	     <<"Total blips: "<<totalBlips<<" blips\n"
	     <<"Number of true signal blips: "<<numTrueSigBlips<<" blips\n"
	     <<"Number of selected signal blips: "<<signal_blips<<" blips\n"
	     <<"Total selected blips: "<<totalSelectedBlips<<" blips\n"
	     <<"Number of neutrons tagged by selected signal blips: "<<neutronIDprBlips<<" blips\n"
	     <<"Total neutrino-induced neutrons: "<<total_nu_Ns<<" neutrons\n"
	     <<"Percent of blips also reconstructed as tracks: "<<numBlipsAlsoTrks/totalBlips*100<<"%\n"
             <<"Percent of signal blips also reconstructed as tracks: "<<numSigBlipsAlsoTrks/numTrueSigBlips*100<<"%\n"
             //<<"Percent of total blips with a true g4id ALSO not protons/electrons: "<<(matchedBlips*100)/totalBlips<<"%\n\n"

	     <<"Percent of selected blips with invalid RNN prediction: "<<(invalid_RNN_blips*100)/totalBlips<<"%\n\n";

    float eff = neutronIDprBlips/total_nu_Ns; float pur = signal_blips/totalSelectedBlips; float evt_pur = nTotalTaggedEvents/nTrueSignalEvents;
    float eff_unc = eff*sqrt( 1/neutronIDprBlips + 1/total_nu_Ns ); float pur_unc = pur*sqrt( 1/signal_blips + 1/totalSelectedBlips );
    //\nOld Efficiency: %5.5f%% +/- %.2f	oldMethodNeutronEff*100/total_nu_Ns, 1/sqrt(oldMethodNeutronEff),
    auto [belEff, belPur, numBel] = effPurBelowEnergy( h_foundNeutronE, h_trueNeutronE, h_selectedsMomE, 100 );
    auto [abEff, abPur, numAb] = effPurAboveEnergy( h_foundNeutronE, h_trueNeutronE, h_selectedsMomE, 100 );
    auto [ab40Eff, ab40Pur, numAb40] = effPurAboveEnergy( h_foundNeutronE, h_trueNeutronE, h_selectedsMomE, 40 );
    auto [ab50Eff, ab50Pur, numAb50] = effPurAboveEnergy( h_foundNeutronE, h_trueNeutronE, h_selectedsMomE, 50 );
    auto [ab75Eff, ab75Pur, numAb75] = effPurAboveEnergy( h_foundNeutronE, h_trueNeutronE, h_selectedsMomE, 75 );
    printf( "UNDER 100 MeV efficiency and purity after selection: %5.3f%% eff., %5.3f%% pur. (%i)   \nOVER 100 MeV efficiency and purity after selection: %5.3f%% eff., %5.3f%% pur. (%i)      \nOVER 75 MeV: %5.3f%% eff., %5.3f%% pur. (%i)      \nOVER 50 MeV: %5.3f%% eff., %5.3f%% pur. (%i)      \nOVER 40 MeV: %5.3f%% eff., %5.3f%% pur. (%i)\n\n",
	     belEff, belPur, numBel, abEff, abPur, numAb, ab75Eff, ab75Pur, numAb75, ab50Eff, ab50Pur, numAb50, ab40Eff, ab40Pur, numAb40
	  );
    //\nOVER 40 MeV and UNDER 250 MeV efficiency after selection: %5.3f%%   \n\n", effAboveBelowEnergy( h_foundNeutronE, h_trueNeutronE, 40, 250 )
    //(BEFORE selection: %5.3f%% eff.)->effPurBelowEnergy( h_preFoundNeutronE, h_trueNeutronE, 100 ),

    printf( "Efficiency: %5.3f%% +/- %.3f%%   \nPurity: %9.3f%% +/- %.2f%%   \nE*P: %13.7f +/- %.4f   \nEvent-Level Purity: %9.3f%% +/- %.2f%%\n\n",   //   \n2*(P*E)/(P+E): %7.5f\n\n",
	     eff*100, eff_unc*100,
	     pur*100, pur_unc*100,
	     eff*pur, sqrt( TMath::Sq(pur*eff_unc) +  TMath::Sq(eff*pur_unc) ),
	     evt_pur*100, evt_pur*sqrt( 1/nTotalTaggedEvents + 1/nTrueSignalEvents )*100
	     //2*(pur*eff)/(pur+eff)
	     //eff*(pur*pur), 1/sqrt(neutronIDprBlips + signal_blips*signal_blips)
          );

   	// could also do (max - min)/2 to estimate statistical uncertainty

    if( 1 ) { 
      //h_shrBlipE_data->Scale(1.0 / h_shrBlipE_data->Integral("width")); h_shrBlipE_MC->Scale(1.0 / h_shrBlipE_MC->Integral("width"));
      //merge_plot2(h_shrBlipE_data, "Data", h_shrBlipE_MC, "MC", "Energy Spectrum of Blips Inside Showers");
      //merge_plot2(h_prBlipsPerE, "Blips", h_prTrksPerE, "Tracks", "title");
      merge_plot2( h_prBlipsPerE, TString::Format( "%s (%.1f%%)", "Blips", (h_prBlipsPerE->GetEntries()/numTrueSigBlips)*100 ),
		   h_prTrksPerE, TString::Format( "%s (%.1f%%)", "Tracks", (h_prTrksPerE->GetEntries()/numTrueSigBlips)*100 ), "");

      //TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
      //h_ShrE_vs_OpAng->DrawCopy("hist");

      //TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
      ////gr_cosmic_E_ds->SetMarkerColor(kSpring-5);
      ////gr_cosmic_E_ds->SetMarkerStyle(20);
      ////gr_cosmic_E_ds->SetMarkerSize(0.4);
      ////gr_cosmic_E_ds->Draw("AP");

      //gr_elec_E_ds->SetTitle(""); gr_elec_E_ds->SetStats(kFALSE);
      //gr_elec_E_ds->SetMarkerColor(kViolet-5);
      //gr_elec_E_ds->SetMarkerStyle(20); gr_elec_E_ds->SetMarkerSize(0.4);
      //gr_elec_E_ds->Draw("AP");

      //gr_signal_E_ds->SetMarkerColor(kOrange+6);
      //gr_signal_E_ds->SetMarkerStyle(20); gr_signal_E_ds->SetMarkerSize(0.4);
      //gr_signal_E_ds->Draw("P same");

      //cutoff_graph->SetLineWidth(5);
      //cutoff_graph->Draw("L same");

      //// Clone for legend with larger markers
      //TGraph* leg_elec_clone = (TGraph*)gr_elec_E_ds->Clone("leg_elec_clone");
      //leg_elec_clone->SetMarkerSize(1.5);
      //
      //TGraph* leg_signal_clone = (TGraph*)gr_signal_E_ds->Clone("leg_signal_clone");
      //leg_signal_clone->SetMarkerSize(1.5);

      //TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      ////leg->AddEntry(gr_cosmic_E_ds, "Cosmics", "p");
      //leg->AddEntry(leg_elec_clone, "Photons, e+, and e-", "p");
      //leg->AddEntry(leg_signal_clone, "Signal Protons", "p");
      //leg->AddEntry(cutoff_graph, "#splitline{~90\% of electron blips}{are below this line}", "L");
      //leg->Draw();

      // Simple variable bin sizes above and below 250 MeV
      // 10 MeV bins below 250
      //for(double e = 0; e <= 250; e += 10) { binEdges.push_back(e); }
      // 50 MeV bins above 250
      //for(double e = 250; e <= 600; e += 50) { binEdges.push_back(e); } 

      // Define variable bin widths so that each bin contains 80 preselection only blips and 40 final selected blips
      // Need to refine this algorithm particulary so that the 80/40 numbers are statistically based, not arbitrary as they are now
      TCanvas *c5 = new TCanvas("c5", "Canvas", 800, 600);
      std::vector<double> binEdges;   // define variable bin edges
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
      neutronEnergy_plots->cd(); h_preFoundNeutronE = (TH1D*)h_preFoundNeutronE->Rebin(binEdges.size()-1, "", &binEdges[0]);
      outfile->cd(); TH1D* h_pre_Eff_vs_E = (TH1D*)h_preFoundNeutronE->Clone("pre_Eff_vs_E");
      h_pre_Eff_vs_E->SetTitle("Selection Signal Efficiency vs. Energy;True Neutron Energy [MeV];Selection (n, Np) Efficiency");
      h_pre_Eff_vs_E->SetLineColor(kBlack); h_pre_Eff_vs_E->SetMarkerColor(kBlack);
      neutronEnergy_plots->cd(); TH1D* h_preTrueNeutronE = (TH1D*)h_trueNeutronE->Rebin(binEdges.size()-1, "preTrueNeutronE", &binEdges[0]);
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
      outfile->cd(); TH1D* h_selection_Eff_vs_E = (TH1D*)h_foundNeutronE->Clone("selection_Eff_vs_E");
      h_selection_Eff_vs_E->SetTitle("Selection Signal Efficiency vs. Energy;True Neutron Energy [MeV];Selection (n, Np) Efficiency");
      h_selection_Eff_vs_E->SetLineColor(kRed); h_selection_Eff_vs_E->SetMarkerColor(kRed);
      neutronEnergy_plots->cd(); TH1D* h_selTrueNeutronE = (TH1D*)h_trueNeutronE->Rebin(binEdges.size()-1, "selTrueNeutronE", &binEdges[0]);
      h_selection_Eff_vs_E->Divide(h_selTrueNeutronE); //h_selection_Eff_vs_E->DrawCopy("same");
      TH1D* h_sel_drawn = (TH1D*)h_selection_Eff_vs_E->DrawCopy("same");
      
      TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg1->AddEntry(h_pre_drawn, "Preselection");
      leg1->AddEntry(h_sel_drawn, "After Cuts");
      leg1->Draw();

      TCanvas *c6 = new TCanvas("c6", "Canvas", 800, 600);
      TH1D* h_bkg_BDT_score = new TH1D("bkg_BDT_score",";BDT Score;Number of Candidates (Arb. Norm.)",100,0,1);   // removed title: BDT Score for Background Blips
      h_bkg_BDT_score->SetStats(kFALSE);
      TList* bkgList = new TList(); bkgList->Add(h_BDT_score[3]); bkgList->Add(h_BDT_score[4]);
      bkgList->Add(h_BDT_score[5]); bkgList->Add(h_BDT_score[6]); bkgList->Add(h_BDT_score[7]);
      h_bkg_BDT_score->Merge(bkgList); h_bkg_BDT_score->Scale(1.0 / h_bkg_BDT_score->Integral("width"));

      TH1D* h_bkg_drawn = (TH1D*)h_bkg_BDT_score->DrawCopy("hist");
      h_bkg_drawn->SetLineColor(kOrange);
      h_bkg_drawn->SetLineWidth(2);
      h_bkg_drawn->SetFillColorAlpha(kOrange, 0.5);
      //h_bkg_BDT_score->SetLineColor(kOrange); h_bkg_BDT_score->SetFillColorAlpha(kOrange, 0.5);
      //h_bkg_BDT_score->Draw("hist");

      TH1D* h_signal_BDT_score = new TH1D("signal_BDT_score",";BDT Score;Number of Candidates (Arb. Norm.)",100,0,1);   // removed title: BDT Score for Signal and Background Blips
      TList* sigList = new TList(); sigList->Add(h_BDT_score[1]); sigList->Add(h_BDT_score[2]);
      h_signal_BDT_score->Merge(sigList); h_signal_BDT_score->Scale(1.0 / h_signal_BDT_score->Integral("width"));

      TH1D* h_signal_drawn = (TH1D*)h_signal_BDT_score->DrawCopy("hist same");
      h_signal_drawn->SetLineColor(kBlue);
      h_signal_drawn->SetLineWidth(2);
      h_signal_drawn->SetFillColorAlpha(kBlue, 0.5);
      //h_signal_BDT_score->Draw("hist same");

      TLegend *leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg2->AddEntry(h_signal_drawn, "Signal", "lf");
      leg2->AddEntry(h_bkg_drawn, "Background", "lf");
      leg2->AddEntry((TObject*)0, "Area Normalized", "");
      leg2->Draw();

      TCanvas *c7 = new TCanvas("c7", "Canvas", 800, 600);
      h_trueNeutronEvsProtonE->SetStats(kFALSE); h_trueNeutronEvsProtonE->SetTitle("");
      h_trueNeutronEvsProtonE->DrawCopy();

      TCanvas *c8 = new TCanvas("c8", "Canvas", 800, 600);
      h_shrEnergy->DrawCopy();


      // electron lifetime function fit to signal blip efficiency vs. X position
      TCanvas *c9 = new TCanvas("c9", "Canvas", 800, 600);
      TH1D* h_signal_eff_vs_x = (TH1D*)h_signal_X->Clone("signal_eff_vs_x");
      h_signal_eff_vs_x->Sumw2(); h_trueSignalX->Sumw2();
      h_signal_eff_vs_x->Divide(h_trueSignalX);
      TF1 *fit_func = new TF1("fit_func", "[0] * TMath::Exp(-(10.0*x/1.1)/[1])", 0, 250);   // ("", A*e^(-t/T), Xmin, Xmax)
      // seed [0] from histogram: amp0 = efficiency near x~0				    // t = drift time, T = e- lifetime
      double amp0 = h_signal_eff_vs_x->GetBinContent(h_signal_eff_vs_x->FindBin(1.0));
      fit_func->SetParameters(amp0 > 0 ? amp0 : 1.0, 1000.0);
      // add limits so the minimizer can't go unphysical
      fit_func->SetParLimits(1, 10.0, 1e6);   // lifetime must be positive and sane
      fit_func->SetLineColor(kRed); fit_func->SetLineWidth(3);
      TH1D* hdrawn = (TH1D*)h_signal_eff_vs_x->DrawCopy();   // does UnZoom() implicitly to fix the stale range inherited from the pre-divide clone
      hdrawn->Fit(fit_func, "RS", "", 0, 250);
      c9->Update();


      TStyle *defaultStyle = new TStyle(*gStyle);
      gStyle->SetOptTitle(1); gStyle->SetTitleFont(42, "t");
      gStyle->SetTitleFontSize(0.05); gStyle->SetTitleBorderSize(0);
      gStyle->SetTitleFillColor(0);  gStyle->SetTitleStyle(0);
      gStyle->SetTitleX(0.0144);  gStyle->SetTitleW(0.9712);
      gStyle->SetTitleY(0.9878);  gStyle->SetTitleH(0.0574);
      
      TCanvas *c10 = new TCanvas("c10", "Canvas", 800, 600);
      h_trueNeutronEvsCosTheta->SetStats(kFALSE);
      h_trueNeutronEvsCosTheta->SetTitle("True Neutron Energy vs Cos Angle Between Neutrino Direction and Blip Position");
      h_trueNeutronEvsCosTheta->DrawCopy("COLZ");
      c10->Modified();
      c10->Update();
      gROOT->SetStyle(defaultStyle->GetName());
      //gROOT->SetStyle("Plain");//("Default");
      gROOT->ForceStyle();


      TCanvas *c11 = new TCanvas("c11", "Canvas", 800, 600);
      //h_trueSig_cosTheta->Scale(1.0 / h_trueSig_cosTheta->Integral("width"));
      h_trueSig_cosTheta->DrawCopy();
      theta_plots->cd(); TH1D* h_selSig_cosTheta = new TH1D("selSig_cosTheta",";;",100,0,1);
      TList* selSigList = new TList(); selSigList->Add(h_cos_theta[1]); sigList->Add(h_cos_theta[2]);
      h_selSig_cosTheta->Merge(selSigList); //h_selSig_cosTheta->Scale(1.0 / h_selSig_cosTheta->Integral("width"));

      TH1D* h_selSig_drawn = (TH1D*)h_selSig_cosTheta->DrawCopy("hist same");
      h_selSig_drawn->SetLineColor(kBlue);
      h_selSig_drawn->SetLineWidth(2);
      //h_selSig_drawn->SetFillColorAlpha(kBlue, 0.5);

      TLegend *leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg3->AddEntry(h_trueSig_cosTheta, "All Signal", "lf");
      leg3->AddEntry(h_selSig_drawn, "Selected Signal", "lf");
      leg3->Draw();


      TCanvas *c12 = new TCanvas("c12", "Canvas", 800, 600);     
      Normalize(h_nBlips_xp0n); Normalize(h_nBlips_xpxn); Normalize(h_nBlips_0p0n); Normalize(h_nBlips_0pxn);
      
      // --- Style ---
      h_nBlips_xp0n->SetLineColor(kBlue); h_nBlips_xpxn->SetLineColor(kRed);
      h_nBlips_0p0n->SetLineColor(kGreen+1); h_nBlips_0pxn->SetLineColor(kOrange+1);
      
      for (auto h : {h_nBlips_xp0n, h_nBlips_xpxn, h_nBlips_0p0n, h_nBlips_0pxn}) {
        h->SetLineWidth(2); h->SetFillStyle(0);  // hollow, line only
      }
      
      // Find the overall y-max to set range
      double ymax = 0;
      for (auto h : {h_nBlips_xp0n, h_nBlips_xpxn, h_nBlips_0p0n, h_nBlips_0pxn})
          ymax = std::max(ymax, h->GetMaximum());

      h_nBlips_0p0n->SetMaximum(ymax * 1.1);  // draw the tallest one first
      h_nBlips_0p0n->GetXaxis()->SetTitle("Blip Multiplicity"); h_nBlips_0p0n->GetYaxis()->SetTitle("A.U.");
      //h_nBlips_0p0n->GetXaxis()->SetRangeUser(0, 20);

      outfile->cd();
      TH1* h0p0n_drawn = (TH1*)h_nBlips_0p0n->DrawCopy("hist");
      TH1* h0pxn_drawn = (TH1*)h_nBlips_0pxn->DrawCopy("hist same");
      TH1* hxp0n_drawn = (TH1*)h_nBlips_xp0n->DrawCopy("hist same");
      TH1* hxpxn_drawn = (TH1*)h_nBlips_xpxn->DrawCopy("hist same");
      
      TLegend* leg4 = new TLegend(0.55, 0.60, 0.88, 0.88);
      leg4->AddEntry(hxp0n_drawn, "xp0n", "l");
      leg4->AddEntry(hxpxn_drawn, "xpxn", "l");
      leg4->AddEntry(h0p0n_drawn, "0p0n", "l");
      leg4->AddEntry(h0pxn_drawn, "0pxn", "l");
      leg4->SetBorderSize(0); leg4->Draw(); 


      theta_plots->cd();   // could use merge_plot3.C for this
      TCanvas *titlec13 = new TCanvas("c13", "Canvas", 800, 600);
      h_ProtonCosTheta->Sumw2(); h_NeutronCosTheta->Sumw2(); h_PionCosTheta->Sumw2(); //h_OtherCosTheta->Sumw2(); 
      Normalize(h_ProtonCosTheta); Normalize(h_NeutronCosTheta); Normalize(h_PionCosTheta); //Normalize(h_OtherCosTheta);

      h_ProtonCosTheta->SetLineColor(kBlue); h_NeutronCosTheta->SetLineColor(kRed);
      h_PionCosTheta->SetLineColor(kGreen+1); //h_OtherCosTheta->SetLineColor(kOrange+1);
      h_NeutronCosTheta->GetXaxis()->SetTitle("Cos(#Theta)"); h_NeutronCosTheta->GetYaxis()->SetTitle("A.U.");
      h_NeutronCosTheta->SetTitle("Cos(Angle) Between Neutrino Direction and Vertex to Secondary (n,p) Blip Position With Various Primary Ancestors");

      TH1* h_ProtonCosTheta_drawn = (TH1*)h_ProtonCosTheta->DrawCopy("E1 hist");
      TH1* h_NeutronCosTheta_drawn = (TH1*)h_NeutronCosTheta->DrawCopy("E1 hist same");
      TH1* h_PionCosTheta_drawn = (TH1*)h_PionCosTheta->DrawCopy("E1 hist same");
      //TH1* h_OtherCosTheta_drawn = (TH1*)h_OtherCosTheta->DrawCopy("hist same");

      TLegend* leg5 = new TLegend(0.55, 0.60, 0.88, 0.88);
      leg5->AddEntry(h_ProtonCosTheta_drawn, "Proton", "le");
      leg5->AddEntry(h_NeutronCosTheta_drawn, "Neutron", "le");
      leg5->AddEntry(h_PionCosTheta_drawn, "#pi+/-", "le");
      //leg5->AddEntry(h_OtherCosTheta_drawn, "Other", "l");
      leg5->SetBorderSize(0); leg5->Draw();

    }

    //Float_t a;
    //TTree otree("evd_tree", "TTree with a single float branch for total evd");
    //TBranch* branch = otree.Branch("total_event_displays", &a, "total_event_displays/F");
    //a = selectedEvents;
    //otree.Fill();

    //===========================================
    // Write the output and close the file
    outfile->Write();
    outfile->Close();

  }
