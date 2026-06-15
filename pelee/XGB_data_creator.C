
  #include "TROOT.h"
  #include <cmath>
  #include <algorithm>
  #include <iomanip>   // for std::boolalpha

  int maxEvents = -2000;

  int			   n_pfps;
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
  float 		   slice_orig_topo_score = 0;
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
  std::vector<int>         * blip_trkid;

  std::vector<float>	   * trk_sce_start_x_v = 0;
  std::vector<float>	   * trk_sce_start_y_v = 0;
  std::vector<float>	   * trk_sce_start_z_v = 0;
  std::vector<float>	   * trk_sce_end_x_v = 0;
  std::vector<float>	   * trk_sce_end_y_v = 0;
  std::vector<float>	   * trk_sce_end_z_v = 0;
  std::vector<float>	   * trk_llr_pid_score_v = 0;
  //std::vector<float>	   * trk_len = 0;
  std::vector<float>       * trk_score_v = 0;

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
  std::vector<float>	   * pfng2hipfrac;
  std::vector<int>	   * trk_pfp_id_v;
  //std::vector<int>         * trk_id;   // id of leading trk only

  std::vector<int>	   * all_mc_trkid = 0;
  std::vector<int>	   * all_mc_pdg = 0;
  std::vector<float>	   * all_mc_E = 0;
  std::vector<float>	   * all_mc_endx = 0;
  std::vector<float>	   * all_mc_endy = 0;
  std::vector<float>	   * all_mc_endz = 0;
  std::vector<int>	   * all_mc_mother = 0;
  std::vector<std::string> * all_mc_process = 0;

  //Int_t    reco_id[1000];   // these are ints not vectors!

  //bool   nu_isNeutrino;
  float  reco_nu_vtx_x;
  float  reco_nu_vtx_y;
  float  reco_nu_vtx_z;

  // Input file
  //std::string fFileName = "checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root";   // retupled BNB nu
  //std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260224.root";   // BNB nu with more blip vars
  std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root";

  std::string nuTreeName = "nuselection/NeutrinoSelectionFilter";
  //std::string wcTreeName = "wcpselection/T_PFeval";

  std::unordered_map<int, int> map_allTrkID_partIDX;
  std::unordered_map<int, int> map_allTrkID_momTrkID;

  //====================================================================
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

  bool IsInFV_Exact( TVector3 point ) {
    return point.X() >= 0.	&&  point.X() <= 256.  &&
	   point.Y() >= -116.25 &&  point.Y() <= 116.25   &&
	   point.Z() >= 0.	&&  point.Z() <= 1036.8 ;
  }   // from R. Acciarri et al 2017 JINST 12 P02017

  bool IsInFV_Buffer( const TVector3& point ) {
    return point.X() >= 10.     &&  point.X() <= 246.35 &&
           point.Y() >= -106.5  &&  point.Y() <= 106.5  &&
           point.Z() >= 10.     &&  point.Z() <= 1026.8 ;
  }   // Same as Burk and Andy's analysis but may want to loosen since our analysis isn't track focused

  bool IsInFV_BufferAndDeadZone( const TVector3& point ) {
    return point.X() >= 5.      &&  point.X() <= 251.   &&
           point.Y() >= -110.   &&  point.Y() <= 110    &&
           point.Z() >= 20.     &&  point.Z() <= 986.   &&
           (point.Z() < 675.    ||  point.Z() > 775.)   ;
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

  std::vector<float>       blip_E;
  std::vector<float>       dist2Vtx;
  std::vector<float>	   dist2Trk;
  std::vector<float>	   dist2Edge;
  std::vector<float>       blip_size;
  std::vector<float>       blip_density;
  std::vector<float>	   alignment_angle;
  std::vector<float>	   Z_dist2Vtx;
  //std::vector<float>	   dist2ShrVtx;
  //std::vector<float>	   angle_blip_shr;
  //std::vector<float>       blip_X;
  //std::vector<float>       blip_Y;
  //std::vector<float>       blip_Z;
  //std::vector<float>       blip_dX;
  //std::vector<float>       blip_dW;
  //std::vector<float>       blip_idx;
  //std::vector<float>	     event;

  // Helper function to set up branches
  void SetupBranches(TTree* tree) {
    tree->Branch("blip_E",	    &blip_E);
    tree->Branch("dist2Vtx",	    &dist2Vtx);
    tree->Branch("dist2Trk",        &dist2Trk);
    tree->Branch("dist2Edge",       &dist2Edge);
    tree->Branch("blip_size",	    &blip_size);
    tree->Branch("blip_density",    &blip_density);
    tree->Branch("alignment_angle", &alignment_angle);
    tree->Branch("Z_dist2Vtx",	    &Z_dist2Vtx);
    //tree->Branch("dist2ShrVtx",	    &dist2ShrVtx);
    //tree->Branch("angle_blip_shr",  &angle_blip_shr);
  }

  // ********************************************************
  void XGB_data_creator(){

    // Read in the TFile1000 and grab the TTree
    TFile* file = new TFile(("files/"+fFileName).c_str(),"READ");
    TTree* nuTree = (TTree*)file->Get(nuTreeName.c_str());
    //TTree* wcTree = (TTree*)file->Get(wcTreeName.c_str());
    std::cout<<"Reading data file "<<fFileName<<"\n\n";

    // Set branch address for all variables you plan to use
    nuTree->SetBranchAddress("n_pfps",              &n_pfps);
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
    //nuTree->SetBranchAddress("blip_true_energy",    &blip_true_energy);

    //nuTree->SetBranchAddress("blip_pl0_nwires",	    &blip_pl0_nwires);
    //nuTree->SetBranchAddress("blip_pl1_nwires",	    &blip_pl1_nwires);
    //nuTree->SetBranchAddress("blip_pl2_nwires",	    &blip_pl2_nwires);
    nuTree->SetBranchAddress("blip_pl0_bydeadwire", &blip_pl0_bydeadwire);
    nuTree->SetBranchAddress("blip_pl1_bydeadwire", &blip_pl1_bydeadwire);
    nuTree->SetBranchAddress("blip_pl2_bydeadwire", &blip_pl2_bydeadwire);
    nuTree->SetBranchAddress("blip_rnn_dir_isValid",&blip_rnn_dir_isValid);
    nuTree->SetBranchAddress("blip_rnn_dir_x",      &blip_rnn_dir_x);
    nuTree->SetBranchAddress("blip_rnn_dir_y",      &blip_rnn_dir_y);
    nuTree->SetBranchAddress("blip_rnn_dir_z",      &blip_rnn_dir_z);
    //nuTree->SetBranchAddress("blip_trkid",	    &blip_trkid);

    nuTree->SetBranchAddress("trk_sce_start_x_v",   &trk_sce_start_x_v);
    nuTree->SetBranchAddress("trk_sce_start_y_v",   &trk_sce_start_y_v);
    nuTree->SetBranchAddress("trk_sce_start_z_v",   &trk_sce_start_z_v);
    nuTree->SetBranchAddress("trk_sce_end_x_v",     &trk_sce_end_x_v);
    nuTree->SetBranchAddress("trk_sce_end_y_v",     &trk_sce_end_y_v);
    nuTree->SetBranchAddress("trk_sce_end_z_v",     &trk_sce_end_z_v);
    nuTree->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v);
    ////nuTree->SetBranchAddress("trk_len",		    &trk_len);
    nuTree->SetBranchAddress("trk_score_v",	    &trk_score_v);

    //nuTree->SetBranchAddress("n_showers",           &n_showers);
    nuTree->SetBranchAddress("shr_start_x_v",	    &shr_start_x_v);
    nuTree->SetBranchAddress("shr_start_y_v",	    &shr_start_y_v);
    nuTree->SetBranchAddress("shr_start_z_v",	    &shr_start_z_v);
    nuTree->SetBranchAddress("shr_px_v",	    &shr_px_v);
    nuTree->SetBranchAddress("shr_py_v",	    &shr_py_v);
    nuTree->SetBranchAddress("shr_pz_v",	    &shr_pz_v);
    //nuTree->SetBranchAddress("shr_dist_v",	    &shr_dist_v);
    //nuTree->SetBranchAddress("shr_openangle_v",     &shr_openangle_v);
    //nuTree->SetBranchAddress("shr_openangle",       &shr_openangle);
    nuTree->SetBranchAddress("shr_energy_y_v",	    &shr_energy_y_v);
    //nuTree->SetBranchAddress("pfng2hipfrac",	    &pfng2hipfrac);
    //nuTree->SetBranchAddress("trk_pfp_id_v",	    &trk_pfp_id_v);

    nuTree->SetBranchAddress("all_mc_trkid",	    &all_mc_trkid);
    nuTree->SetBranchAddress("all_mc_pdg",	    &all_mc_pdg);
    //nuTree->SetBranchAddress("all_mc_E",	    &all_mc_E);
    //nuTree->SetBranchAddress("all_mc_endx",         &all_mc_endx);
    //nuTree->SetBranchAddress("all_mc_endy",         &all_mc_endy);
    //nuTree->SetBranchAddress("all_mc_endz",         &all_mc_endz);
    nuTree->SetBranchAddress("all_mc_mother",       &all_mc_mother);
    //nuTree->SetBranchAddress("all_mc_process",      &all_mc_process);

    nuTree->SetBranchAddress("nu_pdg",		    &nu_pdg);
    nuTree->SetBranchAddress("ccnc",		    &ccnc);
    nuTree->SetBranchAddress("reco_nu_vtx_sce_x",   &reco_nu_vtx_sce_x);
    nuTree->SetBranchAddress("reco_nu_vtx_sce_y",   &reco_nu_vtx_sce_y);
    nuTree->SetBranchAddress("reco_nu_vtx_sce_z",   &reco_nu_vtx_sce_z);
    nuTree->SetBranchAddress("slice_orig_pass_id",  &slice_orig_pass_id);
    //nuTree->SetBranchAddress("isVtxInFiducial",     &isVtxInFiducial);
    nuTree->SetBranchAddress("slice_orig_topo_score",&slice_orig_topo_score);
    nuTree->SetBranchAddress("crtveto",		    &crtveto);
    nuTree->SetBranchAddress("_closestNuCosmicDist",&_closestNuCosmicDist);
    //nuTree->SetBranchAddress("bdt_numuCC",	    &bdt_numuCC);
    //nuTree->SetBranchAddress("reco_nu_vtx_x",       &reco_nu_vtx_x);
    //nuTree->SetBranchAddress("reco_nu_vtx_y",       &reco_nu_vtx_y);
    //nuTree->SetBranchAddress("reco_nu_vtx_z",       &reco_nu_vtx_z);

    //wcTree->SetBranchAddress("reco_id",		    &reco_id);

    // make output file and set its trees & branches
    TFile* BDT_data = new TFile("BDT_data.root", "recreate");
    TTree* sig_train = new TTree("sig_train","sig_train");
    TTree* bkg_train = new TTree("bkg_train","bkg_train");
    TTree* sig_test = new TTree("sig_test","sig_test");
    TTree* bkg_test = new TTree("bkg_test","bkg_test");

    SetupBranches(sig_train); SetupBranches(bkg_train);
    SetupBranches(sig_test); SetupBranches(bkg_test);

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

    int   category;

    //============================== EVENT LOOP ==============================//
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      //if( iEvent > 100000 ) continue;
      nuTree->GetEntry(iEvent);

      int all_mc_length = static_cast<int>(all_mc_trkid->size());

      if(iEvent < 10 || iEvent % 10000 == 0) {
	std::cout<<"\n=======================================\n";
	std::cout<<"Reading event "<<iEvent<<" out of "<<numEvents<<"\n";
	std::cout<<"NBlips: "<<nblips_saved<<"   NParticles (all_mc): "<<all_mc_length<<"\n";
      }
      selectedEvents++;

      map_allTrkID_partIDX.clear();
      map_allTrkID_momTrkID.clear();

      TVector3 reco_nu_vtx( reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z );

      //============================== NUMUCC INCLUSIVE SELECTION ==============================//
      // Preselection
                       if( slice_orig_pass_id != 1 ||   // replaces nslice != 1
      IsInFV_BufferAndDeadZone( reco_nu_vtx ) != 1 ||
                     slice_orig_topo_score <= 0.06 ||   //replaces topological_score
                                      crtveto == 1 ||
                         _closestNuCosmicDist <= 5 ) { continue; }   // not SCE in ubana code

      //============================== PFP LOOP ==============================//
      std::vector<std::tuple<TVector3, TVector3, float>> shr_vtx_P_E; bool numuCC = false;
      for(size_t i_pfp = 0; i_pfp < n_pfps; i_pfp++){

	// Muon Tracks   (could we assume this is always the leading track?)
	if( trk_score_v->at(i_pfp) > 0.8 ) {   // 0.45 if 0.8 checked elsewhere
	  TVector3 trkStart(trk_sce_start_x_v->at(i_pfp), trk_sce_start_y_v->at(i_pfp), trk_sce_start_z_v->at(i_pfp));
	  TVector3 trkEnd(trk_sce_end_x_v->at(i_pfp), trk_sce_end_y_v->at(i_pfp), trk_sce_end_z_v->at(i_pfp));
	  float trk_length = (trkEnd - trkStart).Mag();

	  if( IsInFV_Buffer(trkStart) && IsInFV_Buffer(trkEnd) &&
	      std::min( (trkStart-reco_nu_vtx).Mag(), (trkEnd-reco_nu_vtx).Mag() ) < 4 &&
	      trk_llr_pid_score_v->at(i_pfp) > 0.2 &&
	      trk_length > 10 )
	  { numuCC = true; }
	}

	//pfng2hipfrac
	//std::cout<<"trk_pfp_id_v: "<<trk_pfp_id_v->at(i_pfp)<<"\n\n";   // sometime large negative (=invalid, data?)
	// Showers
	if( trk_score_v->at(i_pfp) < 0.45 ) {
	  TVector3 shrStart(shr_start_x_v->at(i_pfp), shr_start_y_v->at(i_pfp), shr_start_z_v->at(i_pfp));
	  TVector3 shr_p(shr_px_v->at(i_pfp), shr_py_v->at(i_pfp), shr_pz_v->at(i_pfp));   // already normalized
	  //TVector3 shr_distVec = shr_p * shr_dist_v->at(i_pfp);   // ( |Px|*d, |Py|*d, |Pz|*d )
	  //TVector3 shrEnd = shrStart + shr_distVec;
	  float shrE = shr_energy_y_v->at(i_pfp);
	  shr_vtx_P_E.emplace_back(shrStart, shr_p, shrE);   //push_back(std::make_tuple(shrStart, shr_p, shrE));
	}
      }

      // true numuCC but doesn't pass reco. selection since we use these in our analysis
      // true numuCC has ccnc == 0 (CC) and nu_pdg == 14 (muon neutrino)
      if( numuCC || ccnc != 0 || nu_pdg != 14 ) continue;

      //============================== PARTICLE LOOP ==============================//
      for(int i=0; i<all_mc_length; i++){

	int mc_trkID = all_mc_trkid->at(i); int mc_momTrkID = all_mc_mother->at(i);

	map_allTrkID_partIDX[mc_trkID] = i;
	map_allTrkID_momTrkID[mc_trkID] = mc_momTrkID;

      }//end particle loop

      std::vector<float> sig_blip_E;
      std::vector<float> sig_dist2Vtx;
      std::vector<float> sig_dist2Trk;
      std::vector<float> sig_dist2Edge;
      std::vector<float> sig_blip_size;
      std::vector<float> sig_blip_density;
      std::vector<float> sig_alignment_angle;
      std::vector<float> sig_Z_dist2Vtx;
      //std::vector<float> sig_blip_to_shrVtx;
      //std::vector<float> sig_blip_shr_angle;
      //std::vector<float> sig_blip_X;
      //std::vector<float> sig_blip_Y;
      //std::vector<float> sig_blip_Z;
      //std::vector<float> sig_blip_dX;
      //std::vector<float> sig_blip_dW;
      //std::vector<float> sig_blip_idx;
      //std::vector<float> sig_iEvent;

      std::vector<float> bkg_blip_E; 
      std::vector<float> bkg_dist2Vtx;
      std::vector<float> bkg_dist2Trk;
      std::vector<float> bkg_dist2Edge;
      std::vector<float> bkg_blip_size;
      std::vector<float> bkg_blip_density;
      std::vector<float> bkg_alignment_angle;
      std::vector<float> bkg_Z_dist2Vtx;
      //std::vector<float> bkg_blip_to_shrVtx;
      //std::vector<float> bkg_blip_shr_angle;

      //============================== BLIP LOOP ==============================//
      for(size_t i=0; i<nblips_saved; i++){

        totalBlips++; category=6;   // default to "other" if not assigned

	float blipX = blip_x->at(i), blipY = blip_y->at(i), blipZ = blip_z->at(i);
        TVector3 blip_loc( blipX, blipY, blipZ ); TVector3 reco_trace_vec( blip_loc - reco_nu_vtx );
	float dist=(blip_loc-reco_nu_vtx).Mag(); float Zdist = blipZ - reco_nu_vtx.Z();

	float blipE = blip_energy->at(i); totalBlipE+=blipE;

	float dx = blip_dx->at(i); float dw = blip_dw->at(i);
	float blipSize = sqrt( dx*dx + dw*dw );   // direct multiplication is fastest ig

	float density = blipE/blipSize; float dist_to_trk = blip_proxtrkdist->at(i);
	float cutoff = cutoff_graph->Eval(blipE);

	int pdg = blip_true_pdg->at(i);	int g4id = blip_true_g4id->at(i);	

	int allIDX = -9; bool isPrimary; 
	int mom_trkID = -9; bool mom_isPrimary = false;
	int mom_allIDX = -9; int mom_pdg = -9;
	auto it = map_allTrkID_partIDX.find(g4id);
	if( it != map_allTrkID_partIDX.end() ) {
	  allIDX = it->second;
	  mom_trkID = all_mc_mother->at( allIDX );
	  isPrimary = mom_trkID == 0;   // if mother == 0, particle is primary
 
	  if( !isPrimary && map_allTrkID_partIDX.count(mom_trkID) != 0 ) {
	    mom_allIDX = map_allTrkID_partIDX[mom_trkID];
	    mom_isPrimary = all_mc_mother->at( mom_allIDX ) == 0; mom_pdg = all_mc_pdg->at( mom_allIDX );
	  }
	}

        //============================== CUTS ==============================//
	if( blip_touchtrk->at(i) || dist_to_trk <= 1 || !IsInFV_Exact(blip_loc) ) continue;
	if( blip_pl0_bydeadwire->at(i) == 1 || blip_pl1_bydeadwire->at(i) == 1 || blip_pl2_bydeadwire->at(i) == 1 ) continue;
        if( dist >= 100 ) continue; //if( density < Fit_Function_PID_Cut(blipE) ) continue; if( blipE <= 3 ) continue;
        //bool isInShower = false; 
        //for( const auto& [shr_vtx, shr_P, E] : shr_vtx_P_E ) { if ( IsWithinShowerCone(shr_vtx, shr_P, blip_loc, 90.0) ) { isInShower = true; break; } }
        //if( isInShower ) continue;
	//double minDist2Shr = 100000; TVector3 blip_to_shrVtx(0,0,0); double blip_shr_angle = 0;
	//for( const auto& [shr_vtx, shr_P, E] : shr_vtx_P_E ) {
	//   if ( IsWithinShowerCone(shr_vtx, shr_P, blip_loc, 90.0) && (blip_loc - shr_vtx).Mag() < minDist2Shr ) {
	//     blip_to_shrVtx = blip_loc - shr_vtx;
	//     blip_shr_angle = blip_to_shrVtx.Angle(shr_P);
	//     minDist2Shr = blip_to_shrVtx.Mag();
	//   }
	//}	
	//if( blip_to_shrVtx.Mag() > 10000 || blip_shr_angle > 10000 )
	//  { std:cout<<"What the helly?!\nMinimum distance to shower: "<<minDist2Shr<<"\nblip2ShrVtx variable: "<<blip_to_shrVtx.Mag()<<"\nangle_blip_shr variable: "<<blip_shr_angle<<"\n\n"; }
	//std::cout<<"Blip trk ID: "<<blip_trkid->at(i)<<"\n\n";

        TVector3 reco_dir(0,0,0); float align_angle=-9;
        if( !blip_rnn_dir_isValid->at(i) ) continue;
	else {
          // RNN model prediction
          reco_dir.SetXYZ( blip_rnn_dir_x->at(i), blip_rnn_dir_y->at(i), blip_rnn_dir_z->at(i) );
          align_angle = reco_dir.Angle(reco_trace_vec);
	  //align_angle = align_angle*TMath::RadToDeg(); if( align_angle > 90 ) align_angle = 180-align_angle; 
        }

	//std::cout<<"Blip predicted direction alignment to vertex (angle): "<<align_angle<<"\n";
	//std::cout<<"Blip location: ("<<blipX<<", "<<blipY<<", "<<blipZ<<")\nDistance to nearest edge: "<<distToNearestEdge(blip_loc)<<"\n";

        // Fill signal and background variables
	if( !isPrimary && pdg == 2212 && mom_pdg == 2112 && mom_isPrimary ) {   // non-primary protons from primary neutron mom
	  sig_blip_E.push_back(blipE);
	  sig_dist2Vtx.push_back(dist);
	  sig_dist2Trk.push_back(dist_to_trk);
	  sig_dist2Edge.push_back(distToNearestEdge(blip_loc));
	  sig_blip_size.push_back(blipSize);
	  sig_blip_density.push_back(density);
	  sig_alignment_angle.push_back(align_angle);
	  sig_Z_dist2Vtx.push_back(Zdist);
	  //sig_blip_to_shrVtx.push_back(blip_to_shrVtx.Mag());
	  //sig_blip_shr_angle.push_back(blip_shr_angle);
	  //sig_blip_X.push_back(blipX);
          //sig_blip_Y.push_back(blipY);
          //sig_blip_Z.push_back(blipZ);
	  //sig_blip_dX.push_back(dx);
	  //sig_blip_dW.push_back(dw);
	  //sig_blip_idx.push_back(i);
	  //sig_iEvent.push_back(iEvent);

	} else {
	  bkg_blip_E.push_back(blipE);
	  bkg_dist2Vtx.push_back(dist);
	  bkg_dist2Trk.push_back(dist_to_trk);
	  bkg_dist2Edge.push_back(distToNearestEdge(blip_loc));
          bkg_blip_size.push_back(blipSize);
          bkg_blip_density.push_back(density);
	  bkg_alignment_angle.push_back(align_angle);
	  bkg_Z_dist2Vtx.push_back(Zdist);
	  //bkg_blip_to_shrVtx.push_back(blip_to_shrVtx.Mag());
	  //bkg_blip_shr_angle.push_back(blip_shr_angle);
	}

      }//end blip loop

      // Fill trees with blip variable arrays
      if( iEvent % 2 == 0 ) {

	blip_E = sig_blip_E;
	dist2Vtx = sig_dist2Vtx;
	dist2Trk = sig_dist2Trk;
	dist2Edge = sig_dist2Edge;
	blip_size = sig_blip_size;
	blip_density = sig_blip_density;
	alignment_angle = sig_alignment_angle;
	Z_dist2Vtx = sig_Z_dist2Vtx;
	//dist2ShrVtx = sig_blip_to_shrVtx;
	//angle_blip_shr = sig_blip_shr_angle;
        //blip_X = sig_blip_X;
        //blip_Y = sig_blip_Y;
        //blip_Z = sig_blip_Z;
	//blip_dX = sig_blip_dX;
	//blip_dW = sig_blip_dW;
	//blip_idx = sig_blip_idx;
	//event = sig_iEvent;
	sig_train->Fill();

	blip_E = bkg_blip_E;
	dist2Vtx = bkg_dist2Vtx;
	dist2Trk = bkg_dist2Trk;
	dist2Edge = bkg_dist2Edge;
	blip_size = bkg_blip_size;
	blip_density = bkg_blip_density;
	alignment_angle = bkg_alignment_angle;
	Z_dist2Vtx = bkg_Z_dist2Vtx;
	//dist2ShrVtx = bkg_blip_to_shrVtx;
	//angle_blip_shr = bkg_blip_shr_angle;
	bkg_train->Fill();

      } else {

	blip_E = sig_blip_E;
	dist2Vtx = sig_dist2Vtx;
	dist2Trk = sig_dist2Trk;
	dist2Edge = sig_dist2Edge;
	blip_size = sig_blip_size;
	blip_density = sig_blip_density;
	alignment_angle = sig_alignment_angle;
	Z_dist2Vtx = sig_Z_dist2Vtx;
	//dist2ShrVtx = sig_blip_to_shrVtx;
	//angle_blip_shr = sig_blip_shr_angle;
	sig_test->Fill();

	blip_E = bkg_blip_E;
	dist2Vtx = bkg_dist2Vtx;
	dist2Trk = bkg_dist2Trk;
	dist2Edge = bkg_dist2Edge;
	blip_size = bkg_blip_size;
	blip_density = bkg_blip_density;
	alignment_angle = bkg_alignment_angle;
	Z_dist2Vtx = bkg_Z_dist2Vtx;
	//dist2ShrVtx = bkg_blip_to_shrVtx;
	//angle_blip_shr = bkg_blip_shr_angle;
	bkg_test->Fill();

     }

	//if( sig_blip_E.count() < 1 ) std::cout<<"sig_blip_E is empty for this event:\nPDG = "<<pdg<<"\n\n";
        //if( sig_blip_E.count() < 1 ) std::cout<<"sig_blip_E is empty for this event:\nPDG = "<<pdg<<"\n\n";

    }//end events loop

    // could also try it this way:
    //  blip_E = sig_blipE;
    //  dist2Vtx = sig_dist2Vtx;
    //  sig_train->Fill();
    //  
    //  blip_E = bkg_blipE;
    //  dist2Vtx = bkg_dist2Vtx;
    //  bkg_train->Fill();

    //} else {

    //  blip_E = sig_blipE;
    //  dist2Vtx = sig_dist2Vtx;
    //  sig_test->Fill();
    //  
    //  blip_E = bkg_blipE;
    //  dist2Vtx = bkg_dist2Vtx; 
    //  bkg_test->Fill();

    // }
  
    //===========================================
    // Write the output and close the file
    BDT_data->Write("", TObject::kOverwrite);   //BDT_data->Write();
    BDT_data->Close();

  }
