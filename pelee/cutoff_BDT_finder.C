
  #include "TROOT.h"
  #include <cmath>
  #include <algorithm>
  #include <iomanip>   // for std::boolalpha
  #include <iostream>
  #include <vector>
  #include <string>
  #include "merge_plot2.C"

  // ROOT / TMVA headers
  #include "TMVA/RReader.hxx"
  #include "TMVA/Reader.h"
  #include "ROOT/RVec.hxx"

  int maxEvents = -2000;

  int			   n_pfps;
  //int			   selected;
  //int			   n_tracks;
  //int			   n_showers;
  float			   reco_nu_vtx_sce_x = 0;
  float			   reco_nu_vtx_sce_y = 0;
  float			   reco_nu_vtx_sce_z = 0;
  int			   slice_orig_pass_id = 0;
  bool			   isVtxInFiducial = 0;
  float			   slice_orig_topo_score=0;
  //float 		   topological_score = 0;
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

  //Int_t    reco_id[1000];   // these are ints not vectors!

  //bool   nu_isNeutrino;
  float  reco_nu_vtx_x;
  float  reco_nu_vtx_y;
  float  reco_nu_vtx_z;

  // Input file
  //std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260224.root";   // BNB nu with more blip vars (full events)
  std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root";

  std::string nuTreeName = "nuselection/NeutrinoSelectionFilter";

  std::unordered_map<int, int> map_allTrkID_partIDX;
  std::unordered_map<int, int> map_allTrkID_momTrkID;

  bool IsInFV_Exact( const TVector3& point ) {
    return point.X() >= 0.      &&  point.X() <= 256.   &&
           point.Y() >= -116.25 &&  point.Y() <= 116.25 &&
           point.Z() >= 0.      &&  point.Z() <= 1036.8 ;
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

  bool IsWithinShowerCone(const TVector3& vtx,      // Shower Vtx
			  const TVector3& dir,      // Shower direction (momentum)
			  const TVector3& blip,     // Blip coordinates
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
  //TGraph* cutoff_graph = (TGraph*)f->Get("cutoff_graph");
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

  // ── Load the model ────────────────────────────────────────────────────────
  //TMVA::Experimental::RReader model("n_p_select_noShow_100k.xml");
  TMVA::Reader *reader = new TMVA::Reader( );

  // Provide the features exactly as used during training
  //std::vector<std::string> featureNames = {"blip_E", "dist2Vtx", "blip_X", "blip_Y", "blip_Z", "blip_dX", "blip_dW", "blip_size", "blip_density"};
  float blip_E, dist2Vtx, dist2Trk, dist2Edge, blip_size, blip_density, alignment_angle, Z_dist2Vtx;//, dist2ShrVtx, angle_blip_shr; 
  //std::vector<float> inputFeatures = { blip_E, dist2Vtx, blip_X, blip_Y, blip_Z, blip_dX, blip_dW, blip_size, blip_density };
  reader->AddVariable("blip_E", &blip_E);
  reader->AddVariable("dist2Vtx", &dist2Vtx);
  reader->AddVariable("dist2Trk", &dist2Trk);
  reader->AddVariable("dist2Edge", &dist2Edge);
  reader->AddVariable("blip_size", &blip_size);
  reader->AddVariable("blip_density", &blip_density);
  reader->AddVariable("alignment_angle", &alignment_angle);
  reader->AddVariable("Z_dist2Vtx", &Z_dist2Vtx);
  //reader->AddVariable("dist2ShrVtx", &dist2ShrVtx);
  //reader->AddVariable("angle_blip_shr", &angle_blip_shr);

  reader->BookMVA("BDT", "models_BDT/n_p_select_truNuZdist.xml");

  // ********************************************************
  void cutoff_BDT_finder(){
    // Read in the TFile and grab the TTree
    TFile* file = new TFile(("files/"+fFileName).c_str(),"READ");
    TTree* nuTree = (TTree*)file->Get(nuTreeName.c_str());
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
    //nuTree->SetBranchAddress("blip_true_energy",    &blip_true_energy);

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
    //nuTree->SetBranchAddress("blip_true_dir_x",	    &blip_true_dir_x);
    //nuTree->SetBranchAddress("blip_true_dir_y",	    &blip_true_dir_y);
    //nuTree->SetBranchAddress("blip_true_dir_z",	    &blip_true_dir_z);

    nuTree->SetBranchAddress("n_pfps",		    &n_pfps);
    //nuTree->SetBranchAddress("pfng2semlabel",	    &pfng2semlabel);

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
    //nuTree->SetBranchAddress("shr_dist_v",	    &shr_dist_v);
    //nuTree->SetBranchAddress("shr_openangle_v",     &shr_openangle_v);
    //nuTree->SetBranchAddress("shr_openangle",       &shr_openangle);
    nuTree->SetBranchAddress("shr_energy_y_v",	    &shr_energy_y_v);

    nuTree->SetBranchAddress("all_mc_trkid",	    &all_mc_trkid);
    nuTree->SetBranchAddress("all_mc_pdg",	    &all_mc_pdg);
    nuTree->SetBranchAddress("all_mc_E",	    &all_mc_E);
    //nuTree->SetBranchAddress("all_mc_endx",         &all_mc_endx);
    //nuTree->SetBranchAddress("all_mc_endy",         &all_mc_endy);
    //nuTree->SetBranchAddress("all_mc_endz",         &all_mc_endz);
    nuTree->SetBranchAddress("all_mc_mother",       &all_mc_mother);
    //nuTree->SetBranchAddress("all_mc_process",      &all_mc_process);

    nuTree->SetBranchAddress("reco_nu_vtx_sce_x",   &reco_nu_vtx_sce_x);
    nuTree->SetBranchAddress("reco_nu_vtx_sce_y",   &reco_nu_vtx_sce_y);
    nuTree->SetBranchAddress("reco_nu_vtx_sce_z",   &reco_nu_vtx_sce_z);
    nuTree->SetBranchAddress("slice_orig_pass_id",  &slice_orig_pass_id);
    nuTree->SetBranchAddress("slice_orig_topo_score",&slice_orig_topo_score);
    //nuTree->SetBranchAddress("isVtxInFiducial",     &isVtxInFiducial);
    //nuTree->SetBranchAddress("topological_score",   &topological_score);
    nuTree->SetBranchAddress("crtveto",		    &crtveto);
    nuTree->SetBranchAddress("_closestNuCosmicDist",&_closestNuCosmicDist);

    int numEvents = nuTree->GetEntries();
    std::cout<<"There are "<<numEvents<<" events in the nu tree.\n";  

    TH1D* h_signal_BDT_score = new TH1D("signal_BDT_score","BDT Score for Signal and Background Blips;Score;Arbitrary Units (Area Normalized)",60,-10,2);
    TH1D* h_bkg_BDT_score = new TH1D("bkg_BDT_score","BDT Score for Signal and Background Blips;Score;Arbitrary Units (Area Normalized)",60,-10,2);
    TGraphErrors* g_rareSearch_vs_BDT_score = new TGraphErrors();
    TGraphErrors* g_ExP_sq_vs_BDT_score = new TGraphErrors();
    TGraphErrors* g_pur_vs_BDT_score = new TGraphErrors();
    TGraphErrors* g_eff_vs_BDT_score = new TGraphErrors();
    TGraphErrors* g_ExP_vs_BDT_score = new TGraphErrors();
    //TGraph* g_effStatErr_vs_BDT_score = new TGraph();
    //g_rareSearch_vs_BDT_score->SetName("g_rareSearch_vs_BDT_score"); g_rareSearch_vs_BDT_score->SetTitle("Figure of Merit vs BDT Score;Score;FOM");
    //g_pur_vs_BDT_score->SetName("g_pur_vs_BDT_score");
    //g_eff_vs_BDT_score->SetName("g_eff_vs_BDT_score");
    //g_ExP_vs_BDT_score->SetName("g_ExP_vs_BDT_score");

    std::unordered_map<int, std::vector<int>> map_n_p;
    std::vector<double> sig_scores, bkg_scores;
    std::vector<double> neutronID_sig_scores, totalSelected_scores;

    int total_nu_Ns=0;

    //============================== EVENT LOOP ==============================//
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      //if( iEvent > 100000 ) continue;
      //if( iEvent < 100000  && iEvent % 2 == 0 ) continue;
      nuTree->GetEntry(iEvent);

      int all_mc_length = static_cast<int>(all_mc_trkid->size());
      std::set<int> tagged_neutron_trkIDs;

      if(iEvent < 10 || iEvent % 10000 == 0) {
	std::cout<<"\n=======================================\n";
	std::cout<<"Reading event "<<iEvent<<" out of "<<numEvents<<"\n";
	std::cout<<"NBlips: "<<nblips_saved<<"   NParticles (all_mc): "<<all_mc_length<<"\n";
      }

      map_allTrkID_partIDX.clear();
      map_allTrkID_momTrkID.clear();
      map_n_p.clear();

      TVector3 reco_nu_vtx( reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z );

      //============================== NUMUCC INCLUSIVE SELECTION ==============================//
      // Preselection
		   if( slice_orig_pass_id != 1 ||   // replaces nslice != 1
  IsInFV_BufferAndDeadZone( reco_nu_vtx ) != 1 ||
		 slice_orig_topo_score <= 0.06 ||   //replaces topological_score
				  crtveto == 1 ||
		     _closestNuCosmicDist <= 5 ) { continue; }

      //============================== PFP LOOP ==============================//
      int n_shrs = 0;
      bool numuCC = false;
      std::vector<std::tuple<TVector3, TVector3, float>> shr_vtx_P_E;
      for(size_t i_pfp = 0; i_pfp < n_pfps; i_pfp++){

	// Muon Track   (could we assume this is always the leading track?)
	if( trk_score_v->at(i_pfp) >= 0.45 ) { 
	  TVector3 trkStart(trk_sce_start_x_v->at(i_pfp), trk_sce_start_y_v->at(i_pfp), trk_sce_start_z_v->at(i_pfp));
	  TVector3 trkEnd(trk_sce_end_x_v->at(i_pfp), trk_sce_end_y_v->at(i_pfp), trk_sce_end_z_v->at(i_pfp));
	  float trk_length = (trkEnd - trkStart).Mag();

	  if( IsInFV_Buffer(trkStart) && IsInFV_Buffer(trkEnd) &&
	      std::min( (trkStart-reco_nu_vtx).Mag(), (trkEnd-reco_nu_vtx).Mag() ) < 4 &&
	      trk_llr_pid_score_v->at(i_pfp) > 0.2 &&
	      trk_length > 10 && 
	      trk_score_v->at(i_pfp) > 0.8)
	  { numuCC = true; }

	// Showers
	} else {
	  TVector3 shrStart(shr_start_x_v->at(i_pfp), shr_start_y_v->at(i_pfp), shr_start_z_v->at(i_pfp));
	  TVector3 shr_p(shr_px_v->at(i_pfp), shr_py_v->at(i_pfp), shr_pz_v->at(i_pfp));   // already normalized
	  //TVector3 shr_distVec = shr_p * shr_dist_v->at(i_pfp);   // ( |Px|*d, |Py|*d, |Pz|*d )
	  //TVector3 shrEnd = shrStart + shr_distVec;
	  float shrE = shr_energy_y_v->at(i_pfp);
	  shr_vtx_P_E.emplace_back(shrStart, shr_p, shrE);   //push_back(std::make_tuple(shrStart, shr_p, shrE));

	  n_shrs++;

	}
      }
      // if the event doesn't pass reco. selection do not include the event
      if( numuCC == false ) continue;

      //============================== PARTICLE LOOP ==============================//
      for(int i=0; i<all_mc_length; i++){

	int mc_trkID = all_mc_trkid->at(i); float mc_E = all_mc_E->at(i);
	int mc_pdg = all_mc_pdg->at(i); int mc_momTrkID = all_mc_mother->at(i);
	bool mc_isPrimary = mc_momTrkID == 0;

	map_allTrkID_partIDX[mc_trkID] = i;
	map_allTrkID_momTrkID[mc_trkID] = mc_momTrkID;

	// Printout this particle information to the screen 
	if( 1 && iEvent < 10 ) {
	  printf("%4i, TrackID: %-4i   PDG: %-12i E=%7.3f MeV,   moth=%5i \n",//, %18s \n",
	    i, mc_trkID, mc_pdg, mc_E, mc_momTrkID//, all_mc_process->at(i).c_str()
	  ); }

        // -------------------------------------------------------------
        // Identify primary neutrons and save their information
        if ( mc_pdg == 2112 ) { if( mc_isPrimary ){ total_nu_Ns++; map_n_p[mc_trkID].clear(); } }

	// Protons from primary neutron mothers
	auto it = map_n_p.find(mc_momTrkID);
        if( !mc_isPrimary && it != map_n_p.end() ) { if( mc_pdg == 2212 ){ it->second.push_back(mc_trkID); } }

      }//end particle loop


      float blipX, blipY, blipZ;
      TVector3 blip_loc, reco_trace_vec;
      float dist, Zdist, blipE, dx, dw, size;
      float density, dist_to_trk, cutoff;
      int pdg, g4id; double XGBoost_score;
      int allIDX, mom_trkID, mom_allIDX, mom_pdg;
      bool isPrimary, mom_isPrimary;
      std::map<int, int> map_blipIDX_score;

      //============================== BLIP LOOP ==============================//
      for( size_t i=0; i<nblips_saved; i++ ){

        blipX = blip_x->at(i), blipY = blip_y->at(i), blipZ = blip_z->at(i);
        TVector3 blip_loc( blipX, blipY, blipZ ); TVector3 reco_trace_vec( blip_loc - reco_nu_vtx );
        dist = reco_trace_vec.Mag(); Zdist = blipZ - reco_nu_vtx.Z();

        blipE = blip_energy->at(i);

        dx = blip_dx->at(i); dw = blip_dw->at(i);
        size = sqrt( dx*dx + dw*dw );

        density = blipE/size; dist_to_trk = blip_proxtrkdist->at(i);
        //cutoff = cutoff_graph->Eval(blipE);

	g4id = blip_true_g4id->at(i); pdg = blip_true_pdg->at(i);

        //============================== CUTS ==============================//
        if( blip_touchtrk->at(i) || dist_to_trk <= 1 || !IsInFV_Exact(blip_loc) ) continue;
        if( blip_pl0_bydeadwire->at(i) == 1 || blip_pl1_bydeadwire->at(i) == 1 || blip_pl2_bydeadwire->at(i) == 1 ) continue;
        if( dist >= 100 ) continue; //if( density < Fit_Function_PID_Cut(blipE) ) continue;
	if( blipE <= 3 ) continue;

	allIDX = -9; mom_trkID = -9; mom_isPrimary = false; mom_allIDX = -9; mom_pdg = -9;
	auto it = map_allTrkID_partIDX.find(g4id);
	if( it != map_allTrkID_partIDX.end() ) {
	  allIDX = it->second;
          mom_trkID = all_mc_mother->at( allIDX );
          isPrimary = mom_trkID == 0;
 
          if( !isPrimary && map_allTrkID_partIDX.count(mom_trkID) != 0 ) {
            mom_allIDX = map_allTrkID_partIDX[mom_trkID];
            mom_isPrimary = all_mc_mother->at( mom_allIDX ) == 0; mom_pdg = all_mc_pdg->at( mom_allIDX );
          }
        }

        //double minDist2Shr = 100000; TVector3 blip_to_shrVtx(0,0,0); double blip_shr_angle = 0;
        //for( const auto& [shr_vtx, shrP, E] : shr_vtx_P_E ) {
        //  if ( IsWithinShowerCone(shr_vtx, shrP, blip_loc, 90.0) && (blip_loc - shr_vtx).Mag() < minDist2Shr ) {
        //    blip_to_shrVtx = blip_loc - shr_vtx;
        //    blip_shr_angle = blip_to_shrVtx.Angle(shrP);
        //    minDist2Shr = blip_to_shrVtx.Mag();
        //  }
        //}

        // RNN model prediction
        double res_angle = -1.0; TVector3 reco_dir; float align_angle, offset;
        if( blip_rnn_dir_isValid->at(i) ) {
          reco_dir.SetXYZ( blip_rnn_dir_x->at(i), blip_rnn_dir_y->at(i), blip_rnn_dir_z->at(i) );
          align_angle = reco_dir.Angle(reco_trace_vec); offset = abs(TMath::Cos(align_angle));   // pointing offset
        } else continue;

	bool isInShower = false;
	for( const auto& [shr_vtx, shr_P, E] : shr_vtx_P_E ) { if ( IsWithinShowerCone(shr_vtx, shr_P, blip_loc, 90.0) ) { isInShower = true; break; } }
	if( isInShower ) continue;

	// BDT prediction; Compute output scores for binary classification and convert 
	blip_E=blipE; dist2Vtx=dist; dist2Trk=dist_to_trk; dist2Edge=distToNearestEdge( blip_loc );
	//blip_X=blipX; blip_Y=blipY; blip_Z=blipZ; blip_dX=dx; blip_dW=dw;
	blip_size=size; blip_density=density; alignment_angle=align_angle; Z_dist2Vtx=Zdist;
	//dist2ShrVtx=blip_to_shrVtx.Mag(); angle_blip_shr=blip_shr_angle;
	double tmva_score = reader->EvaluateMVA( "BDT" ); double base_score = 0.009503577;
	XGBoost_score = 0.5 * TMath::Log( (1 + tmva_score)/(1 - tmva_score) ) - base_score;   // = atanh(tmva_score), sigmoid??
	if( pdg == 2212 && mom_pdg == 2112 && mom_isPrimary ) {   //&& XGBoost_score > -1.97 ) { this was to check P/E values
	  sig_scores.push_back(XGBoost_score);
	  h_signal_BDT_score->Fill(XGBoost_score);
	  if( tagged_neutron_trkIDs.insert(mom_trkID).second ) { neutronID_sig_scores.push_back(XGBoost_score); }
	} else {
	  bkg_scores.push_back(XGBoost_score);
	  h_bkg_BDT_score->Fill(XGBoost_score);
	}

	//if( XGBoost_score > -1.97 ) totalSelected_scores.push_back(XGBoost_score);
	totalSelected_scores.push_back(XGBoost_score);

	//map_blipIDX_score[i] = XGBoost_score;
	//std::transform(bkg_indices.begin(), bkg_indices.end(), 
        //           std::back_inserter(bkg_scores),
        //           [&](const int& key) { return map_blipIDX_score[key]; });

      }//end blip loop

    }//end events loop

    // Efficiency and purity metrics vs. BDT score plots
    int ipoint = 0; float step_size=0.1; float maxExP_val=0.0;
    float maxRareMetric=0, purMaxRare=0, effMaxRare=0, scoreMaxRare=0;
    bool hitOneEff=false; float maxPurNearOneEff=0, effNearOne=0, scoreNearOneEff=0;
    float eff_here=0, pur_here=0, rareSearch=0, ExP_sq=0, effStatErr_here=0;
    float eff_err = 0, pur_err = 0, ExP_err = 0, rareSearch_err = 0, ExP_sq_err = 0;
    // Sort scores for quicker addition
    std::sort(neutronID_sig_scores.begin(), neutronID_sig_scores.end()); std::sort(sig_scores.begin(), sig_scores.end());
    std::sort(totalSelected_scores.begin(), totalSelected_scores.end());
    for( float threshold=-10; threshold<2; threshold+=step_size ) {
      // Add up signal and background at certain thresholds (aka BDT scores)
      auto it = std::lower_bound(neutronID_sig_scores.begin(), neutronID_sig_scores.end(), threshold);
      float num_nuID_sig = std::distance(it, neutronID_sig_scores.end());   // distance to end of vector
      it = std::lower_bound(sig_scores.begin(), sig_scores.end(), threshold);
      float num_sig = std::distance(it, sig_scores.end());
      it = std::lower_bound(totalSelected_scores.begin(), totalSelected_scores.end(), threshold);
      int num_totSelected = std::distance(it, totalSelected_scores.end());
      //std::cout<<"num_nuID_sig: "<<num_nuID_sig<<" num_totSelected: "<<num_totSelected<<"\n";
      //std::cout<<"num_sig: "<<num_sig<<" total_nu_Ns: "<<total_nu_Ns<<"\n";

      if( num_nuID_sig > 0 && num_totSelected > 0 ) {
	// Calculate metrics if valid
	eff_here = num_nuID_sig / total_nu_Ns;
	pur_here = num_sig / num_totSelected;
	rareSearch = num_sig / sqrt( num_sig + (num_totSelected-num_sig) );
	ExP_sq = eff_here*pur_here*pur_here;
	effStatErr_here = sqrt( 1/num_nuID_sig + 1/total_nu_Ns );

	// Binomial errors on eff since total_nu_Ns is fixed (not a source of error)
	eff_err = sqrt( eff_here * (1.f - eff_here) / total_nu_Ns );
        // Poisson errors on purity since the total selected blips does change
	pur_err = pur_here * sqrt( 1/num_sig + 1/num_totSelected );
	// E*P propagated error
	if( eff_here > 0 && pur_here > 0 )
	  ExP_err = eff_here * pur_here * sqrt( pow(eff_err/eff_here, 2) + pow(pur_err/pur_here, 2) );
	// ExP_sq = eff * pur^2  propagated error
	if( eff_here > 0 && pur_here > 0 )
	  ExP_sq_err = ExP_sq * sqrt( pow(eff_err/eff_here, 2) + pow(2*pur_err/pur_here, 2) );
	// rareSearch = s / sqrt(s+b)  Poisson propagated error
	float s = num_sig, b = num_totSelected - num_sig, N = s + b;
	if( N > 0 ) {
	  float dFOM_ds = (0.5f*s + b) / pow(N, 1.5f);
	  float dFOM_db = -0.5f*s    / pow(N, 1.5f);
	  rareSearch_err = sqrt( dFOM_ds*dFOM_ds*s + dFOM_db*dFOM_db*b );
	}

      } else { eff_here = 0.0; pur_here = 0.0; rareSearch = 0.0; ExP_sq = 0.0; } // If invalid, set to 0

      // Fill graphs
      g_pur_vs_BDT_score ->SetPoint     (ipoint, threshold, pur_here);
      g_pur_vs_BDT_score ->SetPointError(ipoint, 0,         pur_err);
      g_eff_vs_BDT_score ->SetPoint     (ipoint, threshold, eff_here * 10);
      g_eff_vs_BDT_score ->SetPointError(ipoint, 0,         eff_err  * 10);
      g_rareSearch_vs_BDT_score->SetPoint(ipoint, threshold, rareSearch / 100);
      g_rareSearch_vs_BDT_score->SetPointError(ipoint, 0,         rareSearch_err / 100);
      g_ExP_sq_vs_BDT_score->SetPoint	(ipoint, threshold, ExP_sq * 100);
      g_ExP_sq_vs_BDT_score->SetPointError(ipoint, 0,         ExP_sq_err * 100);
      float ExP_here = eff_here * pur_here;
      g_ExP_vs_BDT_score ->SetPoint     (ipoint, threshold, ExP_here * 100);
      g_ExP_vs_BDT_score ->SetPointError(ipoint, 0,         ExP_err  * 100);
      //g_effStatErr_vs_BDT_score->SetPoint(ipoint, threshold, effStatErr_here);
      ipoint++;

	std::cout << threshold << "  eff=" << eff_here << "  pur=" << pur_here
           << "  rareSearch=" << rareSearch << "  ExP=" << ExP_here << "\n";

      // Find max E*P across all thresholds
      if( ExP_here > maxExP_val ) maxExP_val = ExP_here;

      // Find the maximum purity achieved before we hit our lowest acceptable efficiency (~0.9%)
      if( eff_here < 0.0085 ) hitOneEff = true; //if( effStatErr_here >= 0.1 ) hitLowEff = true;
      if( !hitOneEff && pur_here > maxPurNearOneEff ) {
	maxPurNearOneEff = pur_here; effNearOne = eff_here;
	scoreNearOneEff = threshold;
      }
      if( rareSearch > maxRareMetric ) {
	purMaxRare = pur_here; effMaxRare = eff_here;
	scoreMaxRare = threshold; maxRareMetric = rareSearch;
      } 
    }

    // Among thresholds within 1% of max E*P, find the one with highest purity
    float bestPur_inExPband = 0; float bestScore_inExPband = 0; float bestExP_inExPband = 0; float bestEff_inExPband = 0;
    for( int ip = 0; ip < g_ExP_vs_BDT_score->GetN(); ip++ ){
        Double_t ExP_here = g_ExP_vs_BDT_score->GetY()[ip];
	Double_t eff_here = g_eff_vs_BDT_score->GetY()[ip] / 10;
        Double_t pur_here = g_pur_vs_BDT_score->GetY()[ip];
        Double_t score_here = g_ExP_vs_BDT_score->GetX()[ip];
        if( ExP_here >= 0.99 * maxExP_val * 100 && pur_here > bestPur_inExPband ){
	    bestEff_inExPband = eff_here;
            bestPur_inExPband = pur_here;
            bestScore_inExPband = score_here;
            bestExP_inExPband = ExP_here;
        }
    }

    printf("\nWithin 1%% of max E*P (%.3fE-03):\nE*P = %.3fE-03\nE: %.3f%%, P: %.3f%%\nscore = %.2f\n\n",
            maxExP_val*1000, bestExP_inExPband*10, bestEff_inExPband*100, bestPur_inExPband*100, bestScore_inExPband);

    printf("Highest purity when efficiency is >= 0.9%%:\nP: %.3f%%\nE*P = %.3fE-03\nscore = %.2f\nE: %.3f%%\n\n",
            maxPurNearOneEff*100, effNearOne*maxPurNearOneEff*1000, scoreNearOneEff, effNearOne*100);

    printf("Max rare search metric (sig/sqrt(sig+bkg)=%.3f):\nE*P = %.3fE-03\nE: %.3f%%, P: %.3f%%\nscore = %.2f\n\n",
	    maxRareMetric, purMaxRare*effMaxRare*1000, effMaxRare*100, purMaxRare*100, scoreMaxRare);

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    // Draw filled bands first (bottom layer)
    auto drawBand = [](TGraphErrors* g, Color_t col, float alpha = 0.25) {
        g->SetFillColorAlpha(col, alpha);
        g->SetLineColorAlpha(col, 0.0);   // invisible line on band pass
        g->Draw("E3 same");
    };
    auto drawLine = [](TGraph* g, Color_t col, int lw = 3) {
        g->SetLineColor(col); g->SetLineWidth(lw);
        g->Draw("L same");
    };

    // Lay down axes with the first graph (no points/line yet)
    g_pur_vs_BDT_score->SetLineColorAlpha(kViolet-5, 0.0);
    g_pur_vs_BDT_score->Draw("AE3");   // A = draw axes
    g_pur_vs_BDT_score->GetXaxis()->SetTitle("BDT Score");
    g_pur_vs_BDT_score->GetYaxis()->SetTitle("Metric Value");
    
    drawBand(g_eff_vs_BDT_score,	kAzure-5);
    drawBand(g_rareSearch_vs_BDT_score,	kBlack,  0.15);
    drawBand(g_ExP_sq_vs_BDT_score,	kSpring-5);
    drawBand(g_ExP_vs_BDT_score,	kOrange-5);
    // pur on top since it's the primary metric
    drawBand(g_pur_vs_BDT_score,	kViolet-5);
    
    // Draw lines on top of all bands
    drawLine(g_pur_vs_BDT_score,	kViolet-5);
    drawLine(g_eff_vs_BDT_score,	kAzure-5);
    drawLine(g_rareSearch_vs_BDT_score,	kBlack);
    drawLine(g_ExP_sq_vs_BDT_score,	kSpring-5);
    drawLine(g_ExP_vs_BDT_score,	kOrange-5);
    
    TLine *cutoffSpot = new TLine(scoreNearOneEff, gPad->GetUymin(), scoreNearOneEff, gPad->GetUymax());
    cutoffSpot->SetLineColor(kRed); cutoffSpot->SetLineWidth(3);
    cutoffSpot->Draw();
    
    // Legend — remove the old effStatErr entry
    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(g_pur_vs_BDT_score,  "Purity (#pm Poisson)");
    leg->AddEntry(g_eff_vs_BDT_score,  "Efficiency #times10 (#pm binomial)");
    leg->AddEntry(g_rareSearch_vs_BDT_score, "S/#sqrt{S+B} /100 (#pm Poisson)");
    leg->AddEntry(g_ExP_sq_vs_BDT_score, "E#timesP^{2} #times100 (#pm propagated)");
    leg->AddEntry(g_ExP_vs_BDT_score,  "E#timesP #times100 (#pm propagated)");
    leg->AddEntry(cutoffSpot,
        TString::Format("P at E ~ 0.9%% = %.2f%%", maxPurNearOneEff*100), "l");
    leg->AddEntry((TObject*)0, TString::Format("BDT cut = %.1f", scoreNearOneEff), "");
    leg->Draw();

    ////Highest purity above 10%% stat. precision in efficiency = %.2f%%", maxPurNearOneEff*100), "l");
    //leg->AddEntry(cutoffSpot, TString::Format("Highest purity when efficiency is at or above 1%% = %.2f%%", maxPurNearOneEff*100), "l");
    //leg->AddEntry((TObject*)0, TString::Format("at BDT cutoff score of %.3f", scoreNearOneEff), "");
    //leg->Draw();

    Double_t *y = g_ExP_vs_BDT_score->GetY();
    Int_t n = g_ExP_vs_BDT_score->GetN();
    Double_t avgExP = TMath::Mean(n, y);
    y = g_pur_vs_BDT_score->GetY(); n = g_pur_vs_BDT_score->GetN();
    Double_t avgP = TMath::Mean(n, y);
    y = g_eff_vs_BDT_score->GetY(); n = g_eff_vs_BDT_score->GetN();
    Double_t avgE = TMath::Mean(n, y);
    std::cout<<"Average E: "<<avgE*10<<"%\nAverage P: "<<avgP*100<<"%\nAverage E*P: "<<avgExP*10<<"E-03\n\n";

    // can try to use std::nth_element for efficiency
    std::sort(bkg_scores.begin(), bkg_scores.end());
    int count = bkg_scores.size(); float BDT_cutoff = bkg_scores[0.995*(count-1)];
    std::cout<<"BDT cutoff: "<<BDT_cutoff<<"\nNumber of bkg blips: "<<count<<"\n\n";

    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    h_bkg_BDT_score->SetLineColor(kOrange); h_bkg_BDT_score->SetFillColorAlpha(kOrange, 0.5);
    h_bkg_BDT_score->Scale(1.0 / h_bkg_BDT_score->Integral("width"));
    h_bkg_BDT_score->DrawCopy("hist");

    h_signal_BDT_score->SetLineColor(kBlue); h_signal_BDT_score->SetFillColorAlpha(kBlue, 0.5);
    h_signal_BDT_score->Scale(1.0 / h_signal_BDT_score->Integral("width"));
    h_signal_BDT_score->DrawCopy("hist same");

    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_signal_BDT_score, "Signal");
    leg1->AddEntry(h_bkg_BDT_score, "Background");
    leg1->AddEntry((TObject*)0, TString::Format("%s %0.3f", "99.5% of background blips are below: ", BDT_cutoff), "");
    leg1->AddEntry((TObject*)0, TString::Format("%s %i", "Number of background blips:", count), "");
    leg1->Draw();
    c1->SaveAs("bkg_sig_BDT_allBlips.jpg");

    //float eff = (float)neutronID_sig_scores.size() / total_nu_Ns;
    //float pur = (float)sig_scores.size() / totalSelected_scores.size();
    // Only count the blips that pass all the cuts including the BDT score cut (implementing BDT cut here)
    BDT_cutoff = std::round(BDT_cutoff*100.0)/100.0;
    float numNeutronsTagged = (float)std::count_if(neutronID_sig_scores.begin(), neutronID_sig_scores.end(), [BDT_cutoff](float s){ return s >= BDT_cutoff; });
    float eff = numNeutronsTagged / total_nu_Ns;
    float numSelectedSigPrBlips = (float)std::count_if(sig_scores.begin(), sig_scores.end(), [BDT_cutoff](float s){ return s >= BDT_cutoff; });
    float totalSelectedBlips = (float)std::count_if(totalSelected_scores.begin(), totalSelected_scores.end(), [BDT_cutoff](float s){ return s >= BDT_cutoff; });
    float pur = numSelectedSigPrBlips / totalSelectedBlips;

    std::cout<<"Number of neutrons tagged: "<<numNeutronsTagged<<"\nTotal neutrino-induced neutrons: "<<total_nu_Ns
	     <<"\nNumber of selected signal proton blips: "<<numSelectedSigPrBlips<<"\nTotal selected blips: "<<totalSelectedBlips<<"\n\n";

    printf( "Efficiency: %5.3f%%   \nPurity: %9.3f%%   \nE*P: %13.3fE-03   \n2*(P*E)/(P+E): %7.5f\n\n",
	     eff*100,
	     pur*100,
	     eff*pur*1000,
	     2*(pur*eff)/(pur+eff)
          );

  }


