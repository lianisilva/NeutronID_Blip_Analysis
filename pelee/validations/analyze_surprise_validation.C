
  #include "TROOT.h"
  #include <cmath>
  #include <algorithm>
  #include "merge_plot.C"

  int maxEvents = -2000;

  int                 selected;
  int                 nblips_saved;
  int		      n_tracks;
  std::vector<float>  * blip_dx;
  std::vector<float>  * blip_dw;
  std::vector<float>  * blip_energy;
  std::vector<float>  * blip_charge;
  std::vector<int>    * blip_nplanes;
  std::vector<float>  * blip_proxtrkdist;
  std::vector<int>    * blip_proxtrkid;
  std::vector<bool>   * blip_touchtrk;
  std::vector<float>  * blip_touchtrkid;
  std::vector<int>    * blip_true_g4id;
  std::vector<int>    * blip_true_pdg;
  std::vector<float>  * blip_true_energy;
  std::vector<float>  * blip_ener;
  std::vector<float>  * blip_x = 0;
  std::vector<float>  * blip_y = 0;
  std::vector<float>  * blip_z = 0;
  std::vector<int>    * blip_pl0_nwires = 0;
  std::vector<int>    * blip_pl1_nwires = 0;
  std::vector<int>    * blip_pl2_nwires = 0;
  std::vector<int>    * blip_pl0_bydeadwire = 0;
  std::vector<int>    * blip_pl1_bydeadwire = 0;
  std::vector<int>    * blip_pl2_bydeadwire = 0;

  std::vector<int>    * all_mc_trkid = 0;
  std::vector<int>    * all_mc_pdg = 0;

  std::vector<int>    * mc_generator_trackid = 0;
  std::vector<int>    * mc_generator_pdg = 0;
  std::vector<float>  * mc_generator_E = 0;

  std::vector<int>    * mc_trackid = 0;
  std::vector<int>    * mc_pdg = 0;
  std::vector<float>  * mc_vx = 0;
  std::vector<float>  * mc_vy = 0;
  std::vector<float>  * mc_vz = 0;
  std::vector<float>  * mc_endx = 0;
  std::vector<float>  * mc_endy = 0;
  std::vector<float>  * mc_endz = 0;
  std::vector<float>  * mc_E = 0;

  //Int_t    reco_id[1000];   // these are ints not vectors!

  //bool   nu_isNeutrino;
  float  reco_nu_vtx_x;
  float  reco_nu_vtx_y;
  float  reco_nu_vtx_z;

  TH1D* h_nblips_vs_vtxDist;
  TH1D* h_blipdx;
  TH1D* h_blipdw;
  TH1D* h_blipE;
  TH1D* h_nBlips;
  TH1D* h_BlipsPerEvent;
  TH1D* h_charge;
  TH1D* h_nPlanes;
  TH1D* h_proxTrkID;
  TH1D* h_proxTrkDist;
  TH1D* h_touchTrkID;
  TH1D* h_maxWires;
  TH1D* h_nwires_pl0;
  TH1D* h_nwires_pl1;
  TH1D* h_nwires_pl2;
  TH1D* h_trueG4ID;
  TH1D* h_truePDG;
  TH1D* h_trueE;
  TH3D* h_blipXYZ;
  TH2D* h_blipZY;
  TH1D* h_blipX;
  TH1D* h_blipY;
  TH1D* h_blipZ;
  TH1D* h_trueMatchX;
  TH1D* h_trueMatchY;
  TH1D* h_trueMatchZ;
  TH2D* h_trueMatchZY;
  TH2D* h_withMatchZY;
  TH2D* h_nonTrueMatchZY;
  TH2D* h_trueMatchDistZY;
  TH2D* h_withMatchDistZY;
  TH2D* h_nonTrueMatchDistZY;
  TH1D* h_blipXdistToVtx;
  TH1D* h_blipYdistToVtx;
  TH1D* h_blipZdistToVtx;
  TH2D* h_blipZYdistToVtx;

  // Input file
  std::string fFileName = "checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root";   // retupled BNB nu
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";   // BNB nu
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_11_BNB_beam_on_surprise_reco2_hist.root";   // BNB beam-on data
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root";   // BNB beam-off
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_09_BNB_nue_overlay_surprise_reco2_hist.root";   // nue
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root";   // dirt
  //std::string fFileName = "MCC9.10_Run4b_v10_04_07_09_BNB_NC_pi0_overlay_surprise_reco2_hist.root";   // NC pi0
  //std::string fFileName = "MCC9.10_Run4b_NuMI_beam_on_RHC_data_surprise_v10_04_07_09_reco2_hist.root"; // NUMI beam-on
  //std::string fFileName = "MCC9.10_Run4b_NuMI_beam_off_RHC_data_surprise_v10_04_07_09_reco2_hist.root"; // NUMI beam-off
  //std::string fFileName = "MCC9.10_Run4b_NuMI_RHC_nu_overlay_surprise_v10_04_07_09_reco2_hist.root"; // NUMI nu
  //std::string fFileName = "MCC9.10_Run4b_NuMI_RHC_nue_overlay_surprise_v10_04_07_09_reco2_hist.root"; // NUMI nue

  std::string nuTreeName = "nuselection/NeutrinoSelectionFilter";
  //std::string wcTreeName = "wcpselection/T_PFeval";

  // ********************************************************
  void analyze_surprise_validation(){

    // Read in the TFile1000 and grab the TTree
    TFile* file = new TFile(("files/"+fFileName).c_str(),"READ");
    TTree* nuTree = (TTree*)file->Get(nuTreeName.c_str());
    //TTree* wcTree = (TTree*)file->Get(wcTreeName.c_str());
    std::cout<<"Reading data file "<<fFileName<<"\n";

    // Set branch address for all variables you plan to use
    nuTree->SetBranchAddress("selected",	    &selected);
    nuTree->SetBranchAddress("n_tracks",	    &n_tracks);
    nuTree->SetBranchAddress("nblips_saved",	    &nblips_saved);
    nuTree->SetBranchAddress("blip_true_g4id",	    &blip_true_g4id);
    nuTree->SetBranchAddress("blip_x",		    &blip_x);
    nuTree->SetBranchAddress("blip_y",		    &blip_y);
    nuTree->SetBranchAddress("blip_z",		    &blip_z);
    nuTree->SetBranchAddress("blip_dx",             &blip_dx);
    nuTree->SetBranchAddress("blip_dw",             &blip_dw);
    nuTree->SetBranchAddress("blip_energy",	    &blip_energy);
    nuTree->SetBranchAddress("blip_charge",         &blip_charge);
    nuTree->SetBranchAddress("blip_nplanes",        &blip_nplanes);
    nuTree->SetBranchAddress("blip_proxtrkdist",    &blip_proxtrkdist);
    nuTree->SetBranchAddress("blip_proxtrkid",      &blip_proxtrkid);
    nuTree->SetBranchAddress("blip_touchtrkid",     &blip_touchtrkid);
    nuTree->SetBranchAddress("blip_touchtrk",       &blip_touchtrk);
    nuTree->SetBranchAddress("blip_true_g4id",      &blip_true_g4id);
    nuTree->SetBranchAddress("blip_true_pdg",       &blip_true_pdg);
    nuTree->SetBranchAddress("blip_true_energy",    &blip_true_energy);

    nuTree->SetBranchAddress("blip_pl0_nwires",	    &blip_pl0_nwires);
    nuTree->SetBranchAddress("blip_pl1_nwires",	    &blip_pl1_nwires);
    nuTree->SetBranchAddress("blip_pl2_nwires",	    &blip_pl2_nwires);
    nuTree->SetBranchAddress("blip_pl0_bydeadwire", &blip_pl0_bydeadwire);
    nuTree->SetBranchAddress("blip_pl1_bydeadwire", &blip_pl1_bydeadwire);
    nuTree->SetBranchAddress("blip_pl2_bydeadwire", &blip_pl2_bydeadwire);

    nuTree->SetBranchAddress("all_mc_trkid",	    &all_mc_trkid);
    nuTree->SetBranchAddress("all_mc_pdg",	    &all_mc_pdg);

    nuTree->SetBranchAddress("mc_generator_trackid",&mc_generator_trackid);
    nuTree->SetBranchAddress("mc_generator_pdg",    &mc_generator_pdg);
    nuTree->SetBranchAddress("mc_generator_E",	    &mc_generator_E);

    nuTree->SetBranchAddress("mc_trackid",	    &mc_trackid);
    nuTree->SetBranchAddress("mc_pdg",		    &mc_pdg);
    nuTree->SetBranchAddress("mc_vx",               &mc_vx);
    nuTree->SetBranchAddress("mc_vy",               &mc_vy);
    nuTree->SetBranchAddress("mc_vz",               &mc_vz);
    nuTree->SetBranchAddress("mc_endx",             &mc_endx);
    nuTree->SetBranchAddress("mc_endy",             &mc_endy);
    nuTree->SetBranchAddress("mc_endz",             &mc_endz);
    nuTree->SetBranchAddress("mc_E",		    &mc_E);

    //tree->SetBranchAddress("nu_isNeutrino",       &nu_isNeutrino);
    nuTree->SetBranchAddress("reco_nu_vtx_x",       &reco_nu_vtx_x);
    nuTree->SetBranchAddress("reco_nu_vtx_y",       &reco_nu_vtx_y);
    nuTree->SetBranchAddress("reco_nu_vtx_z",       &reco_nu_vtx_z);
    nuTree->Draw("blip_touchtrk >> touchtrk", "", "goff");   // ROOT created hist is auto named htemp
    TH1F *h_touchtrk = (TH1F*)gDirectory->Get("touchtrk");
    nuTree->Draw("blip_pl0_bydeadwire >> nearDeath0", "", "goff"); TH1F *h_nearDeath0 = (TH1F*)gDirectory->Get("nearDeath0");
    nuTree->Draw("blip_pl1_bydeadwire >> nearDeath1", "", "goff"); TH1F *h_nearDeath1 = (TH1F*)gDirectory->Get("nearDeath1");
    nuTree->Draw("blip_pl2_bydeadwire >> nearDeath2", "", "goff"); TH1F *h_nearDeath2 = (TH1F*)gDirectory->Get("nearDeath2");

    //wcTree->SetBranchAddress("reco_id",		    &reco_id);

    // make output file to before declaring plots
    TFile* outfile = new TFile("test_out.root", "recreate");
  
    // Make histograms
    h_BlipsPerEvent     = new TH1D("BlipsPerEvent","Number of Blips per Event;# Blips",250,0,500);
    //h_blipXYZ		= new TH3D("blipXYZ","Blip XYZ;x;y;z",80,0,256,80,-116,116,80,0,1036);
    h_blipX             = new TH1D("blipX","Blip X;x",45,-100,350);   // 256
    h_blipY             = new TH1D("blipY","Blip Y;y",28,-140,140);   // -/+ 116
    h_blipZ             = new TH1D("blipZ","Blip Z;z",115,-50,1100);   // 1036
    h_blipZY            = new TH2D("blipZY","Blip Z vs Y;z;y",224,-40,1080,50,-150,150);  // 0/1036, -116/116
    h_trueMatchX        = new TH1D("trueMatchX","Truth Matched Blip X;x",45,-100,350);
    h_trueMatchY        = new TH1D("trueMatchY","Truth Matched Blip Y;y",28,-140,140);
    h_trueMatchZ        = new TH1D("trueMatchZ","Truth Matched Blip Z;z",115,-50,1100);
    h_trueMatchZY       = new TH2D("trueMatchZY","Truth Matched Blip Z vs Y;z;y",224,-40,1080,50,-150,150);
    h_withMatchZY       = new TH2D("withMatchZY","Blip Matched to MC Particle Z vs Y;z;y",224,-40,1080,50,-150,150);
    h_nonTrueMatchZY    = new TH2D("nonTrueMatchZY","Non-Truth Matched Blip Z vs Y;z;y",224,-40,1080,50,-150,150);   //blips without trueG4ID
    //h_blipXY            = new TH2D("blipXY","Blip XY;x;y",80,0,256,80,-116,116);
    //h_blipXZ            = new TH2D("blipXZ","Blip XZ;x;z",80,0,256,80,0,1036);
    h_blipdx		= new TH1D("blipdx","Blip dX;dX",25,0,5);
    h_blipdw            = new TH1D("blipdw","Blip dW;dW",25,0,5);
    h_blipE             = new TH1D("blipE","Reconstructed Blip Energy;E [MeVee]",200,-10,30);
    h_charge            = new TH1D("charge","Blip Charge;Charge [e-]",500,0,500000);
    h_nPlanes           = new TH1D("nPlanes","Number of Planes a Blip Spans;# Planes",4,0,4);
    h_proxTrkDist       = new TH1D("proxTrkDist","Proximity to Nearest Track;Distance [cm]",520,-20,500);
    h_proxTrkID         = new TH1D("proxTrkID","Nearest Track ID;Pandora track ID",65,-5,60);
    h_touchTrkID        = new TH1D("touchTrkID","Touching Track ID;Pandora track ID",35,-10,25);
    h_nwires_pl0	= new TH1D("nwires_pl0","Number of Wires per Blip on Plane 0;# Wires;# Blips",18,0,18);
    h_nwires_pl1	= new TH1D("nwires_pl1","Number of Wires per Blip on Plane 1;# Wires;# Blips",18,0,18);
    h_nwires_pl2	= new TH1D("nwires_pl2","Number of Wires per Blip on Plane 2;# Wires;# Blips",18,0,18);
    h_maxWires        	= new TH1D("maxWires","Max Wires per Blip on one Plane;# Wires",18,0,18);
    h_trueG4ID          = new TH1D("trueG4ID","True Blip G4ID;G4ID",4200,-200,4000);
    h_truePDG           = new TH1D("truePDG","True Blip PDG;PDG",18800000,-9400000,9400000);
    h_trueE             = new TH1D("trueE","True Blip Energy;E [MeV]",200,-10,30); 

    h_blipXdistToVtx    = new TH1D("blipXdistToVtx","Blip X Distance to Reco Vertex;x dist",125,-250,250);  // seem to go to 220, X & Y, Z goes all the way ;)
    h_blipYdistToVtx    = new TH1D("blipYdistToVtx","Blip Y Distance to Reco Vertex;y dist",125,-250,250);
    h_blipZdistToVtx    = new TH1D("blipZdistToVtx","Blip Z Distance to Reco Vertex;z dist",125,-250,250);
    h_blipZYdistToVtx   = new TH2D("blipZYdistToVtx","Blip Z vs. Y Distance to Reco Vertex;Z Distance [cm];Y Distance [cm];",224,-40,1080,50,-150,150);
    h_trueMatchDistZY   = new TH2D("trueMatchDistZY","Truth Matched Blip Z vs. Y Distance to Reco Vertex;Z Distance [cm];Y Distance [cm];",224,-40,1080,50,-150,150);
    h_withMatchDistZY   = new TH2D("withMatchDistZY","Blip Matched to MC Particle Z vs. Y Distance to Reco Vertex;Z Distance [cm];Y Distance [cm];",224,-40,1080,50,-150,150);
    h_nonTrueMatchDistZY    = new TH2D("nonTrueMatchDistZY","Non-Truth Matched Blip Z vs. Y Distance to Reco Vertex;Z Distance [cm];Y Distance [cm];",224,-40,1080,50,-150,150);

    h_nblips_vs_vtxDist = new TH1D("nblips_vs_vtxDist",";Distance to reco vtx [cm]; Nblips / r^2",250,0,250);

    //h_blipXYZ->SetOption("colz");
    //h_blipXY->SetOption("colz");
    h_blipZY->SetOption("colz");
    //h_blipXZ->SetOption("colz");

    float minX, minY, minZ = 0;
    float maxX, maxY, maxZ = 0;
 
    int numEvents = nuTree->GetEntries();
    std::cout<<"There are "<<numEvents<<" events in the nu tree.\n";
    //std::cout<<"There are "<<wcTree->GetEntries()<<" events in the WC tree.\n";  

    int totalSampleBlips=0;
    float blipsWg4id=0;
    float non_p_e_blipsWg4id=0;
    float blipsMatchAll=0;
    float blipsMConly=0;
    float blipsGenOnly=0;
    float blipsAllOnly=0;
    //float blipsWCmatch=0;
    float blipsMatch=0;
    float blipsNoMatch=0;
    float totalBlipE=0;
    int nonSelected=0;
    int selectedEvents=0;

    float dist_sum=0;
    //std::vector<float> numTrks;
    std::map<int, int> map_allTrkID_partIDX;
    std::map<int, int> map_MCtrkID_partIDX;
    std::map<int, int> map_genTrkID_partIDX; 
    //std::map<int, int> map_WCtruthID_partIDX;

    //====================================================
    // Loop through the events. 
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      //if( iEvent > 10000 ) continue;
      nuTree->GetEntry(iEvent);
      //wcTree->GetEntry(iEvent);

      map_allTrkID_partIDX.clear();
      map_MCtrkID_partIDX.clear();
      map_genTrkID_partIDX.clear();
      //map_WCtruthID_partIDX.clear();

      int mc_length = static_cast<int>(mc_trackid->size());
      for(int i=0; i<mc_length; i++){ map_MCtrkID_partIDX[mc_trackid->at(i)] = i; }

      int all_mc_length = static_cast<int>(all_mc_trkid->size());
      for(int i=0; i<all_mc_length; i++){ map_allTrkID_partIDX[all_mc_trkid->at(i)] = i; }

      int mc_gen_length = static_cast<int>(mc_generator_trackid->size());
      for(int i=0; i<mc_gen_length; i++){ map_genTrkID_partIDX[mc_generator_trackid->at(i)] = i; }

      //int WC_length = 1000;
      //for(int i=0; i<WC_length; i++){ map_WCtruthID_partIDX[reco_id[i]] = i; }

      //TVector3 reco_nu_vtx( reco_nu_vtx_x, reco_nu_vtx_y, reco_nu_vtx_z );
      //bool hasRecoVtx = false;
      //if( fabs(reco_nu_vtx_x)<1000 && fabs(reco_nu_vtx_y)<1000 && reco_nu_vtx_z>-100 ) hasRecoVtx = true;
      //if( hasRecoVtx == false || nblips_saved<=0 ) { nonSelected++; continue; }
      selectedEvents++;

      if(iEvent % 1000 == 0) {
	std::cout<<"\n=======================================\n";
	std::cout<<"Reading event "<<iEvent<<" out of "<<numEvents<<"\n";
	std::cout<<"NBlips: "<<nblips_saved<<"   NParticles (mc_var): "<<mc_length<<"\n"; 
      }

      for(size_t i=0; i<nblips_saved; i++){

	float blipX = blip_x->at(i), blipY = blip_y->at(i), blipZ = blip_z->at(i);

	if( minX > blipX ) minX=blipX; if( minY > blipY ) minY=blipY; if( minZ > blipZ ) minZ=blipZ;
	if( maxX < blipX ) maxX=blipX; if( maxY < blipY ) maxY=blipY; if( maxZ < blipZ ) maxZ=blipZ;
        //TVector3 blip_loc( blipX, blipY, blipZ );
	//float dist=(blip_loc-reco_nu_vtx).Mag();
	float dx = blipX - reco_nu_vtx_x; float dy = blipY - reco_nu_vtx_y; float dz = blipZ - reco_nu_vtx_z;
	float dist = sqrt(dx*dx + dy*dy + dz*dz); dist_sum+=(float)dist;
	bool hasMCmatch = false; bool hasGenMatch = false; bool hasAllMatch = false; //bool hasWCmatch = false;
	float blipE = blip_energy->at(i);
	if(i % 10 == 0 || 1) {
	  //h_nblips_vs_vtxDist->Fill( dist, 1./pow(dist,2) );   //((4./3.)*M_PI*pow(dist,3)) );   // normalize by 1/r^2, make spherical dist. flat
	  //h_blipE->Fill(blipE);
	  totalBlipE+=blipE;

	  //h_nwires_pl0->Fill( blip_pl0_nwires->at(i) < 1 ? 0 : blip_pl0_nwires->at(i));
	  //h_nwires_pl1->Fill( blip_pl1_nwires->at(i) < 1 ? 0 : blip_pl1_nwires->at(i)); 
	  //h_nwires_pl2->Fill( blip_pl2_nwires->at(i) < 1 ? 0 : blip_pl2_nwires->at(i)); 
	  //int maxWires = std::max({blip_pl0_nwires->at(i), blip_pl1_nwires->at(i), blip_pl2_nwires->at(i)});

	  //h_maxWires->Fill(maxWires); //h_blipXYZ->Fill(blipX,blipY,blipZ);
	  //h_blipZY->Fill(blipZ,blipY); h_blipX->Fill(blipX); h_blipY->Fill(blipY); h_blipZ->Fill(blipZ);
	  //h_blipXdistToVtx->Fill(blipX - reco_nu_vtx_x); h_blipYdistToVtx->Fill(blipY - reco_nu_vtx_y); h_blipZdistToVtx->Fill(blipZ - reco_nu_vtx_z); h_blipZYdistToVtx->Fill(blipZ - reco_nu_vtx_z,blipY - reco_nu_vtx_y);

	  //h_blipdx->Fill(blip_dx->at(i)); h_blipdw->Fill(blip_dw->at(i));
          //h_charge->Fill(blip_charge->at(i)); h_nPlanes->Fill(blip_nplanes->at(i));
          //h_proxTrkDist->Fill(blip_proxtrkdist->at(i)); h_proxTrkID->Fill(blip_proxtrkid->at(i)); h_touchTrkID->Fill(blip_touchtrkid->at(i));
          //h_trueG4ID->Fill(blip_true_g4id->at(i)); h_truePDG->Fill(blip_true_pdg->at(i)); h_trueE->Fill(blip_true_energy->at(i));
	}

	int blip_g4id = blip_true_g4id->at(i);
	if( blip_g4id>=0 ) {   // choose whether to fill histogram
	  //if( blip_true_pdg->at(i) != all_mc_pdg->at(i) ) std::cout<<"WOAH, HOLD UR HORSES! blip_true_pdg &&  all_mc_pdg don't match!!"<<endl;
	  blipsWg4id++;   // add up truth matched blips 
	  if( abs(blip_true_pdg->at(i)) != 11 && blip_true_pdg->at(i) != 22 ) non_p_e_blipsWg4id++;
	  //h_trueMatchX->Fill(blipX);
          //h_trueMatchY->Fill(blipY);
          //h_trueMatchZ->Fill(blipZ);
          //h_trueMatchZY->Fill(blipZ, blipY);
          //h_trueMatchDistZY->Fill(blipZ - reco_nu_vtx_z,blipY - reco_nu_vtx_y);
  	  if( map_MCtrkID_partIDX.find( blip_g4id ) != map_MCtrkID_partIDX.end() ) { hasMCmatch=true;
	    //h_withMatchZY->Fill(blipZ, blipY); h_withMatchDistZY->Fill(blipZ - reco_nu_vtx_z, blipY - reco_nu_vtx_y);
	  }
	  if( map_genTrkID_partIDX.find( blip_g4id ) != map_genTrkID_partIDX.end() ) { hasGenMatch=true; }
	  if( map_allTrkID_partIDX.find( blip_g4id ) != map_allTrkID_partIDX.end() ) { hasAllMatch=true; }
          //if( map_WCtruthID_partIDX.find( blip_g4id ) != map_WCtruthID_partIDX.end() ) { hasWCmatch=true; blipsWCmatch++; }

          if( hasMCmatch || hasGenMatch || hasAllMatch ) {
	    blipsMatch++;
	    if( hasMCmatch && hasGenMatch && hasAllMatch ) { blipsMatchAll++; }
	    else if( hasMCmatch && (!hasAllMatch && !hasGenMatch) ) { blipsMConly++; }
	    else if( hasGenMatch && (!hasAllMatch && !hasMCmatch) ) { blipsGenOnly++; }
	    else if( hasAllMatch && (!hasGenMatch && !hasMCmatch) ) { blipsAllOnly++; }
	  } else { blipsNoMatch++; }


	} //else { h_nonTrueMatchZY->Fill(blipZ, blipY); h_nonTrueMatchDistZY->Fill(blipZ - reco_nu_vtx_z, blipY - reco_nu_vtx_y); }

      }//end blip loop

      h_BlipsPerEvent->Fill(nblips_saved);
      totalSampleBlips+=nblips_saved;
 
    }//end events loop

    // normalize for nBlips per bin per event
    h_BlipsPerEvent->Scale(1.0 / numEvents);
    h_blipX->Scale(1.0 / numEvents); h_blipY->Scale(1.0 / numEvents); h_blipZ->Scale(1.0 / numEvents);
    h_blipdx->Scale(1.0 / numEvents); h_blipdw->Scale(1.0 / numEvents);
    h_blipE->Scale(1.0 / numEvents); h_charge->Scale(1.0 / numEvents);
    h_nPlanes->Scale(1.0 / numEvents);
    h_proxTrkDist->Scale(1.0 / numEvents); h_proxTrkID->Scale(1.0 / numEvents); h_touchTrkID->Scale(1.0 / numEvents);
    h_nwires_pl0->Scale(1.0 / numEvents); h_nwires_pl1->Scale(1.0 / numEvents); h_nwires_pl2->Scale(1.0 / numEvents);
    h_maxWires->Scale(1.0 / numEvents);
    h_trueG4ID->Scale(1.0 / numEvents); h_truePDG->Scale(1.0 / numEvents); h_trueE->Scale(1.0 / numEvents);

    h_blipXdistToVtx->Scale(1.0 / numEvents); h_blipYdistToVtx->Scale(1.0 / numEvents); h_blipZdistToVtx->Scale(1.0 / numEvents);
    h_blipZYdistToVtx->Scale(1.0 / numEvents); h_nblips_vs_vtxDist->Scale(1.0 / numEvents);

    if( nblips_saved == 0) std::cout<<"No blips saved";

    printf("\nBlip multiplicity: %8.3f +/- %.3f\n",
	   double(totalSampleBlips)/double(selectedEvents), h_BlipsPerEvent->GetRMS(11) );

    std::cout<<"\nTotal Events: "<<numEvents<<" events\n"
           <<"Selected Events: "<<selectedEvents<<" events\n"
	   <<"Selected events without blips: "<<nonSelected<<" events\n"
           <<"Percent of selected events without reco nu vtx: "<<nonSelected*100/(selectedEvents+nonSelected)<<"%\n\n"

	   <<"Total blips: "<<totalSampleBlips<<" blips\n"
           //<<"Number of blips with a true g4id: "<<blipsWg4id<<" blips\n"
	   //<<"Number of blips matched to EITHER (mc_ or gen_) trackids: "<<blipsMatch<<" blips\n"
           //<<"Number of blips ONLY matched to an mc_trackid: "<<blipsMConly<<" blips\n"
           //<<"Number of blips ONLY matched to an mc_GEN_trackid: "<<blipsGenOnly<<" blips\n"
           //<<"Number of blips matched to WC reco_id: "<<blipsWCmatch<<" blips\n"
           //<<"Number of blips matched to BOTH (mc_ and gen_) trackids: "<<blipsMatchBoth<<" blips\n"
           //<<"Number of blips NOT matched to either (mc_ and gen_) trackids: "<<blipsNoMatch<<" blips\n"
	   <<"Percent of total blips with a true g4id: "<<(blipsWg4id*100)/totalSampleBlips<<"%\n"
	   <<"Percent of true g4id blips matched to ANY trackid: "<<(blipsMatch*100)/blipsWg4id<<"%\n"
           <<"Percent of true g4id blips ONLY matched to an MC_trackid: "<<(blipsMConly*100)/blipsWg4id<<"%\n"
           <<"Percent of true g4id blips ONLY matched to an mc_GEN_trackid: "<<(blipsGenOnly*100)/blipsWg4id<<"%\n"
           <<"Percent of true g4id blips ONLY matched to an ALL_mc_trackid: "<<(blipsAllOnly*100)/blipsWg4id<<"%\n"
           <<"Percent of true g4id blips matched to EVERY trackid: "<<(blipsMatchAll*100)/blipsWg4id<<"%\n"
           <<"Percent of true g4id blips NOT matched to ANY trackid: "<<(blipsNoMatch*100)/blipsWg4id<<"%\n"
           <<"Percent of total blips with a true g4id ALSO not protons/electrons: "<<(non_p_e_blipsWg4id*100)/totalSampleBlips<<"%\n"
	   <<"Percent of true g4id blips that are NON-p/e: "<<(non_p_e_blipsWg4id*100)/blipsWg4id<<"%\n"
           <<"Percent of non-p/e true g4id blips matched to ANY trackid: "<<(blipsMatch*100)/non_p_e_blipsWg4id<<"%\n\n"

	   <<"Maximum blip (x, y, z): ("<<maxX<<", "<<maxY<<", "<<maxZ<<")\n"
           <<"Minimum blip (x, y, z): ("<<minX<<", "<<minY<<", "<<minZ<<")\n"
	   <<"Mean blip (x, y, z): ("<<h_blipX->GetMean()<<", "<<h_blipY->GetMean()<<", "<<h_blipZ->GetMean()<<")\n\n"

           <<"Mean dist to vtx: "<<h_nblips_vs_vtxDist->GetMean()<<"\n"
           <<"Mean reconstructed energy: "<<h_blipE->GetMean()<<"\n"
           <<"Mean #blips by dead wires plane 0: "<<h_nearDeath0->GetMean()<<"\n"
           <<"Mean #blips by dead wires plane 1: "<<h_nearDeath1->GetMean()<<"\n"
           <<"Mean #blips by dead wires plane 2: "<<h_nearDeath2->GetMean()<<"\n"
	   <<"Mean #blips touching tracks: "<<h_touchtrk->GetMean()<<"\n\n";

    //make this^ have more decimals, something with cout...
  
    //TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    //h_numTrks->DrawCopy();
  
    merge_plot(h_nwires_pl0, "plane 0", h_nwires_pl1, "plane 1", h_nwires_pl2, "plane 2");
  
    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    h_BlipsPerEvent->DrawCopy("hist");
  
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
