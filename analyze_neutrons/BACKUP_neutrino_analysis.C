
  //#include "variables.h"
  //#include "../rnn_inference/rnn_inference"
  #include "../rnn_inference/blip_utils.h"
  #include "accessory_macros/set_branches.cpp"
  //#include "accessory_macros/ProcessHandler.h"
  #include "accessory_macros/merge_plot.C"
  #include "accessory_macros/Visualizer3D.C"
  #include <iostream>
  #include <cmath>
  #include <vector>
  #include <map>
  #include <unordered_set>
  #include <iterator>
  #include <algorithm>
  #include "TMath.h"
  #include "TROOT.h"
  #include "TVectorT.h"
  #include "TDirectory.h"
  #include "TArrow.h"
 
  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;
  
  // Input file name
  //std::string fFileName = "blipana_protons_0to200MeV_overlay_20250826.root";
  //std::string fFileName = "blipana_MCC9_bnb_nu_overlay_mcc9_run3_20250904.root";   // BNB nu
  std::string fFileName = "blipana_prodgenie_bnb_nu_overlay_mcc9_run3_numuCCinclusive_selection_20251028.root"; //incl.
  //std::string fFileName = "remains/all_blips.root"; //inclusive remains
  bool isOutTree = false; const char* model_path = "../rnn_inference/model_traced.pt";

  TFile* outfile;
  TFile* file;
  TTree* tree; 
  int numEvents;
  TDirectory* blipE_plots;
  TDirectory* recoVtxDist_plots;
  TDirectory* cosAlignment_plots;
  TDirectory* alignmentAngle_plots;
  TDirectory* resolution_plots;
  TDirectory* res_vs_E_plots;
  TDirectory* Edensity_plots;
  TDirectory* proxTrkDist_plots;
  TDirectory* byDeadWire_plots;
  TDirectory* cosmicDistToVtx_plots;
  TDirectory* dx_plots;
  TDirectory* dw_plots;
  TDirectory* size_plots;
  TDirectory* Nhits_plots;

  const int numHists=10;
  TH1D*   h_all_nu_neutrons;	    // neutrino-induced neutrons (NINs)
  TH1D*   h_trueNeutronE;
  TH1D*   h_preFoundNeutronE;
  TH1D*   h_foundNeutronE;
  TH1D*   h_NC_Neutrons;	    // Neutral current NINs per neutrino
  TH1D*   h_CC_Neutrons;	    // Charged current NINs per neutrino
  TH1D*   h_maxWires;
  TH2D*   h_blipZY;
  TH1D*   h_numTrks;
  TH1D*   h_trackLength;
  TH1D*   h_trackEnergy;
  TH1D*   h_prTrksPerE;
  TH1D*   h_prBlipsPerE;
  //TH1D*   h_trueProtonE;
  //TH1D*   h_truePrimPrE;
  TH1D*   h_trueBlipE;
  TH1D*   h_neutronTrueE;
  TH1D*   h_gammaTrueE;
  TH1D*   h_signalBlipE;
  TH2D*   h_trueNeutronEvsProtonE;  // neutron kinetic energy vs resultant proton KE
  TH1D*   h_trueRecoVtxAgreement;
  //TH1D*   h_edeps_vs_vtxDist;
  //TH1D*   h_edep_energy_dist;
  TH1D*   h_blipRecoEff_vs_vtxDist;
  TH2D*   h_blipRecoEff_energy_dist;
  TH1D*   h_allPrBlips_vs_vtxDist;
  TH1D*   h_truePr_vs_vtxDist;
  TH1D*   h_truePrCosAlign;
  TH1D*   h_cosAlignment[numHists];
  TH1D*   h_alignmentAngle[numHists];
  TH2D*   h_cos_vs_energy;
  TH2D*   h_pr_cos_vs_energy;
  TH2D*   h_pr_E_ds;
  TGraph* gr_signal_E_ds;
  TH2D*   h_elec_E_ds;
  TGraph* gr_elec_E_ds;
  TGraph* gr_cosmic_E_ds;
  TH2D*   h_Np_vs_Nn;
  TH2D*   h_res_vs_trueE;
  TH2D*   h_pr_res_vs_trueE;
  TGraph* gr_pr_res_vs_trueE;
  TH2D*   h_res_vs_E[numHists];
  TH1D*   h_resolution[numHists];
  TH1D*   h_allRes;
  TH1D*   h_cosmicRes;
  TH1D*   h_neutronRes;
  TH1D*   h_gammaRes;
  TH1D*   h_prRes;
  TH1D*   h_prNonPrimRes;
  TH1D*   h_primPrRes;
  TH1D*   h_nonNeutRes;
  TH1D*   h_blipE[numHists];
  TH1D*   h_recoVtxDist[numHists];
  TH1D*   h_Edensity[numHists];
  TH1D*   h_proxTrkDist[numHists];
  TH1D*   h_byDeadWire[numHists];
  TH1D*   h_cosmicDistToVtx[numHists];
  TH1D*   h_dx[numHists];
  TH1D*   h_dw[numHists];
  TH1D*   h_size[numHists];
  TH1D*   h_Nhits[numHists];

  std::map<int, int> map_g4id_mother;
  std::map<int, int> map_g4id_pdg; 
  std::map<int, int> map_g4id_index;
  // keep track of neutrons and map primary neutrons to their daughters via INDEX
  std::map<int, std::vector<int>> map_n_p;
  std::map<int, std::vector<int>> map_n_n;
  std::map<int, int> map_nu_N_g4id_daught_idx; 
  std::vector<float> E_bins;

  //reset_vars();
  void initialize() {
    // Read in the TFile and grab the TTree
    file = new TFile(fFileName.c_str(),"READ");
    tree = (TTree*)file->Get("blipana/anatree");
    //tree = (TTree*)file->Get("what_remains");
    numEvents = tree->GetEntries();
    std::cout<<"Reading data file "<<fFileName<<"\n";
    // make output file to store plots
    outfile = new TFile("xSecs/temp_out.root", "recreate");
    blipE_plots		  = outfile->mkdir("blipE_plots");
    recoVtxDist_plots	  = outfile->mkdir("recoVtxDist_plots");
    cosAlignment_plots	  = outfile->mkdir("cosAlignment_plots");
    alignmentAngle_plots  = outfile->mkdir("alignmentAngle_plots");
    resolution_plots	  = outfile->mkdir("resolution_plots");
    res_vs_E_plots	  = outfile->mkdir("res_vs_E_plots");
    Edensity_plots	  = outfile->mkdir("Edensity_plots");
    proxTrkDist_plots	  = outfile->mkdir("proxTrkDist_plots");
    byDeadWire_plots	  = outfile->mkdir("byDeadWire_plots");
    cosmicDistToVtx_plots = outfile->mkdir("cosmicDistToVtx_plots");
    dx_plots		  = outfile->mkdir("dx_plots");
    dw_plots		  = outfile->mkdir("dw_plots");
    size_plots		  = outfile->mkdir("size_plots");
    Nhits_plots           = outfile->mkdir("Nhits_plots");
    get_branches(tree, isOutTree);   // does this need to go after outfile?
  }

  //====================================================================
  // Function to determine if a particle descended from another particle.
  // Allows option to break lineage at photons for contiguous parentage.
  bool IsAncestorOf(int particleID, int ancestorID, bool breakAtPhots = false){
    if( particleID == ancestorID  )       return true;
    if( particleID < ancestorID   )       return false;
    if( map_g4id_mother.count(ancestorID) == 0 )  return false;

    int iterations = 0; const int MAX_ITERATIONS = 1000; // prevent infinite loops

    while( particleID > ancestorID && iterations < MAX_ITERATIONS ) {
      iterations++;
      if( map_g4id_mother.count(particleID) == 0 ) { return false; }

      int mother = map_g4id_mother[particleID];
      if( mother == ancestorID ) { return true; } if( mother == 0 ) { return false; }
      //if( !map_g4id_pdg.count(mother) ) { return false; }
      //if ( breakAtPhots == true &&  map_g4id_pdg[mother] == 22 )      { return false; }
      particleID = mother;
    }
    //if( iterations >= MAX_ITERATIONS ) { std::cout << "\nWarning: Hit max iterations in IsAncestorOf - infinite loop?!"; } 
    return false;
  }

  TFile* f = new TFile("accessory_macros/fit_cut.root", "READ"); 
  TGraph* lowCutoff_graph = (TGraph*)f->Get("lowCutoff_graph");
  TGraph* highCutoff_graph = (TGraph*)f->Get("highCutoff_graph");
  TTree* param_tree = (TTree*)f->Get("param_tree");
 
  float highA, highB, highC, highD, lowA, lowB, lowC, lowD, lowE;
  param_tree->SetBranchAddress("highParam_A", &highA);
  param_tree->SetBranchAddress("highParam_B", &highB);
  param_tree->SetBranchAddress("highParam_C", &highC);
  param_tree->SetBranchAddress("highParam_D", &highD);
  param_tree->SetBranchAddress("lowParam_A", &lowA);
  param_tree->SetBranchAddress("lowParam_B", &lowB);
  param_tree->SetBranchAddress("lowParam_C", &lowC);
  param_tree->SetBranchAddress("lowParam_D", &lowD);
  param_tree->SetBranchAddress("lowParam_E", &lowE);
  param_tree->GetEntry(0);  // Only one entry with the parameters

  f->Close();

  float Fit_Function_PID_Cut_MC( float blip_recoE ){
    float x = blip_recoE;
    //if( blip_recoE > 3 )
    return highA * TMath::Log( highB * x ) + highC*x + highD ;
    //else return -lowA*x*x*x + lowB*x*x + lowC*x + lowD/x + lowE;
  }

  //============================== EDEP FUNCTION ==============================//
  void edep_blip_eff_plotter() {   // if you reimplement this, NEED TO ADD BRANCH ADDRESS SETTING
    initialize();
    h_blipRecoEff_vs_vtxDist = new TH1D("blipRecoEff_vs_vtxDist","Blip Distance to Reco Neutrino Vertex;Distance [cm];Nblips/Nparticles",101,-1,201);
    TH1D* h_edeps_vs_vtxDist = (TH1D*)h_blipRecoEff_vs_vtxDist->Clone("edeps_vs_vtxDist");   // temporary hist, not globally declared
    h_blipRecoEff_energy_dist= new TH2D("blipRecoEff_energy_dist","Blip Energy vs Distance to Reco Neutrino Vertex;Distance [cm];Energy [MeV]",40,0,200,100,0,10);
    TH2D* h_edep_energy_dist = (TH2D*)h_blipRecoEff_energy_dist->Clone("edep_energy_dist");
    TVector3 edep_pos; int edepDist; float edepE;
    if( maxEvents < 0 ) maxEvents = numEvents;
    for( int iEvent = 0; iEvent < maxEvents; iEvent++ ) {
      tree->GetEntry(iEvent); if( mctruth_nu_ccnc != 0 || mctruth_nu_pdg != 14 ) continue; 
      TVector3 truth_nu_vtx( mctruth_nu_vtx_x, mctruth_nu_vtx_y, mctruth_nu_vtx_z );
      for( int i=0;i<nedeps;i++ ) {
	edep_pos.SetXYZ( edep_x[i], edep_y[i], edep_z[i] ); edepDist = (edep_pos-truth_nu_vtx).Mag(); edepE = edep_energy[i];
	if( !edep_isPrimary[i] ) { h_edeps_vs_vtxDist->Fill(edepDist); h_edep_energy_dist->Fill(edepDist, edepE);
	  if( edep_blipid[i] >= 0 ) { h_blipRecoEff_vs_vtxDist->Fill(edepDist); h_blipRecoEff_energy_dist->Fill(edepDist, edepE); }
	}
      } //end edep loop
    }
    h_blipRecoEff_vs_vtxDist->Divide( h_edeps_vs_vtxDist ); delete h_edeps_vs_vtxDist;
    h_blipRecoEff_energy_dist->Divide( h_edep_energy_dist ); delete h_edep_energy_dist;
    h_blipRecoEff_energy_dist->GetZaxis()->SetRangeUser(0, 1); h_blipRecoEff_energy_dist->SetOption("colz");
    outfile->Write(); outfile->Close();
  }

  void BACKUP_neutrino_analysis(){

    // Create output file and tree
    //TFile* remainFile = new TFile("remaining_blips.root", "RECREATE");
    //TTree* outTree = new TTree("what_remains", "What Remains When it's All Over");

    // load the shared library
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

    initialize();

    // Make histograms declared above
    h_all_nu_neutrons	     = new TH1D("all_nu_neutrons","Number of neutrino-induced neutrons; # neutrons",23,0,23);
    h_trueNeutronE	     = new TH1D("trueNeutronE","True Neutron Energy;MeVee",60,0,600);h_trueNeutronE->Sumw2();
    h_preFoundNeutronE	     = new TH1D("preFoundNeutronE","Preselection Identified Neutron Energy;MeVee",60,0,600);h_preFoundNeutronE->Sumw2();
    h_foundNeutronE	     = new TH1D("foundNeutronE","Identified Neutron Energy;MeVee",60,0,600);h_foundNeutronE->Sumw2();
    h_NC_Neutrons    	     = new TH1D("NC_Neutrons","Number of NC neutrons per neutrino; # neutrons",15,0,15);
    h_CC_Neutrons     	     = new TH1D("CC_Neutrons","Number of CC neutrons per neutrino; # neutrons",15,0,15);
    h_allPrBlips_vs_vtxDist  = new TH1D("allPrBlips_vs_vtxDist","All Proton Blips Dist to Reco Vtx;Dist to reco vtx [cm]; Nblips / r^2",26,-1,51);
    //h_truePr_vs_vtxDist    = new TH1D("truePr_vs_vtxDist","Proton Dist to Reco Vertex;Dist to reco vtx [cm]; Nblips / r^2",104,-1,51);
    h_truePrCosAlign	     = new TH1D("truePrCosAlign","(n, 1p) Cos(Angle) B/w Reco Inferred Neutron Direction and True Momentum;cos(angle)",50,0,1);
    h_maxWires		     = new TH1D("maxWires","Max Wires per Blip on one Plane;# Wires",18,0,18);
    h_blipZY                 = new TH2D("blipZY","Blip Z vs Y;z;y",224,-40,1080,50,-150,150);
    h_numTrks                = new TH1D("numTrks","Number of Tracks;# tracks",20,0,20);
    h_trackLength	     = new TH1D("trackLengths","Track Lengths;length [cm]",150,-10,300);
    h_prTrksPerE	     = new TH1D("prTrksPerE","Number of Protons Made Into Tracks vs Energy;True proton kinetic energy [MeV]",100,0,200);
    h_prBlipsPerE            = new TH1D("prBlipsPerE","Number of Protons Made Into Blips vs Energy;True proton kinetic energy [MeV]",100,0,200);
    //h_trueProtonE	     = new TH1D("trueProtonE","Neutron Inelastic Proton Energy; E (MeV)",204,-1,101);
    //h_truePrimPrE            = new TH1D("truePrimPrE","Energy of Protons with Primary Neutron Mothers; E (MeV)",204,-1,101);
    h_trueBlipE		     = new TH1D("trueBlipE","True Blip Energy; E (MeV)",204,-1,101);
    h_neutronTrueE           = new TH1D("neutronTrueE","Neutron Descended Blip True Energy; E (MeV)",204,-1,101);
    h_gammaTrueE             = new TH1D("gammaTrueE","Gamma Blip True Energy; E (MeV)",204,-1,101);
    h_signalBlipE            = new TH1D("signalBlipE","Signal Proton Blip from Primary Neutron Reco Energy; E [MeVee]",110,0,55);
    h_trueNeutronEvsProtonE  = new TH2D("trueNeutronEvsProtonE","Proton Energy vs Neutron Energy;Neutron kinetic energy [MeV];FS proton kinetic energy [MeV]",60,0,300,60,0,300);
    h_trueRecoVtxAgreement   = new TH1D("trueRecoVtxAgreement","Difference Between True and Reco Vertices;Distance [cm]",204,-1,101);
    // cos of the angle b/w  proton direction and inferred neutron direction (vector between nu vtx and proton location)
    //h_cos_prRecoDirection  = new TH1D("cos_prRecoDirection","Cos(angle) b/w Reco. e- Direction and Proton Direction;cos(ang);# blips",51,-1.01,1.01);
    //h_cos_prComparison     = new TH1D("cos_prComparison","Initial e- Momentum Cos Minus e- Direction Cos;Difference in Cos(ang);# blips",20,0,2);
    h_cos_vs_energy	     = new TH2D("cos_vs_energy","True Cos(angle) vs True Energy;E [MeV];cos(ang)",200,0,200,102,-1.02,1.02);
    h_pr_cos_vs_energy       = new TH2D("pr_cos_vs_energy","True Cos(angle) vs Proton Energy;proton E [MeV];cos(ang)",200,0,200,102,-1.02,1.02);
    h_pr_E_ds		     = new TH2D("pr_E_ds","Proton Blip Energy Deposition Density;Reco blip E [MeVee]; Eblip/ds [MeV/cm]",500,0,10,250,0,5);//150,0,15,120,0,12);
    gr_signal_E_ds = new TGraph(); gr_signal_E_ds->SetNameTitle("signal_E_ds", "Proton Blip Energy Deposition Density;Reco blip E [MeVee];Eblip/ds [MeV/cm]");
    h_elec_E_ds              = new TH2D("elec_E_ds","Electron Blip Energy Deposition Density;Reco blip E [MeVee]; Eblip/ds [MeV/cm]",500,0,10,250,0,5);
    gr_elec_E_ds = new TGraph(); gr_elec_E_ds->SetNameTitle("gr_elec_E_ds", "Electron Blip Energy Deposition Density;Reco blip E [MeVee];Eblip/ds [MeV/cm]");
    gr_cosmic_E_ds = new TGraph(); gr_cosmic_E_ds->SetNameTitle("cosmic_E_ds", "Cosmic Blip Energy Deposition Density;Reco blip E [MeVee];Eblip/ds [MeV/cm]");
    h_Np_vs_Nn		     = new TH2D("Np_vs_Nn", "Proton Daughters vs Neutron Daughters ;Nn;Np",16,0,16,12,0,12);
    h_res_vs_trueE	     = new TH2D("res_vs_trueE","Angular Resolution vs True Blip Energy;blipE;Angle [degrees]",200,0,200,90,0,90);
    h_pr_res_vs_trueE	     = new TH2D("pr_res_vs_trueE","(n, 1p) Angular Resolution vs True Proton Kinetic Energy;True Proton KE [MeV];Angle [degrees]",50,0,100,45,0,90);
    gr_pr_res_vs_trueE = new TGraph(); gr_pr_res_vs_trueE->SetNameTitle("gr_pr_res_vs_trueE", "Signal Proton Angular Resolution vs True Blip Energy;blipE;Angle [degrees]");
    blipE_plots->cd();
    h_blipE[0]		     = new TH1D("blipE","Reconstructed Blip Energy;Energy [MeVee]",30,0,90);//20,0,100);
    for(int i=1; i<numHists; i++) { h_blipE[i] = (TH1D*)h_blipE[0]->Clone(("blipE_" + std::to_string(i)).c_str()); }
    recoVtxDist_plots->cd();
    h_recoVtxDist[0]	     = new TH1D("recoVtxDist","Blip Distance to Reconstructed Neutrino Vertex;Distance [cm];Nblips",20,0,20); //20,0,100);
    for(int i=1; i<numHists; i++) { h_recoVtxDist[i] = (TH1D*)h_recoVtxDist[0]->Clone(("recoVtxDist_" + std::to_string(i)).c_str()); }
    cosAlignment_plots->cd();
    h_cosAlignment[0]        = new TH1D("cosAlignment","Cos(Alignment Angle) Between Reconstructed Direction and Neutrino Vertex;Cos",20,0,1);//40,-1,1);
    for(int i=1; i<numHists; i++) { h_cosAlignment[i] = (TH1D*)h_cosAlignment[0]->Clone(("cosAlignment_" + std::to_string(i)).c_str()); }
    alignmentAngle_plots->cd();
    h_alignmentAngle[0]      = new TH1D("alignmentAngle","Alignment Angle Between Reconstructed Direction and Neutrino Vertex;Cos",30,0,90);
    for(int i=1; i<numHists; i++) { h_alignmentAngle[i] = (TH1D*)h_alignmentAngle[0]->Clone(("alignmentAngle_" + std::to_string(i)).c_str()); }
    //res_vs_E_plots->cd();
    //h_res_vs_E[0]	     = new TH2D("res_vs_E","Angular Resolution vs. Blip Reco Energy;reco E [MeVee];angle [degrees]",20,0,100,30,0,90);
    //for(int i=1; i<numHists; i++) { h_res_vs_E[i] = (TH2D*)h_res_vs_E[0]->Clone(("res_vs_E_" + std::to_string(i)).c_str()); }
    resolution_plots->cd();
    h_resolution[0]          = new TH1D("resolution","Blip Angular Resolution;angle [degrees]",30,0,90);//46,-1,91);
    for(int i=1; i<numHists; i++) { h_resolution[i] = (TH1D*)h_resolution[0]->Clone(("resolution_" + std::to_string(i)).c_str()); }
    Edensity_plots->cd();
    h_Edensity[0]	     = new TH1D("Edensity","Blip Energy Deposition Density;E/ds [MeV/cm]",30,0,15);
    for(int i=1; i<numHists; i++) { h_Edensity[i] = (TH1D*)h_Edensity[0]->Clone(("Edensity_" + std::to_string(i)).c_str()); }
    proxTrkDist_plots->cd();
    h_proxTrkDist[0]         = new TH1D("proxTrkDist","Blip Distance to Nearest Track;Reconstructed Distance [cm]",21,0,105);   //28,-1,55);   //
    for(int i=1; i<numHists; i++) { h_proxTrkDist[i] = (TH1D*)h_proxTrkDist[0]->Clone(("proxTrkDist_" + std::to_string(i)).c_str()); }
    byDeadWire_plots->cd();
    h_byDeadWire[0]          = new TH1D("byDeadWire","Blip By Dead Wire",2,0,2);
    for(int i=1; i<numHists; i++) { h_byDeadWire[i] = (TH1D*)h_byDeadWire[0]->Clone(("byDeadWire" + std::to_string(i)).c_str()); }
    cosmicDistToVtx_plots->cd();
    h_cosmicDistToVtx[0]      = new TH1D("cosmicDistToVtx","Reco. Neutrino Vertex Distance to Nearest Cosmic Track;Reconstructed Distance [cm]",41,0,510);
    for(int i=1; i<numHists; i++) { h_cosmicDistToVtx[i] = (TH1D*)h_cosmicDistToVtx[0]->Clone(("cosmicDistToVtx_" + std::to_string(i)).c_str()); }
    dx_plots->cd();
    h_dx[0]		      = new TH1D("dx","Reconstructed Blip Width in X;dx [cm]",30,0,9);
    for(int i=1; i<numHists; i++) { h_dx[i] = (TH1D*)h_dx[0]->Clone(("dx_" + std::to_string(i)).c_str()); }
    dw_plots->cd();
    h_dw[0]		      = new TH1D("dw","Reconstructed Blip Width Across Wires;dw [cm]",22,0,6.6);
    for(int i=1; i<numHists; i++) { h_dw[i] = (TH1D*)h_dw[0]->Clone(("dw_" + std::to_string(i)).c_str()); }
    size_plots->cd();
    h_size[0]		      = new TH1D("size","Reconstructed Blip Size;Size [cm]",30,0,9);
    for(int i=1; i<numHists; i++) { h_size[i] = (TH1D*)h_size[0]->Clone(("size_" + std::to_string(i)).c_str()); }
    Nhits_plots->cd();
    h_Nhits[0]                 = new TH1D("Nhits","Number of Hits Per Blip;Nhits",55,0,55);
    for(int i=1; i<numHists; i++) { h_Nhits[i] = (TH1D*)h_Nhits[0]->Clone(("Nhits_" + std::to_string(i)).c_str()); }
    Nhits_plots->cd();
   // h_Nhits[0]                 = new TH1D("Nhits","Number of Hits Per Blip;Nhits",55,0,55);
   // for(int i=1; i<numHists; i++) { h_Nhits[i] = (TH1D*)h_Nhits[0]->Clone(("Nhits_" + std::to_string(i)).c_str()); }
   // outfile->cd();

    // seems to be adding only to TTree?? (outfile)
    h_Np_vs_Nn->SetOption("colz text");
    h_Np_vs_Nn->SetMarkerSize(1.5);  // Text size
    h_Np_vs_Nn->SetMarkerColor(kBlack); // Text color

    int   total_nus=0;
    float total_Ns=0;
    float total_nu_Ns=0;
    int   total_N_prs=0;
    int   totalBlips=0;
    float totalSelectedBlips=0;
    int   trueBlips=0;
    int   total_N_blips=0;
    int   total_gammas=0;
    int   total_gamma_blips=0;
    int   total_N_gammas=0;   // these (and below) are gamma BLIPS
    int   total_nCap_gammas=0;
    int   total_electrons=0;
    float n_p_blips=0;
    float signal_blips=0;
    float single_pr=0;
    float num_n_p=0;
    int otherNotGamma=0;
    int otherIsGamma=0;
    int   exitWoutScat=0; int   hasScat=0;
    int   getsCaptured=0; int   Nother=0;
    // do I need to declare these here?
    int   p, p_g4id, p_pdg, mom_g4id, mom_idx, grandmom_g4id, grandmom_idx;
    float trueE, blipE;
    int   category;
    std::vector<float> cos_allBlips;
    std::vector<float> pr_energy;
    //ProcessHandler handler( h_cos_nIne_gamma, h_cos_nCap_gamma, h_cos_eBrem_gamma, h_cos_nOther_gamma, total_nCap_gammas );

    //============================== EVENT LOOP ==============================//
    if( maxEvents < 0 ) maxEvents = numEvents;
    for( int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      //reset_vars();   //clear_branch_vars();
      tree->GetEntry(iEvent);
      //if( mctruth_nu_ccnc != 0 || mctruth_nu_pdg != 14 ) continue; 
      total_nus++;   // one neutrino per event

      // set blipID map
      std::map<int,std::vector<int>> map_blipID_hitID;
      for(int ih=0; ih<nhits; ih++){ if(hit_blipid[ih] >= 0) { map_blipID_hitID[hit_blipid[ih]].push_back(ih); } }

      std::unordered_set<int> neutron_g4id; std::unordered_set<int> n_p_g4id_set;
      std::vector<int> ncapture_products_g4id; //std::vector<VisParticle> particles_for_vis;
      std::vector<int> g4ids_for_vis;   std::vector<int> otherPrimaryG4ID;
      float selectedBlips=0;
      int   contained_nu_Ns=0;
      int   nu_Ns=0;
      int   NC_Ns=0; int   CC_Ns=0;
      bool  create_display = false;   // criteria for creating 3D visualization

      TVector3 reco_nu_vtx( nu_reco_vtx_x, nu_reco_vtx_y, nu_reco_vtx_z );
      TVector3 truth_nu_vtx( mctruth_nu_vtx_x, mctruth_nu_vtx_y, mctruth_nu_vtx_z );
      float cosmicProx = nu_crt_closestCosDist; bool nearCosmic=false;
      //if( cosmicProx < 20 ) nearCosmic=true;

      map_g4id_mother.clear(); map_g4id_pdg.clear();
      map_g4id_index.clear(); neutron_g4id.clear();
      map_nu_N_g4id_daught_idx.clear();
      map_n_p.clear(); map_n_n.clear();
      n_p_g4id_set.clear();

      //============================== PARTICLE LOOP ==============================//
      for(int i=0; i<nparticles; i++){
	int pdg = part_pdg[i], g4id = part_g4id[i], mom = part_mother[i];
	map_g4id_mother[g4id] = mom; map_g4id_pdg[g4id] = pdg; map_g4id_index[g4id] = i;

	TVector3 part_loc( part_startPointx[i], part_startPointy[i], part_startPointz[i] );
	float dist = (part_loc-reco_nu_vtx).Mag(); float trueDist = (part_loc-truth_nu_vtx).Mag();
	TVector3 part_P( part_Px[i], part_Py[i], part_Pz[i] ); TVector3 trace_vec( part_loc-truth_nu_vtx );
        float ang_blip = part_P.Angle(trace_vec); float cos = TMath::Cos(ang_blip);

        // Printout this particle information to the screen 
        if( 0 && iEvent==2251 ){   //&& part_isContained[i]
          printf("%4i G4ID: %-4i PDG: %-10i dL=%6.1f cm, KE0=%7.3f MeV, Edep=%7.3f MeV, moth=%5i, %16s, NDaughts=%3i, start(X,Y,Z)=(%1.2f, %1.2f, %1.2f)\n",
            i, part_g4id[i], part_pdg[i], part_pathlen[i], part_KE[i], part_depEnergy[i], part_mother[i],
	    part_process->at(i).c_str(), part_nDaughters[i], part_loc.X(), part_loc.Y(), part_loc.Z()
          );
        }

        // -------------------------------------------------------------
        // Identify any neutrons and save their information
        if ( part_pdg[i]==2112 ) {
	  total_Ns++; neutron_g4id.insert(g4id);
	  if( part_isPrimary[i] ) {   // primary neutrons
	    total_nu_Ns++; nu_Ns++; h_trueNeutronE->Fill(part_KE[i]);
	    map_n_p[g4id].clear(); map_n_n[g4id].clear();
	  }
	  if(mctruth_nu_ccnc) NC_Ns++; else CC_Ns++;
	}  //if( part_isContained[i] ) { contained_nu_Ns++; }

	if( part_isPrimary[i] && pdg != 2112 && pdg != 2212 ) { otherPrimaryG4ID.push_back(g4id); }

	if( part_process->at(i) == "nCapture" ) { ncapture_products_g4id.push_back(g4id); }

        if( !part_isPrimary[i] && map_n_p.count(mom) ) {   // from primary neutron mothers
	  map_nu_N_g4id_daught_idx[mom] = i;
          if( pdg == 2212 ) { num_n_p++; n_p_g4id_set.insert(g4id); map_n_p[mom].push_back(g4id);   // protons
	    h_trueNeutronEvsProtonE->Fill(part_KE[map_g4id_index[mom]],part_KE[i]); }
          if( pdg == 2112 ) map_n_n[mom].push_back(g4id);   // neutrons
        }
	if( pdg == 2212 && map_g4id_pdg[mom] == 2112 ) total_N_prs++;   // may cause future indexing errors

        // gammas (PDG=22), electrons (PDG=11)
        if( part_pdg[i] == 22 ) { total_gammas++; continue; } if( part_pdg[i] == 11 ) { total_electrons++; continue; }

      }//end loop over all MCParticles


      //============================== TRACK LOOP ==============================//
      h_NC_Neutrons->Fill(NC_Ns); h_CC_Neutrons->Fill(CC_Ns); h_all_nu_neutrons->Fill(nu_Ns);
      h_numTrks->Fill(ntrks); for( int i=0; i<ntrks; i++ ) {
	h_trackLength->Fill(trk_length[i]);
	if( n_p_g4id_set.count(trk_g4id[i]) ) h_prTrksPerE->Fill(part_KE[map_g4id_index[trk_g4id[i]]]);
      }//end track loop


      if(0) {
	exitWoutScat += nu_Ns - map_nu_N_g4id_daught_idx.size();
	for ( const auto& [Ng4id, daughtIDX] : map_nu_N_g4id_daught_idx ) {
	  bool capture = false; float daughtX = part_startPointx[daughtIDX],
				      daughtY = part_startPointy[daughtIDX],
				      daughtZ = part_startPointz[daughtIDX];
	  if( 21.50 <= daughtX && daughtX <= 234.85 && -95.0 <= daughtY && daughtY <= 95.00 && 21.50 <= daughtZ && daughtZ <= 966.80 ) {   // if mother && daughter inside FV (defined here: MICROBOONE-NOTE-1099-PUB)
	    for( int i=0; i<ncapture_products_g4id.size(); i++ ) {
	      if( IsAncestorOf( ncapture_products_g4id[i], Ng4id, false ) ) { capture=true; getsCaptured++; break; }
	    }
	    if( capture == false ) {
	      if( part_process->at(daughtIDX) == "neutronInelastic" ) hasScat++;
	      else Nother++;
	    } 
	  } else exitWoutScat++;
	} }

      //============================== BLIP LOOP ==============================//
      for(int i=0; i<nblips; i++) {
	//std::cout<<"Edepid: "<<blip_edepid[i]<<", p_index: "<<map_g4id_index[blip_g4id[i]]<<"\n";
	//if( blip_edepid[i] < 0 ) continue;
	//if( map_blipID_hitID.find(blip_id[i]) == map_blipID_hitID.end() ) continue;

        // RNN model selection, if added back, need to put Nhit business up here
	//std::vector<float> recoDirComponents = blipDirection(
        //  blip_id[i], map_blipID_hitID, file, hit_rms, hit_amp,
        //  hit_wire, hit_plane, hit_driftT, model_path
        //);
	//TVector3 reco_dir( recoDirComponents[0], recoDirComponents[1], recoDirComponents[2] ); 
	//if( reco_dir.X() == -9 ) continue;
	//if( map_g4id_index.count(blip_g4id[i]) > 0 ) {
	//  p = map_g4id_index.at(blip_g4id[i]);
        //  //if( !(part_pdg[p] == 2212) || !part_isContained[p] || !part_isPrimary[p] ) continue;
	//  TVector3 true_P(part_Px[p], part_Py[p], part_Pz[p]);
	//  float res_angle = reco_dir.Angle(true_P)*TMath::RadToDeg();
	//  if( res_angle > 90 ) res_angle = 180-res_angle;
	//  h_pr_res_vs_trueE->Fill(part_KE[p], res_angle);
	//} else p=-9;
	//totalBlips++; continue;

	if(1){ 
	category=-9;   totalBlips++;
	
	float blipX=blip_x[i], blipY=blip_y[i], blipZ=blip_z[i];
	TVector3 blip_loc( blipX, blipY, blipZ ); TVector3 reco_trace_vec( blip_loc - reco_nu_vtx );
	bool isProton = true; int Nhits = map_blipID_hitID.at(blip_id[i]).size();
	float dist = (blip_loc-reco_nu_vtx).Mag(); float dist2Trk = blip_proxtrkdist[i]; bool nearDeath=false;
	float blipE = blip_energy[i]; float density = blipE/blip_size[i]; float touchTrk = blip_touchtrk[i];
	float vtxAlign = abs((reco_nu_vtx-truth_nu_vtx).Mag()); //int maxWires = std::round(blip_dw[i]/0.3); h_maxWires->Fill(maxWires);   // in cm, convert to Nwires by /0.3 cm
	h_trueRecoVtxAgreement->Fill(vtxAlign);

	// setting particle index & vars ( may cause trouble if not truth matched? )
	if( map_g4id_index.count(blip_g4id[i]) > 0 ) p = map_g4id_index[blip_g4id[i]]; else p = -9;
	bool isPrimary = part_isPrimary[p];
	p_g4id = part_g4id[p]; p_pdg = part_pdg[p]; trueE = part_KE[p];
	if ( isPrimary ) { mom_g4id = mom_idx = -9; }
	else { mom_g4id = part_mother[p]; mom_idx = map_g4id_index[mom_g4id];
	  if( part_isPrimary[mom_idx] ) { grandmom_g4id = grandmom_idx = -9; }
	  else { grandmom_g4id = part_mother[mom_idx]; grandmom_idx = map_g4id_index[grandmom_g4id]; }
	}
	TVector3 true_P(part_Px[p], part_Py[p], part_Pz[p]);

        if( blip_pl2_bydeadwire[i] != 1 ) { if ( blip_pl2_bydeadwire[i] != 0 ) std::cout<<"Blip by dead wire (int) does not equal 0 or 1 ! It equals: "<<blip_pl2_bydeadwire[i]<<std::endl; }

	if( !isPrimary && map_n_p.count(mom_g4id) && p_pdg == 2212 ) h_prBlipsPerE->Fill(trueE);

	//============================== CUTS ==============================//

	if( touchTrk || dist2Trk <= 1 ) continue;   // || vtxAlign >= 1.5 
	if( blip_pl0_bydeadwire[i] == 1 || blip_pl1_bydeadwire[i] == 1 || blip_pl2_bydeadwire[i] == 1 ) nearDeath=true;
	if( dist >= 100  ) continue;
	if( p_pdg == 2212 && part_pdg[mom_idx] == 2112 && part_isPrimary[mom_idx] && !nearCosmic )   // signal blips
	  { h_preFoundNeutronE->Fill(part_KE[mom_idx]); }

        if( dist <= 5 ) continue;
	if( blipE <= 3 ) continue;
	if( density < Fit_Function_PID_Cut_MC( blipE )) isProton=false; //if( density < cutoff_graph->Eval(blipE) ) continue;
	//if( map_blipID_hitID.count(blip_id[i]) == 0 || Nhits <= 2 ) continue;   // hits cut, only necessary if blipE < 3

	double res_angle = -1.0; TVector3 reco_dir; float align_angle, offset;
	if(1) {
          // Call the blipDirection function from the library
          std::vector<float> recoDirComponents = blipDirection(
	    blip_id[i], map_blipID_hitID, file, hit_rms, hit_amp,
	    hit_wire, hit_plane, hit_driftT, model_path
          );

	  // RNN model selection
	  reco_dir.SetXYZ( recoDirComponents[0], recoDirComponents[1], recoDirComponents[2] );
	  align_angle = reco_dir.Angle(reco_trace_vec); offset = abs(TMath::Cos(align_angle));   // pointing offset
	  align_angle = align_angle*TMath::RadToDeg(); res_angle = reco_dir.Angle(true_P)*TMath::RadToDeg();
	  if( align_angle > 90 ) align_angle = 180-align_angle;
          if( res_angle > 90 ) res_angle = 180-res_angle;

	  //if( offset < 0.25 ) continue;

	  if( isProton && !nearDeath && !nearCosmic ) {
	    h_cosAlignment[0]->Fill(offset); h_alignmentAngle[0]->Fill(align_angle);
	    if( blip_edepid[i] < 0 ) { h_cosAlignment[9]->Fill(offset); h_alignmentAngle[9]->Fill(align_angle); }
	    else { h_resolution[0]->Fill(res_angle); h_res_vs_trueE->Fill(trueE, res_angle); } //h_res_vs_E[0]->Fill(blipE, res_angle); }
	  }
	}

	if( isProton ) {
	  if( !nearDeath ) {
	    if( !nearCosmic ) {
	      selectedBlips++; h_blipE[0]->Fill(blipE); h_recoVtxDist[0]->Fill(dist); h_Edensity[0]->Fill(density); h_proxTrkDist[0]->Fill(dist2Trk);
	      h_dx[0]->Fill(blip_dx[i]); h_dw[0]->Fill(blip_dw[i]); h_size[0]->Fill(blip_size[i]); h_Nhits[0]->Fill(Nhits);
	      //set_branch_vars( isOutTree, i, map_g4id_index, 0, 0 );
	    } h_cosmicDistToVtx[0]->Fill(cosmicProx);
	  } h_byDeadWire[0]->Fill(nearDeath);
	}

	// cosmic separation
	if( blip_edepid[i] < 0 ) {
	  if( !nearDeath ) {
	    if( blipE <= 10 ) { gr_cosmic_E_ds->SetPoint(gr_cosmic_E_ds->GetN(), blipE, density); }
	    if( isProton ) {
	      if( !nearCosmic ) {
		h_blipE[9]->Fill(blipE); h_recoVtxDist[9]->Fill(dist);
		h_Edensity[9]->Fill(density); h_proxTrkDist[9]->Fill(dist2Trk);
 		h_dx[9]->Fill(blip_dx[i]); h_dw[9]->Fill(blip_dw[i]);
		h_size[9]->Fill(blip_size[i]); h_Nhits[9]->Fill(Nhits);
	      } h_cosmicDistToVtx[9]->Fill(cosmicProx);
	    }
	  } if( isProton ) { h_byDeadWire[9]->Fill(nearDeath); } continue;
	} trueBlips++;

	TVector3 part_loc(part_startPointx[p], part_startPointy[p], part_startPointz[p]); //std::cout<<"|part_loc-blip_loc|: "<<(part_loc-blip_loc).Mag()<<"\n";
	TVector3 true_trace_vec( part_loc - truth_nu_vtx ); float true_ang_blip = true_P.Angle(true_trace_vec);
	float true_cos = TMath::Cos(true_ang_blip); cos_allBlips.push_back(true_cos); h_cos_vs_energy->Fill(trueE, true_cos);
	float reco_ang_blip = true_P.Angle(reco_trace_vec); float reco_cos = TMath::Cos(reco_ang_blip);

	if( p_pdg == 2212 ) { h_pr_E_ds->Fill(blipE, density); if( part_isPrimary[p] ) { category=4; } }   // primary protons

	bool is_n_descent = std::any_of(neutron_g4id.begin(), neutron_g4id.end(), [&](int n_g4id) { return IsAncestorOf(p_g4id, n_g4id); });

	// NEUTRON descended blips
	if( is_n_descent ) {
	  // add up the total (primary) neutron-induced blips
          if( abs(p_pdg) == 11 && part_pdg[mom_idx] == 22 && map_g4id_pdg[part_mother[mom_idx]] == 2112 && part_isPrimary[map_g4id_index[part_mother[mom_idx]]] )
	  { total_N_blips++; } else if ( part_pdg[mom_idx] == 2112 && part_isPrimary[mom_idx] ) total_N_blips++;
	  // all other neutron-induced blips
          if( part_pdg[mom_idx] == 2112 && p_pdg != 2212 ) { h_neutronTrueE->Fill(trueE); category=6; }

	  //if( p_pdg == 2212 && part_isPrimary[p] && map_n_p.count(mom_g4id) && map_n_p.at(mom_g4id).size() == 1 )
	  //{ h_pr_res_vs_trueE->Fill(trueE, res_angle); }

	  // PROTONS
 	  if( p_pdg == 2212 && !part_isPrimary[p] ) {   // any proton
	    if( part_pdg[mom_idx] == 2112 ) {   // from any neutron mom (should all be nIne)
	      n_p_blips++; h_allPrBlips_vs_vtxDist->Fill( dist );
	      if( part_isPrimary[mom_idx] ) { // *primary* neutron mom
		if( blipE <= 10 ) gr_signal_E_ds->SetPoint(gr_signal_E_ds->GetN(), blipE, density);
		if( !nearCosmic && isProton ) { if( !nearDeath ) {
		    signal_blips++; h_foundNeutronE->Fill(part_KE[mom_idx]);
		    h_signalBlipE->Fill(blipE); 
		  } //else { std::cout<<run<<" "<<subrun<<" "<<iEvent<<std::endl; }
		}
		if( map_n_p.count(mom_g4id) && map_n_p.at(mom_g4id).size() == 1 ) {   // single protons only
		  if( !nearCosmic && !nearDeath && isProton ) { single_pr++;
		    h_truePrCosAlign->Fill(true_cos); h_pr_res_vs_trueE->Fill(trueE, res_angle); }
		  category=1;   // 1p from primary neutron moms
		  h_pr_cos_vs_energy->Fill(trueE, true_cos); pr_energy.push_back(trueE);
		  //if( true_cos < 0.2 && trueE > 60 ) std::cout<<run<<" "<<subrun<<" "<<iEvent<<std::endl;//create_display = true;
		  //g4ids_for_vis.push_back(mom_g4id); g4ids_for_vis.push_back(p_g4id);
		} else { category=2; }   // Np from primary neutron moms
	      } else { category=3; }   // non-primary neutron-induced protons
	    }
	  }

	  // GAMMAS (neutron-descended)
	  if( abs(p_pdg) == 11 ) { if( part_pdg[mom_idx] == 22 && part_pdg[grandmom_idx] == 2112 ) { total_N_gammas++; category=5; } }

	}//end neutron check

	// ELECTRONS
	if( abs(p_pdg) == 11 ) { h_elec_E_ds->Fill(blipE, density); if( blipE <= 10 ) gr_elec_E_ds->SetPoint(gr_elec_E_ds->GetN(), blipE, density);
	  if( part_pdg[mom_idx] == 22 && map_g4id_pdg[part_mother[mom_idx]] == 2112 ) { h_gammaTrueE->Fill(trueE); }
	}

	// all "OTHER" (non-neutron induced) particles
	if( part_pdg[mom_idx] != 2112 && !(abs(p_pdg) == 11 && part_pdg[mom_idx] == 22) && !(part_isPrimary[p] && p_pdg == 2212) ) {
	  category=7;
	  if( isProton && !nearDeath && !nearCosmic )
            { otherNotGamma++; }//printf("NON-gamma 'other' process: %20s     PDG: %-10i\n", part_process->at(p).c_str(), part_pdg[p]); }
        }
	// GAMMAS
	else if( abs(p_pdg) == 11 && part_pdg[mom_idx] == 22 ) { total_gamma_blips++;
	  if( part_pdg[grandmom_idx] != 2112 ) {
	    if( part_process->at(mom_idx) == "eBrem" ) { category=8; }   //std::cout<<"BREM\n";
	    else {
	      category=7;
	      if( isProton && !nearDeath && !nearCosmic ) { otherIsGamma++; }
		//printf("'Other' gamma process: %24s     Mother PDG: %-10i    Event: %i     G4ID: %i\n", part_process->at(mom_idx).c_str(), part_pdg[grandmom_idx], iEvent, mom_g4id);
	    }
	  }
	}

	for( int oG4ID : otherPrimaryG4ID ) {
	  if( IsAncestorOf(p_g4id, oG4ID, false) && category != 3 && category != 5 && category != 6 ) { category=8; }
	}

	if( isProton ) {
	  if( !nearDeath ) {
	    if( !nearCosmic ) {
	    if( category == -9 ) { std::cout<<"\nOH NAUR!, naur catagaury set ;(\nPDG = "<<p_pdg<<std::endl; continue; }
	      h_Edensity[category]->Fill(density); h_proxTrkDist[category]->Fill(dist2Trk);
	      h_recoVtxDist[category]->Fill(dist); h_blipE[category]->Fill(blipE);
	      h_dx[category]->Fill(blip_dx[i]); h_dw[category]->Fill(blip_dw[i]);
	      h_size[category]->Fill(blip_size[i]); h_Nhits[category]->Fill(Nhits);
	      h_cosAlignment[category]->Fill(offset); //h_res_vs_E[category]->Fill(blipE, res_angle);
	      h_alignmentAngle[category]->Fill(align_angle); h_resolution[category]->Fill(res_angle);
	    } h_cosmicDistToVtx[category]->Fill(cosmicProx);
	  } h_byDeadWire[category]->Fill(nearDeath);
	}

	h_trueBlipE->Fill(trueE); h_blipZY->Fill(blipZ,blipY);
	//outTree->Fill();

      }}//end blip loop

      // fill Np vs Nn plot
      for (const auto& n_pair : map_n_n) {
	int n_g4id = n_pair.first; int Nn = n_pair.second.size();
	auto p_it = map_n_p.find(n_g4id);
	if (p_it != map_n_p.end()) {
	  int Np = p_it->second.size(); h_Np_vs_Nn->Fill(Nn, Np);
	}
      }

      auto visualizeEvent = [&]() {
        Visualizer3D* vis = new Visualizer3D("EventDisplay", Form("Event %d - 3D Particle Visualization", iEvent));
        vis->SetNeutrinoVertex(nu_reco_vtx_x, nu_reco_vtx_y, nu_reco_vtx_z);
        int n_protons = 0, n_neutrons = 0;
        for(auto g4id : g4ids_for_vis) {
          if(g4id < 0) continue;   // do I need this?
          int idx = map_g4id_index[g4id];
          if(abs(map_g4id_pdg[g4id]) == 2212) n_protons++;   // why abs??
          else if(abs(map_g4id_pdg[g4id]) == 2112) n_neutrons++;

          vis->AddParticle(g4id, map_g4id_pdg[g4id], 
                           part_Px[idx], part_Py[idx], part_Pz[idx], 
                           part_KE[idx], part_isPrimary[idx], map_g4id_mother[g4id]);
        }
        vis->Draw(200);   TString summary = Form("Event Summary: %d protons, %d neutrons", n_protons, n_neutrons);
	vis->GetCanvas()->SetTitle(summary);   vis->GetCanvas()->Update();
      };
      if(create_display) { visualizeEvent(); } 

      totalSelectedBlips+=selectedBlips;
      //set_branch_vars( isOutTree, 0, map_g4id_index, 1, selectedBlips );
      //define_branches(outTree, isOutTree); //outTree->SetAutoSave(0);  // disable autosave

    }//end loop over events

    //handler.printEnergyStatistics();

    merge_plot(h_prBlipsPerE, "Blips", h_prTrksPerE, "Tracks");
    //merge_plot(h_NC_Neutrons, "NC", h_CC_Neutrons, "CC");   // can't run this twice in same macro b/c only 1 canvas (c1)
    //merge_plot(h_cos_pr_all, "Np", h_cos_pr, "1p");   // from primary n-inelastic and contained 
    //merge_plot(h_trueProtonE, "all neutron inelastic protons", h_truePrimPrE, "only with primary neutron mothers");

    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600); 
    h_pr_E_ds->DrawCopy("hist");
    highCutoff_graph->Draw("L same");
    lowCutoff_graph->Draw("L same");

    TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
    gr_cosmic_E_ds->SetMarkerColor(kSpring-5);
    gr_cosmic_E_ds->SetMarkerStyle(20);
    gr_cosmic_E_ds->SetMarkerSize(0.4);
    gr_cosmic_E_ds->Draw("AP");

    gr_elec_E_ds->SetMarkerColor(kViolet-5);
    gr_elec_E_ds->SetMarkerStyle(20);
    gr_elec_E_ds->SetMarkerSize(0.4);
    gr_elec_E_ds->Draw("P same");
 
    gr_signal_E_ds->SetMarkerColor(kOrange+6);
    gr_signal_E_ds->SetMarkerStyle(20);
    gr_signal_E_ds->SetMarkerSize(0.4);
    gr_signal_E_ds->Draw("P same");

    highCutoff_graph->SetLineWidth(5);
    highCutoff_graph->Draw("L same");

    //lowCutoff_graph->SetLineWidth(5);
    //lowCutoff_graph->Draw("L same");

    TLegend *leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(gr_cosmic_E_ds, "Cosmics", "p");
    leg->AddEntry(gr_elec_E_ds, "Electrons", "p");
    leg->AddEntry(gr_signal_E_ds, "Signal Protons", "p");
    //leg->AddEntry(cutoff_graph, "90% of e- blips fall below this line", "L");
    leg->Draw();

    TCanvas *c4 = new TCanvas("c4", "Canvas", 800, 600);
    h_cosAlignment[1]->DrawCopy("hist");

    TCanvas *c5 = new TCanvas("c5", "Canvas", 800, 600);
    h_all_nu_neutrons->DrawCopy("hist");
    //h_Np_vs_Nn->DrawCopy("hist");

    TCanvas *c6 = new TCanvas("c6", "Canvas", 800, 600);
    TH1D* h_pre_Eff_vs_E = (TH1D*)h_preFoundNeutronE->Clone("pre_Eff_vs_E");
    h_pre_Eff_vs_E->SetTitle("(pre)selection Signal Efficiency vs. Energy;True Neutron Energy [MeV];Selection (n, Np) Efficiency");
    h_pre_Eff_vs_E->SetLineColor(kBlack); h_pre_Eff_vs_E->SetMarkerColor(kBlack);
    h_pre_Eff_vs_E->Divide(h_trueNeutronE); h_pre_Eff_vs_E->DrawCopy();
    TH1D* h_pre_drawn = (TH1D*)h_pre_Eff_vs_E->DrawCopy();  // Capture the drawn copy

    TH1D* h_selection_Eff_vs_E = (TH1D*)h_foundNeutronE->Clone("selection_Eff_vs_E");
    h_selection_Eff_vs_E->SetTitle("Selection Signal Efficiency vs. Energy;True Neutron Energy [MeV]");
    h_selection_Eff_vs_E->SetLineColor(kRed); h_selection_Eff_vs_E->SetMarkerColor(kRed);
    h_selection_Eff_vs_E->Divide(h_trueNeutronE); h_selection_Eff_vs_E->DrawCopy("same");
    TH1D* h_sel_drawn = (TH1D*)h_selection_Eff_vs_E->DrawCopy("same");  // Capture the drawn copy

    TLegend *leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg1->AddEntry(h_pre_drawn, "Preselection");
    leg1->AddEntry(h_sel_drawn, "After Cuts");
    leg1->Draw();

    //TH1D* h_cosmicDistCos = (TH1D*)h_cosmicDistToVtx[9]->Clone("cosmicDistCos");
    //h_cosmicDistCos->SetTitle("Cosmic Distance to Vertex;Distance [cm];Cosmics/nBlips");
    //h_cosmicDistCos->Divide(h_cosmicDistToVtx[0]); h_cosmicDistCos->DrawCopy("hist");

    //TCanvas *c8 = new TCanvas("c8", "Canvas", 800, 600);
    //TH1D* h_cosmicDistPr = (TH1D*)h_cosmicDistToVtx[1]->Clone("cosmicDistPr");
    //h_cosmicDistPr->SetTitle("Cosmic Distance to Vertex;Distance [cm];Signal/nBlips");
    //h_cosmicDistPr->Divide(h_cosmicDistToVtx[0]); h_cosmicDistPr->DrawCopy("hist");

    //TH1D* h_efficiency = (TH1D*)h_numerator->Clone("efficiency_name");
    //h_efficiency->Divide(h_numerator, h_denominator, 1.0, 1.0, "B");

    float resVal;
    for( int i=0; i<199; i++ ){
      TH1D* tempProj = h_pr_res_vs_trueE->ProjectionY("tempProj",i,i+1);
      //std::cout<<"Projection integral: "<<tempProj->Integral()<<", count: "<<tempProj->GetEntries()<<"\n";
      resVal = tempProj->GetMean();//tempProj->Integral();//tempProj->GetEntries();
      gr_pr_res_vs_trueE->SetPoint(gr_pr_res_vs_trueE->GetN(),i,resVal);
    }
    TCanvas *c7 = new TCanvas("c7", "Canvas", 800, 600);
    gr_pr_res_vs_trueE->Draw("AP");

    // Get projection of y-axis bins [ybin_low, ybin_high] onto x-axis
    int total_N_intrxns = h_Np_vs_Nn->GetEntries();
    TH1D* proj_1p = h_Np_vs_Nn->ProjectionX("proj1p", 2, 2);
    TH1D* proj_2p = h_Np_vs_Nn->ProjectionX("proj2p", 3, 3);
    double sum_1p = proj_1p->Integral(); double sum_2p = proj_2p->Integral();
    delete proj_1p; delete proj_2p; //sum_1p /= total_N_intrxns; sum_2p /= total_N_intrxns;

    float pr_energy_mean = TMath::Mean(pr_energy.begin(), pr_energy.end());
    float avg_cos_allBlips = TMath::Mean(cos_allBlips.begin(), cos_allBlips.end());
    float totalOtherInterxns = h_blipE[7]->GetEntries()+h_blipE[8]->GetEntries();

    std::cout<<"\nNumber of events in sample: "<<numEvents<<
    " events\nNumber of neutrinos in sample: "<<total_nus<<
    " neutrinos\nNumber of neutrons in sample: "<<total_Ns<<
    " neutrons\nNumber of neutrino-induced neutrons (primary) in sample: "<<total_nu_Ns<<
    " neutrons\nPercent of neutrinos that produce neutrons (primary): "<<total_nu_Ns/total_nus*100<<
    "%\nPercent of neutrons that produce protons: "<<total_N_prs/total_Ns*100<<
    "%\nPercent of primary neutrons that produce protons: "<<num_n_p/total_nu_Ns*100<<

    "%\n\nNumber of neutron inelastic protons in sample: "<<total_N_prs<<
    " protons\nNumber of primary nu-ind (n, p) in sample: "<<num_n_p<<
    " protons\nMean proton blip energy: "<<pr_energy_mean<<
    " MeV\nTotal number of (primary) N-ind proton blips: "<<signal_blips<<
    " blips\nNumber of (primary) N-ind SINGLE proton blips: "<<single_pr<<
    " blips\nFraction of neutron interactions that are (1p, Nn) = "
	<<sum_1p<<"/"<<total_N_intrxns<<" = "<<sum_1p/total_N_intrxns<<
    "\nFraction of neutron interactions that are (2p, Nn) = "
        <<sum_2p<<"/"<<total_N_intrxns<<" = "<<sum_2p/total_N_intrxns<<
    "\nExit detector without inelastic scatter: "<<exitWoutScat/total_nu_Ns*100<<
    //"%\nExit detector after inelastic scatter: "<<hasScat/total_nu_Ns*100<<
    //"%\nCaptured by argon: "<<getsCaptured/total_nu_Ns*100<<
    //"%\nOther: "<<Nother/total_nu_Ns*100<<

    "%\n\nTotal gammas in file: "<<total_gammas<<
    " gammas\nTotal electrons in file: "<<total_electrons<<
    " electrons\nTotal gamma blips: "<<total_gamma_blips<<
    " gammas\nTotal N-ind gamma blips: "<<total_N_gammas<<
    " gammas\nFraction of 'other' gammas that are eBrems: "<<
	h_blipE[8]->GetEntries()/totalOtherInterxns<<

    " gammas\nFraction of 'other' interactions that aren't gammas: "<<otherNotGamma/totalOtherInterxns<<
    "\nFraction of 'other' interactions that ARE gammas: "<<otherIsGamma/totalOtherInterxns<<

    " gammas\n\nTotal blips: "<<totalBlips<<
    " blips\nTotal blips that passed all cuts: "<<totalSelectedBlips<<
    " blips\nTotal truth matched blips passed cuts: "<<trueBlips<<
    " blips\nTotal primary N-induced blips: "<<total_N_blips<<
    " blips\nAverage cosine ALL blips: "<<avg_cos_allBlips<<

    "\n\nEfficiency: "<<signal_blips/total_nu_Ns*100<<
    "%\nPurity: "<<signal_blips/totalSelectedBlips*100<<"%\n\n";

    Float_t a, b;
    TTree otree("evd_tree", "TTree with a single float branch for total evd"); // two branches
    TBranch* a_branch = otree.Branch("totalEvents", &a, "totalEvents/F");   //total_blips
    //TBranch* b_branch = otree.Branch("signal_blips", &b, "signal_blips/F"); b = signal_blips;
    a = numEvents; otree.Fill();

    //===========================================
    // Write the outputs and close the files

    // output tree file
    //remainFile->cd();
    //outTree->Write();
    //remainFile->Close();

    // input file
    //resFile->Close();
    //file->Close();
  
    // histogram file
    outfile->cd();
    outfile->Write();
    outfile->Close();
    
  }
