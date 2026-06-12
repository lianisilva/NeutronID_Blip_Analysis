  #include "../variables.h"
 //#include "~/root/pelee/out"

  #include <iostream>
  #include <cmath>
  #include <vector>
  #include "TROOT.h"
  #include "TMath.h"
  #include "TVectorT.h"
  #include <algorithm>
 
  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;
  
  // Input file name
  //std::string fFileName = "../blipana_genie_fullosc_overlay_run3_20250304.root";
  std::string fFileName = "../blipana_MCC9_bnb_nu_overlay_mcc9_run3_20250904.root";
  //std::string fFileName = "MCC9.10_Test_Samples_v10_04_07_04_Run4b_pandora_standalone_reco2_BNB_nu_overlay_reco2_reco2_blipana_8016evts.root";
  //std::string fFileName = "MCC9.10_Test_Samples_v10_04_07_04_Run4b_pandora_standalone_reco2_BNB_nu_overlay_reco2_reco2_blipana_test.root";
  //std::string fFileName = "MCC9.10_Test_Samples_v10_04_07_04_Run4b_pandora_standalone_reco2_BNB_nu_overlay_reco2_reco2_blipana.root";
  //std::string fFileName = "SURPRISE_Test_Samples_v10_04_07_05_Run4b_pandora_standalone_reco2_BNB_nu_overlay_reco2_blipana.root";

  // only deal with contained neutrons in calculation and plotting
  TH1D*   h_ContainedNeutrons;      // neutrino-induced neutrons (NINs) per neutrino
  TH1D*   h_NIblips;		    // neutrino-induced blips
  TH1D*   h_NC_Neutrons;	    // Neutral current NINs per neutrino
  TH1D*   h_CC_Neutrons;	    // Charged current NINs per neutrino
  TH1D*   h_BlipsPerEvent;
  TH1D*   h_blipE;
  TH1D*   h_nblips_vs_vtxDist;
  TH1D*   h_blipXdistToVtx;
  TH1D*   h_blipYdistToVtx;
  TH1D*   h_blipZdistToVtx;
  TH1D*   h_maxWires;
  TH2D*   h_blipZY;
  TH1D*   h_numTrks;
  TH1D*   h_trackLength;


  void validation(){
 
    // Read in the TFile and grab the TTree
    TFile* file = new TFile(fFileName.c_str(),"READ");
    TTree* tree = (TTree*)file->Get("blipana/anatree");
    std::cout<<"Reading data file "<<fFileName<<"\n";

    // Set branch address for all variables you plan to use
    // (if a branch address isn't set, the variables will 
    // not be filled when each event is read). For a description 
    // of each variable, see the comments in 'variables.h'
    tree->SetBranchAddress("event",           &event);
    tree->SetBranchAddress("run",             &run);
    tree->SetBranchAddress("subrun",          &subrun);
    tree->SetBranchAddress("nparticles",      &nparticles);
    tree->SetBranchAddress("part_process",    &part_process);
    tree->SetBranchAddress("part_isPrimary",  &part_isPrimary);
    tree->SetBranchAddress("part_isContained",&part_isContained);
    tree->SetBranchAddress("part_g4id",       &part_g4id);
    tree->SetBranchAddress("part_pdg",        &part_pdg);
    tree->SetBranchAddress("part_mother",     &part_mother);
    tree->SetBranchAddress("part_KE",         &part_KE);
    tree->SetBranchAddress("part_pathlen",    &part_pathlen);
    tree->SetBranchAddress("part_depEnergy",  &part_depEnergy);
    tree->SetBranchAddress("part_nDaughters", &part_nDaughters);
    tree->SetBranchAddress("nblips",          &nblips);
    tree->SetBranchAddress("blip_x",          &blip_x);
    tree->SetBranchAddress("blip_y",          &blip_y);
    tree->SetBranchAddress("blip_z",          &blip_z);
    tree->SetBranchAddress("blip_size",       &blip_size);
    tree->SetBranchAddress("blip_energy",     &blip_energy);
    tree->SetBranchAddress("blip_edepid",     &blip_edepid);
    tree->SetBranchAddress("blip_g4id",       &blip_g4id);
    tree->SetBranchAddress("blip_dw",         &blip_dw);
    tree->SetBranchAddress("blip_dw",         &blip_dw);
    tree->SetBranchAddress("ntrks",           &ntrks);
    tree->SetBranchAddress("trk_length",      &trk_length);
    tree->SetBranchAddress("mctruth_nu_ccnc", &mctruth_nu_ccnc);
    tree->SetBranchAddress("nu_reco_vtx_x",   &nu_reco_vtx_x);
    tree->SetBranchAddress("nu_reco_vtx_y",   &nu_reco_vtx_y);
    tree->SetBranchAddress("nu_reco_vtx_z",   &nu_reco_vtx_z);


    // make output file to store plots
    TFile* outfile = new TFile("out.root", "recreate");

    // Make histograms declared above
    h_ContainedNeutrons     = new TH1D("ContainedNeutrons","Number of neutrino-induced neutrons; # neutrons",15,0,15);
    h_blipE		    = new TH1D("blipE","Blip Energy; E (MeV)",50,0,25);
    h_BlipsPerEvent	    = new TH1D("BlipsPerEvent","Number of Blips;# blips;events",250,0,500);
    h_nblips_vs_vtxDist	    = new TH1D("nblips_vs_vtxDist",";Distance to reco vtx [cm]; Nblips / r^2",124,0,248);
    h_blipXdistToVtx        = new TH1D("blipXdistToVtx","Blip X Distance to Reco Vertex;x dist",125,-250,250);
    h_blipYdistToVtx        = new TH1D("blipYdistToVtx","Blip Y Distance to Reco Vertex;y dist",125,-250,250);
    h_blipZdistToVtx        = new TH1D("blipZdistToVtx","Blip Z Distance to Reco Vertex;z dist",125,-250,250);
    h_maxWires		    = new TH1D("maxWires","Max Wires per Blip on one Plane;# Wires",18,0,18);
    h_blipZY                = new TH2D("blipZY","Blip Z vs Y;z;y",224,-40,1080,50,-150,150);
    h_numTrks               = new TH1D("numTrks","Number of Tracks;# tracks",200,0,200);
    h_trackLength	    = new TH1D("trackLengths","Track Lengths,length [cm]",201,-2,400);

    int numEvents = tree->GetEntries();

    int   total_nus=0;
    int   total_nu_Ns=0;
    int   total_contained=0;
    int   totalBlipE=0;
    int   totalBlips=0;
    int   totalSampleBlips=0; 

    float minX=0, minY=0, minZ=0;
    float maxX=0, maxY=0, maxZ=0;

    int selectedEvents=0;   // those with nblips > 0 and containing a reco nu vtx
    int nonSelected=0;

    std::vector<float> distToVtx;
    std::vector<float> energy;
    std::vector<float> X;
    std::vector<float> Y;
    std::vector<float> Z;

    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      tree->GetEntry(iEvent);

      total_nus++;   // one neutrino per event
      int   close_nu_blips = 0;

      h_numTrks->Fill(ntrks);
      for( int i=0; i<ntrks; i++){ h_trackLength->Fill(trk_length[i]); }

      //bool hasRecoVtx=false;
      //if( fabs(reco_nu_vtx_x)<1000 && fabs(reco_nu_vtx_y)<1000 && reco_nu_vtx_z>-100 ) hasRecoVtx = true;
      //if( hasRecoVtx == false ) { nonSelected++; continue; }
      selectedEvents++;


      TVector3 reco_nu_vtx( nu_reco_vtx_x, nu_reco_vtx_y, nu_reco_vtx_z );
      int goodBlips = 0;

      for(int i=0; i<nblips; i++){
        totalBlips++;

        float blipX=blip_x[i], blipY=blip_y[i], blipZ=blip_z[i];
	X.push_back(blipX); Y.push_back(blipY); Z.push_back(blipZ);

        TVector3 blip_loc( blipX, blipY, blipZ );
        float dist = (blip_loc-reco_nu_vtx).Mag();
        //if( (unsigned)blip_x[i]<=250 && (unsigned)(blip_y[i]+200)<=400 && (unsigned)blip_z[i]<=1000 )   // (number-lower) <= (upper-lower)

        //if( dist < 250 ) {

          if( minX > blipX ) minX=blipX; if( minY > blipY ) minY=blipY; if( minZ > blipZ ) minZ=blipZ;
          if( maxX < blipX ) maxX=blipX; if( maxY < blipY ) maxY=blipY; if( maxZ < blipZ ) maxZ=blipZ;

          if( blip_energy[i] >= 0 ) energy.push_back(blip_energy[i]);
          distToVtx.push_back(dist);

	  goodBlips++;
	  close_nu_blips++;
	  totalSampleBlips++;

	  h_blipE->Fill(blip_energy[i]);
	  h_nblips_vs_vtxDist->Fill( dist, 1./pow(dist,2) );   // normalize by 1/r^2, make spherical dist. flat
       	  h_blipXdistToVtx->Fill(blipX - nu_reco_vtx_x);
          h_blipYdistToVtx->Fill(blipY - nu_reco_vtx_y);
          h_blipZdistToVtx->Fill(blipZ - nu_reco_vtx_z);
          h_blipZY->Fill(blipZ,blipY);

	  int maxWires = std::round(blip_dw[i]/0.3);
	  h_maxWires->Fill(maxWires);   // in cm so convert to Nwires by /0.3 cm

        //} 

        totalBlipE+=blip_energy[i];

      }//end blip loop

      if( goodBlips<=0 ) std::cout<<"\nNO GOOD BLIPS!\n";
      //if( goodBlips<=0 ) { nonSelected++; --selectedEvents; }
      h_BlipsPerEvent->Fill(close_nu_blips);

    }//end loop over events

    h_BlipsPerEvent->Scale(1.0 / numEvents);

    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    h_numTrks->DrawCopy("hist");

    TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
    h_trackLength->DrawCopy("hist");

    //merge_plot(h_nblips_vs_vtxDist, "Run 3", h_v_nblips_vs_vtxDist, "Run 4b");   // classic Blip Reco vs. nTuple Blips

    float dist_sum = std::reduce(distToVtx.begin(), distToVtx.end());
    float dist_size = static_cast<int>(distToVtx.size());
    float dist_mean = dist_sum/dist_size;

    float energy_sum = std::reduce(energy.begin(), energy.end());
    float energy_size = static_cast<int>(energy.size());
    float energy_mean = energy_sum/energy_size;
    //float energy_mean = TMath::Mean(energy.begin(), energy.end());
    //float energy_mean = totalBlipE/totalBlips;

    float X_sum = std::reduce(X.begin(), X.end()), X_size = static_cast<int>(X.size()), X_mean = X_sum/X_size;
    float Y_sum = std::reduce(Y.begin(), Y.end()), Y_size = static_cast<int>(Y.size()), Y_mean = Y_sum/Y_size;
    float Z_sum = std::reduce(Z.begin(), Z.end()), Z_size = static_cast<int>(Z.size()), Z_mean = Z_sum/Z_size;

    printf("\nBlip multiplicity: %8.3f +/- %.3f\n",
	   double(totalSampleBlips)/double(selectedEvents), h_BlipsPerEvent->GetRMS(11)   );


    std::cout<<"\nNumber of sample events: "<<selectedEvents<<
    " events\nNumber of sample blips (<250cm from nu vtx): "<<totalSampleBlips<<
    " events\nNumber of neutrinos in sample: "<<total_nus<<
    " neutrinos\nNumber of blips in sample: "<<totalBlips<<

    " blips\n\nMaximum blip (x, y, z): ("<<maxX<<", "<<maxY<<", "<<maxZ<<")"<<
    "\nMinimum blip (x, y, z): ("<<minX<<", "<<minY<<", "<<minZ<<")\n"<<
    "Mean blip (x, y, z): ("<<X_mean<<", "<<Y_mean<<", "<<Z_mean<<")\n\n"<<

    "Mean dist to vtx: "<<dist_mean<<"\n"<<
    "Mean energy: "<<energy_mean<<"\n\n"<<

    "Events without blips: "<<nonSelected<<
    " events\nPercent of sample events w/out reco nu vtx: "<<nonSelected*100/(selectedEvents+nonSelected)<<"%\n\n";

    Float_t a;
    TTree otree("evd_tree", "TTree with a single float branch for total evd");
    TBranch* branch = otree.Branch("total_event_displays", &a, "total_event_displays/F");
    a = selectedEvents;   //containedEvents;   //numEvents;
    otree.Fill();

    //===========================================
    // Write the output and close the file
    outfile->Write();
    outfile->Close();

  }
