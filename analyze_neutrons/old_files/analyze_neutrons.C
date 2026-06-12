  ////////////////////////////////////////////////////////
  //  Blip Analysis ROOT Macro
  //
  //  This ROOT/C++ analysis code reads in a TTree produced
  //  by the MicroBooNE 'BlipAna' module and makes some plots
  //  of different blip-related metrics.
  //
  //  This macro must be in the same directory as the header file
  //  "variables.h" as well as a BlipAna ROOT file. To run, type
  //  this command in the terminal:
  //  
  //  > root -l analyze_neutrons.C
  //
  //  W. Foreman
  //  Jan 2025
  ////////////////////////////////////////////////////////
  
 
  // We first include this header file, which contains all the 
  // variables stored in the ROOT TTree we will be reading from.
  #include "variables.h"
  #include "merge_plot.C"

  #include "TMath.h"

  #include <iostream>
  #include <cmath>
  #include <vector>
  #include "TROOT.h"
  #include "TVectorT.h"
  #include <algorithm>
 
  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;
  
  // Input file name
  //std::string fFileName = "blipana_protons_0to200MeV_overlay_20250313.root";
  //std::string fFileName = "blipana_MCC9_bnb_nu_overlay_mcc9_run3_20250904.root";
  std::string fFileName = "blipana_protons_0to200MeV_overlay_20250826.root";

  // Tell ROOT about the histograms we'll be making later
  TH1D*   h_TrueFSGammaE;           // individual neutron capture final-state gamma energies
  TH1D*   h_TrueFSGammaESum;        // summed energies of capture gammas
  TH1D*   h_TrueProtonE;
  TH1D*   h_ProtonBlipE;
  TH1D*   h_TrueProtonLength;	    // path length of secondary protons
  TH1D*   h_TrackLength;	    // reconstructed proton track length
  TH2D*   h_TrueNeutronEvsProtonE;  // neutron kinetic energy vs resultant proton KE
  TH2D*   h_TrueProtonEvsLength;    // proton energy vs pathlength
  TH2D*   h_TrueLengthvsTrkLength;  // proton length vs track length
  TH2D*   h_ProtonEvsTrkLength;     // proton energy vs track length
  TH2D*   h_TrueLengthvsBlipE;      // ! change x & y !
  TH2D*   h_BlipEvsTrkLength;       // ^^^ 
  TH1D*   h_NumTrksPerE;
  TH1D*   h_NumBlipsPerE;
  TH1D*   h_TrkEfficiency;
  TH1D*   h_BlipEfficiency;
  TH1D*   h_BlipOrTrackEff;
  TH1D*   h_BlipXorTrack;
  TH2D*   h_TrueLengthvsBlipSize;
  TH2D*   h_BlipEvsAvgWires;	    // blip (true proton) energy vs avg number of wires (over all planes) in the cluster that makes up the reco blip
  TH2D*   h_BlipEvsAvgWires_0;
  TH2D*   h_BlipEvsAvgWires_1;
  TH2D*   h_BlipEvsAvgWires_2;
  TH2D*   h_BlipEvsNumWires;
  TH2D*   h_ProtonEvsTimespan;
  TH1D*   h_2and2wires;
  TH1D*   h_2and3wires;
  TH1D*   h_3and3wires;
  TH1D*   h_3and4wires;

  // ********************************************************
  // This is the main function that is called when the call is
  // made to ROOT to process the macro. This must match the name
  // of the file ('analyze_neutrons.C')
  // *********************************************************
  void analyze_neutrons(){
  
    // Read in the TFile and grab the TTree
    TFile* file = new TFile(fFileName.c_str(),"READ");
    TTree* tree = (TTree*)file->Get("blipana/anatree");
    std::cout<<"Reading data file "<<fFileName<<"\n";
    
    // Set branch address for all variables you plan to use
    // (if a branch address isn't set, the variables will 
    // not be filled when each event is read). For a description 
    // of each variable, see the comments in 'variables.h'
    tree->SetBranchAddress("event",           &event);
    tree->SetBranchAddress("run",	      &run);
    tree->SetBranchAddress("subrun",          &subrun);
    tree->SetBranchAddress("nparticles",      &nparticles);
    tree->SetBranchAddress("part_isPrimary",  &part_isPrimary);
    tree->SetBranchAddress("part_isContained",&part_isContained);
    tree->SetBranchAddress("part_g4id",       &part_g4id);
    tree->SetBranchAddress("part_pdg",        &part_pdg);
    tree->SetBranchAddress("part_mother",     &part_mother);
    tree->SetBranchAddress("part_KE",         &part_KE);
    tree->SetBranchAddress("part_pathlen",    &part_pathlen);
    tree->SetBranchAddress("part_depEnergy",  &part_depEnergy);
    tree->SetBranchAddress("part_process",    &part_process);
    tree->SetBranchAddress("nedeps",          &nedeps);
    tree->SetBranchAddress("blip_g4id",       &blip_g4id);
    tree->SetBranchAddress("edep_pdg",        &edep_pdg);
    tree->SetBranchAddress("edep_blipid",     &edep_blipid);
    tree->SetBranchAddress("edep_energy",     &edep_energy);
    tree->SetBranchAddress("edep_x",          &edep_x);
    tree->SetBranchAddress("edep_y",          &edep_y);
    tree->SetBranchAddress("edep_z",          &edep_z);
    tree->SetBranchAddress("ntrks",           &ntrks);
    tree->SetBranchAddress("trk_g4id",        &trk_g4id);
    tree->SetBranchAddress("trk_length",      &trk_length);
    tree->SetBranchAddress("trk_startx",      &trk_startx);
    tree->SetBranchAddress("trk_starty",      &trk_starty);
    tree->SetBranchAddress("trk_startz",      &trk_startz);
    tree->SetBranchAddress("trk_endx",        &trk_endx);
    tree->SetBranchAddress("trk_endy",        &trk_endy);
    tree->SetBranchAddress("trk_endz",        &trk_endz);
    tree->SetBranchAddress("nblips",          &nblips);                       
    tree->SetBranchAddress("blip_x",          &blip_x);                        
    tree->SetBranchAddress("blip_y",          &blip_y);                        
    tree->SetBranchAddress("blip_z",          &blip_z);
    tree->SetBranchAddress("blip_size",	      &blip_size);                       
    tree->SetBranchAddress("blip_dw",	      &blip_dw); 
    tree->SetBranchAddress("blip_energy",     &blip_energy);
    tree->SetBranchAddress("blip_edepid",     &blip_edepid);              
    tree->SetBranchAddress("blip_nplanes",    &blip_nplanes);
    tree->SetBranchAddress("blip_pl0_clustid",&blip_pl0_clustid);
    tree->SetBranchAddress("blip_pl1_clustid",&blip_pl1_clustid); 
    tree->SetBranchAddress("blip_pl2_clustid",&blip_pl2_clustid);
    tree->SetBranchAddress("clust_nwires",    &clust_nwires);
    tree->SetBranchAddress("clust_timespan",  &clust_timespan);
    tree->SetBranchAddress("clust_bydeadwire",&clust_bydeadwire);
    tree->SetBranchAddress("clust_ismatch",   &clust_ismatch);
    tree->SetBranchAddress("nhits",           &nhits);
    tree->SetBranchAddress("hit_blipid",      &hit_blipid);
    tree->SetBranchAddress("hit_plane",       &hit_plane);
    tree->SetBranchAddress("hit_wire",        &hit_wire); 
    tree->SetBranchAddress("hit_amp",         &hit_amp);
    tree->SetBranchAddress("hit_peakT",	      &hit_peakT);
    tree->SetBranchAddress("hit_driftT",      &hit_driftT);
    tree->SetBranchAddress("hit_rms",	      &hit_rms); 

    // make output file to store plots
    TFile* outfile = new TFile("out.root", "recreate");
  
    // Make histograms declared above
    h_TrueFSGammaE          = new TH1D("TrueFSGammaE","Neutron capture final-state gammas;true energy [MeV]",100,0,30);
    h_TrueFSGammaESum       = new TH1D("TrueFSGammaESum","Neutron capture final-state gammas;summed true energy [MeV]",200,0,15);
    h_TrueProtonE           = new TH1D("TrueProtonE","True Proton Energy [MeV]",100,0,200);
    h_ProtonBlipE	    = new TH1D("ProtonBlipE","Blip energy [MeV]",100,0,100);
    h_TrueProtonLength	    = new TH1D("TrueProtonLength","True Proton Path Length [cm]",100,0,7);
    h_TrackLength	    = new TH1D("TrackLength","Reconstructed Proton Track Length;Track length [cm]",100,0,14);
    h_TrueNeutronEvsProtonE = new TH2D("TrueNeutronEvsProtonE","Proton Energy vs Neutron Energy;Neutron kinetic energy [MeV];FS proton kinetic energy [MeV]",100,0,205,100,0,200);
    h_TrueProtonEvsLength   = new TH2D("TrueProtonEvsLength","True Proton Pathlength vs Energy;Proton energy [MeV];True pathlength [cm]",100,0,200,100,0,14);
    h_TrueLengthvsTrkLength = new TH2D("TrueLengthvsTrkLength","True Proton Path Length vs Reconstructed Track Length;True pathlength [cm];Track length [cm]",100,0,7,100,0,14);
    h_ProtonEvsTrkLength    = new TH2D("ProtonEvsTrkLength","True Proton Energy vs Reco Track Length;Track length [cm];Proton energy [MeV]",100,0,65,100,0,200);
    h_TrueLengthvsBlipE	    = new TH2D("TrueLengthvsBlipE","True Proton Path Length vs Blip Energy;Blip energy [MeV];True pathlength [cm]",100,0,50,100,0,14);
    h_BlipEvsTrkLength      = new TH2D("BlipEvsTrkLength","Reco Track Length vs Blip Energy;Blip energy [MeV];Track length [cm]",100,0,185,100,0,18);
    h_NumTrksPerE	    = new TH1D("NumTrksPerE","Number of Protons Made Into Tracks vs Energy;Proton energy [MeV];Number of Tracks",100,0,200);
    h_NumBlipsPerE	    = new TH1D("NumBlipsPerE","Number of Protons Made Into Blips vs Energy;Proton energy [MeV];Number of Blips",100,0,200);
    h_TrkEfficiency	    = new TH1D("TrkEfficiency","Reconstruction Efficiency vs Energy;True proton energy [MeV];Track reco efficiency",100,0,200);
    h_BlipEfficiency	    = new TH1D("BlipEfficiency","Reconstruction Efficiency vs Energy;True proton energy [MeV];Blip reco efficiency",100,0,200);
    h_BlipOrTrackEff        = new TH1D("BlipOrTrackEff","Proton Reconstruction as Track or Blip Efficiency vs Energy;True proton energy [MeV];Efficiency",100,0,200);
    h_BlipXorTrack	    = new TH1D("BlipXorTrack","Proton NOT Reconstructed as Track or Blip Probability vs Energy;True proton energy [MeV];Probability",100,0,200);
    h_TrueLengthvsBlipSize  = new TH2D("TrueLengthvsBlipSize","True Proton Path Length vs Blip Size;Blip size [cm];True pathlength [cm]",100,0,4,100,0,4);
    h_BlipEvsAvgWires	    = new TH2D("BlipEvsAvgWires","Blip Energy vs Avg Wires per Blip (Cluster);# wires;Proton blip energy [MeV]",12,0,12,50,0,100);
    h_BlipEvsAvgWires_0     = new TH2D("BlipEvsAvgWires_0","Blip Energy vs Avg Wires per Blip on Plane 0;# wires;Proton blip energy [MeV]",14,0,14,35,0,70);
    h_BlipEvsAvgWires_1     = new TH2D("BlipEvsAvgWires_1","Blip Energy vs Avg Wires per Blip on Plane 1;# wires;Proton blip energy [MeV]",14,0,14,35,0,70);
    h_BlipEvsAvgWires_2     = new TH2D("BlipEvsAvgWires_2","Blip Energy vs Avg Wires per Blip on Plane 2;# wires;Proton blip energy [MeV]",14,0,14,35,0,70);
    h_BlipEvsNumWires       = new TH2D("BlipEvsNumWires","Blip Energy vs Number of Wires per Blip (Cluster);# wires;Proton blip energy [MeV]",16,0,16,100,0,100);
    h_ProtonEvsTimespan     = new TH2D("ProtonEvsTimespan","True Proton Energy vs Max Cluster Timespan;Timespan [ticks];Proton energy [MeV]",50,0,50,25,0,200);
    h_2and2wires	    = new TH1D("2and2wires","Blips Spanning Multiple Wires on 2 Planes;Proton energy [MeV];Number of blips",100,0,100);   // # Blips Spanning >= 2 Wires on 2 Planes
    h_2and3wires            = new TH1D("2and3wires","# Blips >= 2 and >= 3 Wires on 2 Respective Planes;Proton energy [MeV];Number of blips",100,0,100);
    h_3and3wires            = new TH1D("3and3wires","# Blips Spanning >= 3 Wires on 2 Planes;Proton energy [MeV];Number of blips",100,0,100);
    h_3and4wires            = new TH1D("3and4wires","# Blips >= 3 and >= 4 Wires on 2 Respective Planes;Proton energy [MeV];Number of blips",100,0,100);


    // set any 2D histograms to "colz" option (bin values 
    // displayed with rainbow color scale)
    h_TrueNeutronEvsProtonE->SetOption("colz");
    h_TrueProtonEvsLength->SetOption("colz");
    h_TrueLengthvsTrkLength->SetOption("colz");
    h_ProtonEvsTrkLength->SetOption("colz");
    h_BlipEvsAvgWires->SetOption("colz");
    h_ProtonEvsTimespan->SetOption("colz");

    // Count the entries in this TTree file. Each entry corresponds
    // to one readout of the detector, referred to as an "event"
    int numEvents = tree->GetEntries();
    //std::cout<<"There are "<<numEvents<<" events in this file.\n";
    
    int   total_second_prs = 0;
    int   total_pr_trks = 0;
    int   total_pr_notrks = 0;
    int   total_pr_blips = 0;
    int   total_neutrons = 0;
    float min_trk_E = 10;
    float total_trk_E = 0;
    float total_notrk_E = 0;
    float total_blip_E = 0;
    int   blips_by_death = 0;
    int   max_nwires = 0;

    //====================================================
    // Loop through the events. For each one, we will call
    // the "GetEntry" command, which tells ROOT to fill all
    // the variables with set branch addresses (don't worry, 
    // this will all make more sense a bit later on...)
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      tree->GetEntry(iEvent);
      
      //std::cout<<"=======================================\n";
      //std::cout<<"Reading event "<<iEvent<<", which has "<<nblips<<" blips\n";

      // keep track of info about the neutrons
      std::vector<int> neutron_indices; // was -9
      //int   neutron_g4id  = -9;
      std::vector<int> neutron_KEs;
      float neutron_FSgamma_sumE = 0;
      int   neutron_num = 0;

      // keep track of info about the secondary proton (if produced)
      int   proton_num    = 0;
      std::vector<int> proton_indices;


      //==========================================================
      // Loop over the true MCParticles in this event
      //std::cout<<"Event has "<<nparticles<<" MC particles\n";      

      for(int i=0; i<nparticles; i++){

        // -------------------------------------------------------------
        // Printout this particle information to the screen 
	if(0){
          printf("  %5i  G4ID: %-8i PDG: %-10i  dL=%7.1f cm,  KE0=%9.3f MeV,  Edep=%9.3f MeV,  moth=%5i, %12s, NDaughters=%i\n",
            i, part_g4id[i], part_pdg[i], part_pathlen[i], part_KE[i], part_depEnergy[i],
	    part_mother[i], part_process->at(i).c_str(), part_nDaughters[i]
          );
	}

	// -------------------------------------------------------------
	// sum over ALL neutrons in the event and save their information
	
       if ( part_pdg[i]==2112 && part_isPrimary[i] ){
	  neutron_num++;
	  total_neutrons++;
	  neutron_indices.push_back(i);
          //neutron_g4ids = part_trackID[i];
	  neutron_KEs.push_back(part_KE[i]);
	}
      }//end loop over all MCParticles

      for(int i=0; i<nparticles; i++){

        // -------------------------------------------------------------
        // IF the neutron was already found...
        //if( neutron_index >= 0 )
        
          // identify any gammas (PDG=22) that come from that neutron capture
          if( part_pdg[i] == 22 && part_process->at(i).compare("nCapture") == 0 ) {   //  && part_mother[i] == neutron_g4id
            h_TrueFSGammaE->Fill(part_KE[i]);
            neutron_FSgamma_sumE += part_KE[i];
          } 

          // identify any protons produced by neutron inelastic collisions that do not leave the AV 
	  if( part_pdg[i] == 2212 && part_isContained[i] ) { //&& part_process->at(i).compare("neutronInelastic") == 0) { // part_isPrimary[i] ) { 
            proton_num++; total_second_prs++;
            proton_indices.push_back(i);
          }
	
        } //end loop over all MCParticles

      h_TrueFSGammaESum ->Fill(neutron_FSgamma_sumE);

      // set blipID
      std::map<int,std::vector<int>> map_blipID_hitID;
      for(int ih=0; ih<nhits; ih++){
	if (hit_blipid[ih] >= 0) map_blipID_hitID[hit_blipid[ih]].push_back(ih);
       }

      // check hits per blip
      if(0){
	for(int i=0; i<nblips; i++) {
	  if( blip_edepid[i] > 0 ) {
	    //std::cout<<"Hits associated with blip "<<i<<":\n";
	    for( auto hi : map_blipID_hitID[i] ) { std::cout<<"  hit "<<hi<<"\n"; }
	  }
        }
      }

      // loop through all the protons in the event
      for(int i=0; i<proton_num; i++) {
        //std::vector<int> bclust_nwires;   // store the number of wires for each blip (cluster) in an event
	//std::cout<<"\nproton loop #"<<i<<", proton indice: "<<proton_indices[i];
	float max_timespan = 0;
        bool trk = false;
        bool blip = false;
	float proton_KE   = part_KE[proton_indices[i]];
        float proton_len  = part_pathlen[proton_indices[i]];

        h_TrueProtonE->Fill(proton_KE);
	h_TrueProtonLength->Fill(proton_len);
        //h_TrueNeutronEvsProtonE->Fill(neutron_KEs[i],proton_KE);   // can't do this, causes segmentation violation, use g4_id or sumn?
        h_TrueProtonEvsLength->Fill(proton_KE,proton_len);

	int pr_blips = 0;
	bool bydeadwire = false;
        // Check if the proton got reconstructed as a blip
	for(int j=0; j<nblips; j++){
          if( blip_g4id[j] == part_g4id[proton_indices[i]] ) { 
	  blip=true;

	    if( clust_bydeadwire[blip_pl0_clustid[j]] == 1 || clust_bydeadwire[blip_pl1_clustid[j]] == 1 || clust_bydeadwire[blip_pl2_clustid[j]] == 1 ) {
	      bydeadwire = true;
	      blips_by_death++;
	      //continue;
	    }

	    float blip_E = blip_energy[j];
            total_blip_E += proton_KE;
	    total_pr_blips++;
	    pr_blips++;

	    // number of wires in the cluster that makes up the blip
	    int nwires0 = clust_nwires[blip_pl0_clustid[j]];
	    int nwires1 = clust_nwires[blip_pl1_clustid[j]];
	    int nwires2 = clust_nwires[blip_pl2_clustid[j]];
	    // find max wires
	    if( nwires0 > max_nwires ) { max_nwires = nwires0; }
	    if( nwires1 > max_nwires ) { max_nwires = nwires1; }
            if( nwires2 > max_nwires ) { max_nwires = nwires2; }

	    //float timespan = clust_timespan[blip_pl2_clustid[j]];
	    std::vector<int> list_nwires;  // also list blip indices to get E later
	    if( nwires0 > 0 ) { list_nwires.push_back(nwires0); }
	    if( nwires1 > 0 ) { list_nwires.push_back(nwires1); }
	    if( nwires2 > 0 ) { list_nwires.push_back(nwires2); }
	    //if( timespan > max_timespan ) { max_timespan = timespan; }

	    if( (nwires0>1) + (nwires1>1) + (nwires2>1) >= 2 ) {
	      h_2and2wires->Fill(proton_KE);
	      if( nwires0 > 1 ? (nwires1>2||nwires2>2) : nwires1 > 1 ? (nwires2>2||nwires0>2) : nwires2>1 && (nwires2>2||nwires0>2))
		h_2and3wires->Fill(proton_KE);
		if( (nwires0>2) + (nwires1>2) + (nwires2>2) >= 2 )
		  h_3and3wires->Fill(proton_KE);
		  if( nwires0 > 2 ? (nwires1>3||nwires2>3) : nwires1 > 2 ? (nwires2>3||nwires0>3) : nwires2>2 && (nwires2>3||nwires0>3))
		    h_3and4wires->Fill(proton_KE);
	    }

	    //h_BlipEvsNumWires->Fill(num_wires, proton_KE);
	    float avgWires = (nwires0+nwires1+nwires2)/3;
	    h_BlipEvsAvgWires->Fill(avgWires, proton_KE);
	    h_NumBlipsPerE->Fill(proton_KE);
            h_ProtonBlipE->Fill(proton_KE);
	    h_TrueLengthvsBlipE->Fill(proton_KE, proton_len);
	    h_TrueLengthvsBlipSize->Fill(blip_size[j], proton_len);
 
          }
        }//end blip loop


	int num_trks = 0;
        // Check if the proton got reconstructed in as a
        // 3D track (produced by Pandora, arXiv:1708.03135)
        for(int j=0; j<ntrks; j++) {
          if( trk_g4id[j] == part_g4id[proton_indices[i]] ) {
	    trk=true;
	    num_trks++;
	    total_pr_trks++;
	    total_trk_E += proton_KE;

	    h_NumTrksPerE->Fill(proton_KE); 
	    h_TrackLength->Fill(trk_length[j]);
	    h_TrueLengthvsTrkLength->Fill(proton_len, trk_length[j]);
	    h_ProtonEvsTrkLength->Fill(trk_length[j], proton_KE);

	    if( proton_KE < min_trk_E ) {min_trk_E = proton_KE;}

	    for(int k=0; k<nblips; k++){
	      if( blip_g4id[k] == trk_g4id[j] ) {   // && bydeadwire == false
		h_BlipEvsTrkLength->Fill(proton_KE, trk_length[j]);
		if( trk_length[j] > 5.0) {
		  //std::cout<<"\n\n\n"<<run<<" "<<subrun<<" "<<event;
		  for(auto hi : map_blipID_hitID[k] ) {   // if( hit_rms[l] > 10 )
		      //std::coblipana_protons_0to200MeV_overlay_20250826ut<<"\n\nWire#: "<<hit_wire[hi]<<"\nTime tick: "<<2400+hit_driftT[hi]<<"\nPlane "<<hit_plane[hi]<<"\nHit width: "<<hit_rms[hi]<<"\nBlip ID: "<<k;
		  }//end hit loop
		}
	      }
	    }//end blip loop 
          } else { h_TrueLengthvsTrkLength->Fill(proton_len, 0); }
        }//end track loop

	if( trk == false ) {
	  total_pr_notrks++; total_notrk_E += proton_KE;
	  if( blip==false ) { h_BlipXorTrack->Fill(proton_KE); }
	}
	if( trk==true || blip==true ) { h_BlipOrTrackEff->Fill(proton_KE); }

	h_ProtonEvsTimespan->Fill(max_timespan, proton_KE);

      }//end proton loop


    }//end events loop
    
    float p_n_ratio = float(total_second_prs)/total_neutrons;

    h_TrkEfficiency->Divide(h_NumTrksPerE, h_TrueProtonE);
    h_BlipEfficiency->Divide(h_NumBlipsPerE, h_TrueProtonE);
    h_BlipOrTrackEff->Divide(h_TrueProtonE);
    h_BlipXorTrack->Divide(h_TrueProtonE);
    h_2and2wires->Divide(h_2and2wires, h_ProtonBlipE);
    h_2and3wires->Divide(h_2and3wires, h_ProtonBlipE);
    h_3and3wires->Divide(h_3and3wires, h_ProtonBlipE);
    h_3and4wires->Divide(h_3and4wires, h_ProtonBlipE);

    std::cout<<"\n\n\nNumber of secondary protons in sample: "<<total_second_prs<<
    " protons \nProtons produced per neutron: "<<p_n_ratio<<
    " protons/neutrons \nTotal reconstructed proton tracks: "<<total_pr_trks<<
    " tracks \nTotal reconstructed proton blips: "<<total_pr_blips<<
    " blips \nMinimum proton energy producing a track: "<<min_trk_E<<
    " MeV \nAverage proton energy producing a track: "<<total_trk_E/total_pr_trks<<
    " MeV \nAverage proton energy that does NOT produce a track: "<<total_notrk_E/total_pr_notrks<<
    " MeV \nAverage blip energy: "<<total_blip_E/total_pr_blips<<
    " MeV \nBlips by dead wires: "<<blips_by_death<<
    " blips\nMaximum wires per blip: "<<max_nwires<<" wires\n"; 

    merge_plot(h_TrkEfficiency, "Tracks", h_BlipEfficiency, "Blips");   // can't run this twice in same macro b/c only 1 canvas
    //merge_plot(h_NumTrksPerE, "Track Energy", h_ProtonBlipE, "Blip Energy");   // can't run this twice in same macro b/c only 1 canvas

    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    h_BlipOrTrackEff->DrawCopy("hist");

    TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
    h_BlipXorTrack->DrawCopy("hist");

    TCanvas *c4 = new TCanvas("c3", "Canvas", 800, 600);
    h_2and2wires->SetLineColor(kRed-7);
    h_2and2wires->SetFillColor(kRed-7);
    TH1 *h_2x2_copy = h_2and2wires->DrawCopy();
    h_2and3wires->SetLineColor(kAzure-5);
    h_2and3wires->SetFillColor(kAzure-5);
    TH1 *h_2x3_copy = h_2and3wires->DrawCopy("same");
    h_3and3wires->SetLineColor(kTeal-6);
    h_3and3wires->SetFillColor(kTeal-6);
    TH1 *h_3x3_copy = h_3and3wires->DrawCopy("same");
    h_3and4wires->SetLineColor(kViolet-5);
    h_3and4wires->SetFillColor(kViolet-5);
    TH1 *h_3x4_copy = h_3and4wires->DrawCopy("same");
    auto legend = new TLegend(0.42,0.2,0.58,0.4);
    legend->AddEntry(h_2x2_copy,"2 wires each");
    legend->AddEntry(h_2x3_copy,"2 and 3 wires");
    legend->AddEntry(h_3x3_copy,"3 wires each");
    legend->AddEntry(h_3x4_copy,"3 and 4 wires");
    legend->Draw();

    //===========================================
    // Write the output and close the file
    outfile->Write();
    outfile->Close();


}
