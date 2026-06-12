  #include "variables.h"
  #include "merge_plot.C"

  #include "TMath.h"

  #include <iostream>
  #include <cmath>
  #include <vector>
  #include "TROOT.h"
  #include "TVectorT.h"
  #include <algorithm>
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
 
  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;
  
  // Input file name
  std::string fFileName = "blipana_genie_fullosc_overlay_run3_20250304.root";

  void cookUpQuick(){
  
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
    tree->SetBranchAddress("nhits",           &nhits);
    tree->SetBranchAddress("hit_blipid",      &hit_blipid);
    tree->SetBranchAddress("hit_plane",       &hit_plane);
    tree->SetBranchAddress("hit_wire",        &hit_wire); 
    tree->SetBranchAddress("hit_amp",         &hit_amp);
    tree->SetBranchAddress("hit_peakT",	      &hit_peakT);
    tree->SetBranchAddress("hit_driftT",      &hit_driftT);
    tree->SetBranchAddress("hit_rms",	      &hit_rms); 

    // make output file to store plots
    //TFile* outfile = new TFile("CUQ.root", "recreate");
  
    // Make histograms declared above
    TH1D* h_ProtonBlipE	    = new TH1D("ProtonBlipE","Blip energy [MeV]",100,0,100);
    TH1D* h_NumTrksPerE	    = new TH1D("NumTrksPerE","Number of Protons Made Into Tracks vs Energy;Proton energy [MeV];Number of Tracks",100,0,200);

    // Count the entries in this TTree file. Each entry corresponds
    // to one readout of the detector, referred to as an "event"
    int numEvents = tree->GetEntries();
    //std::cout<<"There are "<<numEvents<<" events in this file.\n";
    
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      tree->GetEntry(iEvent);      
      int proton_num=0;
      std::vector<int> proton_indices;

      for(int i=0; i<nparticles; i++){
        // identify any protons produced by neutron inelastic collisions that do not leave the AV 
	if( part_pdg[i] == 2212 && part_isContained[i] && part_process->at(i).compare("neutronInelastic") == 0 ) { // part_isPrimary[i] ) { 
          proton_num++;
          proton_indices.push_back(i);
        }
      } //end loop over all MCParticles

      // loop through all the protons in the event
	for(int i=0; i<proton_num; i++) {
	  float proton_KE = part_KE[proton_indices[i]]; 
	  int proton_g4id = part_g4id[proton_indices[i]];

	  // Check if the proton got reconstructed as a blip
	  for(int j=0; j<nblips; j++){
	      if( blip_g4id[j] == proton_g4id ) { 
	          h_ProtonBlipE->Fill(proton_KE);
	      }
	  }//end blip loop

	  // Check if the proton got reconstructed as a 3D track
	  for(int j=0; j<ntrks; j++) {
	      if( trk_g4id[j] == proton_g4id ) {
		// Add this right before the problematic line:
	  //std::cout << "About to fill h_NumTrksPerE with proton_KE = " << proton_KE << std::endl;
	  //std::cout << "h_NumTrksPerE pointer address: " << h_NumTrksPerE << std::endl;
	  
	  if (h_NumTrksPerE == nullptr) {
	      //std::cout << "ERROR: h_NumTrksPerE is null!" << std::endl;
	      continue;
	  }
	            h_NumTrksPerE->Fill(proton_KE);
	        }
	    }//end track loop
	}//end proton loop

    }//end events loop

    merge_plot(h_NumTrksPerE, "Tracks", h_ProtonBlipE, "Blips");   // can't run this twice in same macro b/c only 1 canvas

    //TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    //h_BlipEvsAvgWires->DrawCopy("hist");

    //TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    //h_BlipEvsTrkLength->DrawCopy("hist");

    //===========================================
    // Write the output and close the file
    //outfile->Write();
    //outfile->Close();


}
