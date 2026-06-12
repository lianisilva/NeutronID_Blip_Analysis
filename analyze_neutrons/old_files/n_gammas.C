  #include "variables.h"
  #include "TChain.h"
  #include "TFile.h"
  #include <cmath>
  #include <cstdlib>
  #include <iomanip>
  #include <fstream>
  #include <iostream>
  #include <vector>

  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;

  // Input file name
  std::string fFileName = "BlipAna_neutrons_0to200MeV_overlay_uboone_20250129.root";

  void n_gammas(){

    // Read in the TFile and grab the TTree
    TFile* file = new TFile(fFileName.c_str(),"READ");
    TTree* tree = (TTree*)file->Get("blipana/anatree");
    std::cout<<"Reading data file "<<fFileName<<"\n";

    // Set branch address for all variables you plan to use
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
    tree->SetBranchAddress("blip_size",       &blip_size);
    tree->SetBranchAddress("blip_dw",         &blip_dw);
    tree->SetBranchAddress("blip_energy",     &blip_energy);
    tree->SetBranchAddress("blip_edepid",     &blip_edepid);
    tree->SetBranchAddress("blip_nplanes",    &blip_nplanes);
    tree->SetBranchAddress("blip_pl0_clustid",&blip_pl0_clustid);
    tree->SetBranchAddress("blip_pl1_clustid",&blip_pl1_clustid);
    tree->SetBranchAddress("blip_pl2_clustid",&blip_pl2_clustid);
    tree->SetBranchAddress("clust_nwires",    &clust_nwires);
    tree->SetBranchAddress("clust_timespan",  &clust_timespan);
    tree->SetBranchAddress("clust_bydeadwire",&clust_bydeadwire);

    // Count the entries in this TTree file.
    int numEvents = tree->GetEntries();

    // make output tree to store data
    TFile* outfile = new TFile("info_file.root", "recreate");
    TTree *T = new TTree("T", "gamma_info");
    std::vector<std::vector<float>> gamma_info;
    T->Branch("gamma_info", &gamma_info);

    //====================================================
    // Loop through the events and call the "GetEntry" command
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){

      tree->GetEntry(iEvent);

      std::vector<int> neutron_indices;
      std::vector<int> gamma_indices;

      for(int i=0; i<nparticles; i++){
        // sum over ALL neutrons in the event and save their information
        if ( part_pdg[i]==2112 && part_isPrimary[i] ){ neutron_indices.push_back(i); }
        // identify any protons produced by neutron capture
        if( part_pdg[i] == 22 && part_process->at(i).compare("nCapture") == 0 ) {
          gamma_indices.push_back(i);
        }
      }

      // loop through all the signal gammas in the event
      for(int i=0; i<static_cast<int>(gamma_indices.size()); i++) {
        for(int j=0; j<nblips; j++){
          if( blip_g4id[j] == part_g4id[gamma_indices[i]] ) {
	    for(int k=0; k<ntrks; k++) {
	      if( trk_g4id[k] == part_mother[gamma_indices[i]] ) {
		gamma_info[0].push_back(trk_endx[k]);
		gamma_info[1].push_back(blip_x[j]);
	      }
	    }//end trk loop
          }
        }//end blip loop
      }//end proton loop
    }//end event loop

    T->Fill();
    T->Print();
    outfile->Write();
    outfile->Close();

  }


