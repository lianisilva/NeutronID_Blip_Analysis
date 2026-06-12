  ////////////////////////////////////////////////////////
  //  Blip Analysis ROOT Macro - Tutorial
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
    
  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;
  
  // Input file name
  std::string fFileName = "BlipAna_neutrons_0to200MeV_overlay_uboone_20250122.root";
  
  // Tell ROOT about the histograms we'll be making later
  TH1D*   h_TrueFSGammaE;           // individual neutron capture final-state gamma energies
  TH1D*   h_TrueFSGammaESum;        // summed energies of capture gammas
  TH2D*   h_TrueNeutronEvsProtonE;  // neutron kinetic energy vs resultant proton KE
  TH2D*   h_TrueProtonEvsLength;    // Proton energy vs length


  // ********************************************************
  // This is the main function that is called when the call is
  // made to ROOT to process the macro. This must match the name
  // of the file ('analyze_neutrons.C')
  // *********************************************************
  void Original(){
  
    // Read in the TFile and grab the TTree
    TFile* file = new TFile(fFileName.c_str(),"READ");
    TTree* tree = (TTree*)file->Get("blipana/anatree");
    std::cout<<"Reading data file "<<fFileName<<"\n";
    
    // Set branch address for all variables you plan to use
    // (if a branch address isn't set, the variables will 
    // not be filled when each event is read). For a description 
    // of each variable, see the comments in 'variables.h'
    tree->SetBranchAddress("nparticles",      &nparticles);
    tree->SetBranchAddress("part_isPrimary",  &part_isPrimary);
    tree->SetBranchAddress("part_trackID",    &part_trackID);
    tree->SetBranchAddress("part_pdg",    &part_pdg);
    tree->SetBranchAddress("part_mother",     &part_mother);
    tree->SetBranchAddress("part_KE",         &part_KE);
    tree->SetBranchAddress("part_pathlen",    &part_pathlen);
    tree->SetBranchAddress("part_depEnergy",  &part_depEnergy);
    tree->SetBranchAddress("part_process",    &part_process);
    tree->SetBranchAddress("nedeps",          &nedeps);
    tree->SetBranchAddress("edep_g4id",       &edep_g4id);
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
    tree->SetBranchAddress("blip_energy",     &blip_energy);
    tree->SetBranchAddress("blip_edepid",     &blip_edepid);              
    tree->SetBranchAddress("blip_nplanes",    &blip_nplanes);
    
    // make output file to store plots
    TFile* outfile = new TFile("out.root", "recreate");
  
    // Make histograms declared above
    h_TrueFSGammaE          = new TH1D("TrueFSGammaE","Neutron capture final-state gammas;true energy [MeV]",100,0,10);
    h_TrueFSGammaESum       = new TH1D("TrueFSGammaESum","Neutron capture final-state gammas;summed true energy [MeV]",200,0,20);
    h_TrueNeutronEvsProtonE = new TH2D("TrueNeutronEvsProtonE",";Neutron kinetic energy [MeV];FS proton kinetic energy [MeV]",100,0,200,100,0,200);
    h_TrueProtonEvsLength   = new TH2D("TrueProtonEvsLength",";Proton energy [MeV];True pathlength [cm]",100,0,200,100,0,100);

    // set any 2D histograms to "colz" option (bin values 
    // displayed with rainbow color scale)
    h_TrueNeutronEvsProtonE->SetOption("colz");
    h_TrueProtonEvsLength->SetOption("colz");

    // Count the entries in this TTree file. Each entry corresponds
    // to one readout of the detector, referred to as an "event"
    int numEvents = tree->GetEntries();
    std::cout<<"There are "<<numEvents<<" events in this file.\n";
   
    int   total_second_prs = 0;
 
    //====================================================
    // Loop through the events. For each one, we will call
    // the "GetEntry" command, which tells ROOT to fill all
    // the variables with set branch addresses (don't worry, 
    // this will all make more sense a bit later on...)
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      
      tree->GetEntry(iEvent);
      
      std::cout<<"=======================================\n";
      std::cout<<"Reading event "<<iEvent<<", which has "<<nblips<<" blips\n";
      

      //==========================================================
      // Loop over the true MCParticles in this event
      std::cout<<"Event has "<<nparticles<<" MC particles\n";
      
      // keep track of info about the primary neutron
      int   neutron_index = -9;
      int   neutron_g4id  = -9;
      int   neutron_KE    = -9;
      float neutron_FSgamma_sumE = 0;

      // keep track of info about the secondary proton (if produced)
      int   proton_num    = 0;
      int   proton_index  = -9;

      for(int i=0; i<nparticles; i++){

        // -------------------------------------------------------------
        // Printout this particle information to the screen. As the analysis
        // progresses past these initial validation/debugging stages,
        // feel free to comment this section out to speed up the event
        // loop process.
        printf("  %5i  G4ID: %-8i PDG: %-10i  dL=%7.1f cm,  KE0=%9.3f MeV,  Edep=%9.3f MeV,  moth=%5i, %12s, NDaughters=%i\n",
          i,
          part_trackID[i],
          part_pdg[i],
          part_pathlen[i],
          part_KE[i],
          part_depEnergy[i],
          part_mother[i],
          part_process->at(i).c_str(),
          part_nDaughters[i]
        );


        // -------------------------------------------------------------
        // identify the primary neutron
        if( part_isPrimary[i] && part_pdg[i]==2112 && neutron_index <= 0){
          neutron_index = i;
          neutron_g4id  = part_trackID[i];
          neutron_KE    = part_KE[i];
        }

        // -------------------------------------------------------------
        // IF the neutron was already found...
        if( neutron_index >= 0 ) {
        
          // identify any gammas (PDG=22) that come from that neutron capture
          if( part_pdg[i] == 22 && part_mother[i] == neutron_g4id ){
            h_TrueFSGammaE->Fill(part_KE[i]);
            neutron_FSgamma_sumE += part_KE[i];
          }

          // identify any protons produced by neutron inelastic collisions
          if( part_pdg[i] == 22 && part_mother[i] == neutron_g4id && part_process->at(i).compare("nCapture") == 0 ){
            proton_num++;
            proton_index  = i;
          }
          
        }

        // sum the number of protons in over all events
        total_second_prs += proton_num;

      }//endloop over all MCParticles
    
      h_TrueFSGammaESum ->Fill(neutron_FSgamma_sumE);
     
      // check that multiple FS protons weren't produced
      if( proton_num == 1 ) {
        float proton_KE   = part_KE[proton_index];
        float proton_len  = part_pathlen[proton_index];
        h_TrueNeutronEvsProtonE->Fill(neutron_KE,proton_KE);
        // if proton doesn't leave the AV, plot its pathlength
        //if( part_isContained[proton_index] ) {
          h_TrueProtonEvsLength->Fill(proton_KE,proton_len);
        //}
      }


      // Check if the proton got reconstructed in as a
      // 3D track (produced by Pandora, arXiv:1708.03135)
      for(int i=0; i<ntrks; i++){
        if( trk_g4id[i] == part_trackID[proton_index] ) {
          std::cout<<"Proton was reconstructed into a track\n";
        }
      }

      
    }//endloop on events
   
    cout<<"Total protons: "<<total_second_prs<<"\n";
 
    //===========================================
    // Write the output and close the file
    outfile->Write();
    outfile->Close();

}
