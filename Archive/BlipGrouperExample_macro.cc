  ////////////////////////////////////////////////////////
  //  Blip Analysis ROOT Macro - Tutorial
  //
  //  This ROOT/C++ analysis code reads in a TTree produced
  //  by the MicroBooNE 'BlipAna' module and makes some plots
  //  of different blip-related metrics. Grouping of unmatched
  //  blips ('hit clusters') on wire planes is also performed.
  //
  //  This macro must be in the same directory as the header file
  //  "variables.h" as well as a BlipAna ROOT file. To run, type
  //  this command in the terminal:
  //  
  //  > root BlipGrouperExample_macro.cc
  //
  //  W. Foreman
  //  June 2023
  ////////////////////////////////////////////////////////
  
  // We first include this header file, which contains all the 
  // variables stored in the ROOT TTree we will be reading from.
  #include "variables.h"
    
  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;
  
  // Input file name
  //std::string fFileName = "BlipAna_20230110_Data_OverlayNull.root";
  std::string fFileName = "../files/BlipAna_20230125_Data_Run3_Unbiased.root";
  
  // Tell ROOT about the histograms we'll be making later
  TH1D*   hBlipEnergy;
  TH1D*   hBlipEnergy_SS; // 'single-site'
  TH1D*   hBlipEnergy_MS; // 'multi-site'
  TH1D*    hBlipEnergy_HotSpot;
  TH1D*   hBlipRegionMult;
  TH1D*   hBlipRegionEnergy;

  
  TH1D*   hBlipDistance2D;
  TH1D*   hBlipDistance2D_EnergyCut;
  
  // ********************************************************
  // This is the main function that is called when the call is
  // made to ROOT to process the macro. This must match the name
  // of the file ('BlipGrouperExample_macro.cc')
  // *********************************************************
  void BlipGrouperExample_macro(){
  
    // Read in the TFile and grab the TTree
    TFile* file = new TFile(fFileName.c_str(),"READ");
    TTree* tree = (TTree*)file->Get("blipanaTrkMask/anatree");
    std::cout<<"Reading data file "<<fFileName<<"\n";
    
    // Set branch address for all variables you plan to use
    // (if a branch address isn't set, the variables will 
    // not be filled when each event is read). For a description 
    // of each variable, see the comments in 'variables.h'
    tree->SetBranchAddress("event",           &event);
    tree->SetBranchAddress("nclusts",         &nclusts);
    tree->SetBranchAddress("clust_plane",     &clust_plane);              
    tree->SetBranchAddress("clust_startwire", &clust_startwire);                
    //tree->SetBranchAddress("clust_endwire",   &clust_endwire);                
    tree->SetBranchAddress("clust_charge",    &clust_charge);            
    tree->SetBranchAddress("clust_time",      &clust_time);                
    //tree->SetBranchAddress("clust_edepid",    &clust_edepid);            
    //tree->SetBranchAddress("clust_blipid",    &clust_blipid);            
    //tree->SetBranchAddress("clust_bydeadwire",&clust_bydeadwire);
    tree->SetBranchAddress("clust_nwires",    &clust_nwires);
    tree->SetBranchAddress("nblips",          &nblips);                       
    //tree->SetBranchAddress("blip_x",          &blip_x);                        
    tree->SetBranchAddress("blip_y",          &blip_y);                        
    tree->SetBranchAddress("blip_z",          &blip_z);                       
    tree->SetBranchAddress("blip_energy",     &blip_energy);
    //tree->SetBranchAddress("blip_charge",     &blip_charge);              
    //tree->SetBranchAddress("blip_edepid",     &blip_edepid);              
    tree->SetBranchAddress("blip_nplanes",    &blip_nplanes);
    //tree->SetBranchAddress("blip_yzcorr",     &blip_yzcorr);
    tree->SetBranchAddress("blip_proxtrkdist",&blip_proxtrkdist);
    //tree->SetBranchAddress("blip_pl0_clustid",&blip_pl0_clustid);
    //tree->SetBranchAddress("blip_pl1_clustid",&blip_pl1_clustid);
    tree->SetBranchAddress("blip_pl2_clustid",&blip_pl2_clustid);
    
    // make output file to store plots
    TFile* outfile = new TFile("out.root", "recreate");
  
    // Make histograms declared above
    hBlipEnergy       = new TH1D("BlipEnergy",      "All blips (trkdist > 30cm);Energy [MeV]",100,0,5);
    hBlipEnergy_SS    = new TH1D("BlipEnergy_SS",   "Single-site blips (trkdist > 30cm);Energy [MeV]",100,0,5);
    hBlipEnergy_MS    = new TH1D("BlipEnergy_MS",   "Multi-site blips (trkdist > 30cm);Energy [MeV]",100,0,5);
    hBlipEnergy_HotSpot    = new TH1D("BlipEnergy_HotSpot",   "Single-site blips (<30cm) in hotspot;Energy [MeV]",100,0,5);
    hBlipRegionMult   = new TH1D("BlipRegionMult",  "All blips (trkdist > 30cm);Num 2D clusters within radius",50,0,50);
    hBlipRegionEnergy = new TH1D("BlipRegionEnergy","All blips (trkdist > 30cm);Summed cluster energy within radius [MeV]",100,0,20);
    
    hBlipDistance2D   = new TH1D("BlipDistance2D",  "Nearest 2D cluster on coll. plane;2D quasi-distance [cm]",150,0,150);
    hBlipDistance2D_EnergyCut   = new TH1D("BlipDistance2D_EnergyCut",  "Nearest 2D cluster on coll. plane (E cut);2D quasi-distance [cm]",150,0,150);

    const float samplePeriod = 0.5; // 0.5 us per time-tick
    const float wireSpacing  = 0.3; // spacing between adjacent wires [cm]
    const float driftSpeed   = 0.1098; // cm/us

    // Count the entries in this TTree file. Each entry corresponds
    // to one readout of the detector, referred to as an "event"
    int numEvents = tree->GetEntries();
    std::cout<<"There are "<<numEvents<<" events in this file.\n";
    
    //====================================================
    // Loop through the events. For each one, we will call
    // the "GetEntry" command, which tells ROOT to fill all
    // the variables with set branch addresses (don't worry, 
    // this will all make more sense a bit later on...)
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      
      tree->GetEntry(iEvent);
      
      if( iEvent % 100 == 0 ){ 
      std::cout<<"=======================================\n";
      std::cout<<"Reading event "<<iEvent<<", which has "<<nblips<<" blips\n";
      }
      

      //==========================================================
      // We will start by looping over ALL the collection-plane
      // hit clusters in this event and calculate the quasi-position 
      // using 'wire-time' coordinates. This is to avoid repeating these
      // calculations hundreds of times later on.
      //
      // We can convert the time and wire number to a proxy for position
      // coordinates by considering:
      //    - electron drift speed (0.1098 cm/microsecond)
      //    - wire-to-wire spacing (0.3 cm)
      //
      std::map<int,float> map_clust_wireCoord;
      std::map<int,float> map_clust_timeCoord;
      for(int j=0; j<nclusts; j++){
        if( clust_plane[j] != 2 ) continue;
        map_clust_wireCoord[j] = wireSpacing * (clust_startwire[j]+clust_nwires[j]);
        map_clust_timeCoord[j] = driftSpeed * samplePeriod * clust_time[j];
      }
      

      //==========================================================
      // We will now look though all the 3D plane-matched 'blips' in this
      // event. A blip is made when clusters of hits on at least
      // two wire planes are matched up in time, indicating that 
      // they are formed by a common group of drifted ionization.
      for(int i=0; i<nblips; i++){
      
        // Let's only consider blips that are at least 30cm away from
        // the nearest 'long' track (>5cm)
        if( blip_proxtrkdist[i] < 30 ) continue;
        //std::cout<<"  * found a blip "<<blip_proxtrkdist[i]<<" cm from nearest trk\n";
        
        // Plot the energy
        hBlipEnergy->Fill(blip_energy[i]);

        // Is this in our representative hotspot?
        bool hotspot = false;
        if( blip_z[i] > 330 && blip_z[i] < 360 && blip_y[i] < -95 ) hotspot = true;

        // Find the wire number and time of the hit cluster on the 
        // collection plane for blip 'i'
        int   clustID   = blip_pl2_clustid[i];

        // use the handy map made earlier to get the 2D wire-time position
        float blip_wireCoord = map_clust_wireCoord[clustID];
        float blip_timeCoord = map_clust_timeCoord[clustID];
        //std::cout<<"   wire= "<<blip_wireCoord<<", time= "<<blip_timeCoord<<"\n";

        //----------------------------------------------------------
        // Now loop through *all* the collection plane clusters and
        // find clusters that are nearby this blip in 2D 'wire-time' 
        // space. We can convert the time and wire number to quasi-distance
        // coordinates by considering:
        //    - electron drift speed (0.1098 cm/microsecond)
        //    - wire-to-wire spacing (0.3 cm)
        // Also keep track of the clust multiplicity and total sum energy.
        int   nclustsRegion = 0;
        float energyRegion  = 0;
        float closestDist2D = 9e9;
        float closestDist2D_cut = 9e9;
        for(int j=0; j<nclusts; j++){
          
          // skip if non-collection plane
          if( clust_plane[j] != 2 ) continue;
          if( j == clustID )        continue;
        
          float clust_energy = clust_charge[j]/0.584 * (23.6e-6);

          // find distance to our 3D matched blip from above
          // 'dw' = wire-based distance [cm]
          // 'dt' = time-based distance [cm]
          float dw = blip_wireCoord - map_clust_wireCoord[j];
          float dt = blip_timeCoord - map_clust_timeCoord[j];
          
          if( dw > 30 ) continue;
          if( dt > 30 ) continue;

          float dist = sqrt( pow(dw,2) + pow(dt,2) );
          //std::cout<<"       -- wire/time "<<map_clust_wireCoord[j]<<" "<<map_clust_timeCoord[j]<<" --> dw/dt "<<dw<<" "<<dt<<" --> dist = "<<dist<<"\n";
          
          if( dist < closestDist2D ) closestDist2D = dist;

          if( blip_energy[i] > 0.5 && clust_energy > 0.5 &&
              dist < closestDist2D_cut ) closestDist2D_cut = dist;

          if( dist < 30 ) {
            // We found a cluster < 30cm from the central blip!
            nclustsRegion++;
            // convert the cluster charge to "energy" by assuming
            // an electron-like recombination survival (R=0.584)
            // and multiplying by 23.6 eV (ionization work function)
            energyRegion += clust_energy;

          }

        }//endloop over the 2D clusters
        
        if( closestDist2D < 1e3 )     hBlipDistance2D           ->Fill(closestDist2D);
        if( closestDist2D_cut < 1e3 ) hBlipDistance2D_EnergyCut ->Fill(closestDist2D_cut);
        
        // fill histograms
        hBlipRegionMult   ->Fill(nclustsRegion);
        hBlipRegionEnergy ->Fill(energyRegion);
        
        // single-site or multi-site blip?
        if( nclustsRegion == 0 ) hBlipEnergy_SS->Fill(blip_energy[i]);
        else                     hBlipEnergy_MS->Fill(blip_energy[i]);

        // hotspot
        if( hotspot && nclustsRegion == 0 ) hBlipEnergy_HotSpot->Fill(blip_energy[i]);
      
      }//endloop on 3D blips
    
    }//endloop on events
    
    //===========================================
    // Area-normalize the energy spectra
    hBlipEnergy_SS->Scale(1./hBlipEnergy_SS->GetEntries());
    hBlipEnergy_MS->Scale(1./hBlipEnergy_MS->GetEntries());
    hBlipEnergy_HotSpot->Scale(1./hBlipEnergy_HotSpot->GetEntries());
    
    hBlipEnergy_SS->SetOption("hist");
    hBlipEnergy_SS->GetYaxis()->SetTitle("Area-normalized entries");
    hBlipEnergy_MS->SetOption("hist");
    hBlipEnergy_MS->GetYaxis()->SetTitle("Area-normalized entries");
    hBlipEnergy_HotSpot->SetOption("hist");
    hBlipEnergy_HotSpot->GetYaxis()->SetTitle("Area-normalized entries");
    
    outfile->Write();
    outfile->Close();

}
