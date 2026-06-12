  #include "variables.h"
  #include "ProcessHandler.h"
  #include "../analyze_neutrons/accessory_macros/merge_plot.C"
  #include <iostream>
  #include <cmath>
  #include <vector>
  #include <map>
  #include <algorithm>
  #include "TMath.h"
  #include "TROOT.h"
  #include "TVectorT.h"
  #include "TArrow.h"
 
  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;

  // Input file name
  //std::string fFileName = "BlipAna_v08_00_00_67_overlay_Tl208g_G10_24struts_3gperStrutperEVD_systematic_variations_LoosenConfig_ALL_STRUTS_merge_PI_PII.root";
  std::string fFileName = "blipana_g10_1strut_100gammas_overlay_lowthresh_20250731.root";
  //std::string fFileName = "blipana_eminus_0to10MeV_overlay_lowthresh_20250731.root";

  // only deal with contained neutrons in calculation and plotting
  TH1D*   h_nblips_vs_vtxDist;
  TH1D*   h_maxWires;
  TH2D*   h_blipZY;
  TH1D*   h_numTrks;
  TH1D*   h_trackLength;
  TH1D*   h_blipE;
  TH1D*   h_gammaBlipE;
  TH1D*   h_gamBlips_vs_vtxDist;
  TH1D*   h_cos_momentum;
  TH1D*   h_cos_recoDirection;
  TH1D*   h_cos_comparison;
  TH2D*   h_cosP_vs_energy;
  TH2D*   h_cosDIR_vs_energy;
  TH2D*   h_angP_vs_energy;
  TH2D*   h_angDIR_vs_energy;
  TH1D*   h_angP;
  TH1D*   h_angDIR;
  TH2D*   h_res_vs_linearity;

  void gamma_analysis(){
 
    // Read in the TFile and grab the TTree
    TFile* file = new TFile(fFileName.c_str(),"READ");
    TTree* tree = (TTree*)file->Get("blipana/anatree");
    std::cout<<"Reading data file "<<fFileName<<"\n";

    // read angle resolution data from Python-generated ROOT file
    TFile* angleFile = new TFile("angles_g.root", "READ");
    TTree* angleTree = (TTree*)angleFile->Get("tree");
    double angle; std::vector<double> angles; // = nullptr;
    angleTree->SetBranchAddress("angle", &angle);

    int nEntries = angleTree->GetEntries();
    for(int i = 0; i < nEntries; i++) {
        angleTree->GetEntry(i);
        angles.push_back(angle);
    }
    std::cout << "Loaded " << angles.size() << " angle values from angles_g.root\n";

    // Set branch address for all variables you plan to use
    // (if a branch address isn't set, the variables will 
    // not be filled when each event is read). For a description 
    // of each variable, see the comments in 'variables.h'
    tree->SetBranchAddress("event",           &event);
    tree->SetBranchAddress("run",             &run);
    tree->SetBranchAddress("subrun",          &subrun);
    tree->SetBranchAddress("nparticles",      &nparticles);
    tree->SetBranchAddress("part_g4id",       &part_g4id);
    tree->SetBranchAddress("part_process",    &part_process);
    tree->SetBranchAddress("part_isPrimary",  &part_isPrimary);
    tree->SetBranchAddress("part_pdg",        &part_pdg);
    tree->SetBranchAddress("part_mother",     &part_mother);
    tree->SetBranchAddress("part_trackID",    &part_trackID);
    tree->SetBranchAddress("part_KE",         &part_KE);
    tree->SetBranchAddress("part_pathlen",    &part_pathlen);
    tree->SetBranchAddress("part_depEnergy",  &part_depEnergy);
    tree->SetBranchAddress("part_nDaughters", &part_nDaughters);
    tree->SetBranchAddress("part_Px",	      &part_Px);
    tree->SetBranchAddress("part_Py",	      &part_Py);
    tree->SetBranchAddress("part_Pz",	      &part_Pz);
    tree->SetBranchAddress("part_startPointx",&part_startPointx);
    tree->SetBranchAddress("part_startPointy",&part_startPointy);
    tree->SetBranchAddress("part_startPointz",&part_startPointz);
    tree->SetBranchAddress("part_endPointx",  &part_endPointx);
    tree->SetBranchAddress("part_endPointy",  &part_endPointy);
    tree->SetBranchAddress("part_endPointz",  &part_endPointz);
    tree->SetBranchAddress("nblips",          &nblips);
    tree->SetBranchAddress("blip_x",          &blip_x);
    tree->SetBranchAddress("blip_y",          &blip_y);
    tree->SetBranchAddress("blip_z",          &blip_z);
    tree->SetBranchAddress("blip_size",       &blip_size);
    tree->SetBranchAddress("blip_energy",     &blip_energy);
    tree->SetBranchAddress("blip_edepid",     &blip_edepid);
    tree->SetBranchAddress("blip_g4id",       &blip_g4id);
    tree->SetBranchAddress("ntrks",           &ntrks);
    tree->SetBranchAddress("trk_length",      &trk_length);
    tree->SetBranchAddress("edep_g4id",	      &edep_g4id);
    tree->SetBranchAddress("edep_g4trkid",    &edep_g4trkid);
    tree->SetBranchAddress("edep_isPrimary",  &edep_isPrimary);
    tree->SetBranchAddress("edep_pdg",	      &edep_pdg);
    tree->SetBranchAddress("nhits",           &nhits);
    tree->SetBranchAddress("hit_plane",       &hit_plane);
    tree->SetBranchAddress("hit_blipid",      &hit_blipid);

    // make output file to store plots
    TFile* outfile = new TFile("gout.root", "recreate");

    // Make histograms declared above
    h_nblips_vs_vtxDist	    = new TH1D("nblips_vs_vtxDist","All Blips Dist to reco vtx;Dist to reco vtx [cm]; Nblips / r^2",25,0,75);
    h_gamBlips_vs_vtxDist   = new TH1D("gamBlips_vs_vtxDist","Gamma Blips Dist to reco vtx;Distance to reco vtx [cm]; Nblips / r^2",25,0,75);
    h_maxWires		    = new TH1D("maxWires","Max Wires per Blip on one Plane;# Wires",18,0,18);
    h_blipZY                = new TH2D("blipZY","Blip Z vs Y;z;y",224,-40,1080,50,-150,150);
    h_numTrks               = new TH1D("numTrks","Number of Tracks;# tracks",20,0,20);
    h_trackLength	    = new TH1D("trackLengths","Track Lengths;length [cm]",150,-10,300);
    h_blipE		    = new TH1D("blipE","Blip Energy; E (MeV)",50,0,1);
    h_gammaBlipE	    = new TH1D("gammaBlipE","Gamma Blip Energy; E (MeV)",50,0,5);
    // cos of the angle b/w proton direction and inferred neutron direction (vector between nu vtx and proton location)
    h_cos_momentum	    = new TH1D("cos_momentum","Cos(angle) b/w Initial e- Momentum and Gamma Direction;cos(ang);# blips",51,-1.01,1.01);
    h_cos_recoDirection     = new TH1D("cos_recoDirection","Cos(angle) b/w Reco. e- Direction and Gamma Direction;cos(ang);# blips",51,-1.01,1.01);
    h_cos_comparison	    = new TH1D("cos_comparison","Initial e- Momentum Cos Minus e- Direction Cos;Difference in Cos(ang);# blips",20,0,2);
    //h_cosP_vs_energy	    = new TH2D("cosP_vs_energy","Cos(angle) Momentum vs Blip Energy;E [MeV];cos(ang)",10,2,3,20,-1.01,1.01);
    //h_cosDIR_vs_energy      = new TH2D("cosDIR_vs_energy","Cos(angle) Direction vs Blip Energy;E [MeV];cos(ang)",10,2,3,20,-1.01,1.01);
    h_angP_vs_energy        = new TH2D("angP_vs_energy","Momentum Angle vs True Blip Energy;E [MeV];angle [degrees]",270,0,2.7,720,0,180);
    h_angDIR_vs_energy      = new TH2D("angDIR_vs_energy","Direction Angle vs True Blip Energy;E [MeV];angle [degrees]",270,0,2.7,720,0,180);
    h_angP		    = new TH1D("angP","Momentum Angle;angle [degrees]",70,2,2.7);
    h_angDIR		    = new TH1D("angDIR","Direction Angle;angle [degrees]",70,2,2.7);
    h_res_vs_linearity      = new TH2D("res_vs_linearity","Angular Resolution vs Linearity;linearity;angle [degrees]",100,0,2,900,0,90); 

    int numEvents = tree->GetEntries();
    int   totalBlipE=0;
    int   totalBlips=0;
    int   totalGammas=0;
    int   total_gamma_blips=0;
    int   total_electrons=0;
    int   rounds=0;
    int   rnnBlipIDX=0;
    std::vector<int> prim_gams_g4trkid;
    //std::map<int, int> map_trkid_mother, map_trkid_pdg, map_trkid_index;
    std::map<int, int> map_g4id_mother, map_g4id_pdg, map_g4id_index;
    //ProcessHandler handler( h_cos_nIne_gamma, h_cos_nCap_gamma, h_cos_eBrem_gamma, h_cos_nOther_gamma, total_nCap_gammas );

    // Loop through all events
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      //if( iEvent != 0 ) continue;
      tree->GetEntry(iEvent);

      // set blipID map     
      std::map<int,std::vector<int>> map_blipID_hitID;
      for(int ih=0; ih<nhits; ih++){ if (hit_blipid[ih] >= 0) map_blipID_hitID[hit_blipid[ih]].push_back(ih); }

      //map_trkid_mother.clear(); map_trkid_pdg.clear(); map_trkid_index.clear();
      map_g4id_mother.clear(); map_g4id_pdg.clear(); map_g4id_index.clear();

      rounds++; 
      //==========================================================
      // Loop over the true MCParticles in this event
      //std::cout<<"Event has "<<nparticles<<" MC particles\n";    
      for(int i=0; i<nparticles; i++){
	int pdg = part_pdg[i], g4id = part_g4id[i], mom = part_mother[i];  // trkid = part_trackID[i],

        // Printout this particle information to the screen 
        if( iEvent==-1 ){   //&& rounds<=20   //part_trackID[i]
          printf("%4i  G4ID: %-3i  PDG: %-3i  dL=%5.1f cm,  KE0=%7.3f MeV,  Edep=%7.3f MeV,  moth=%4i, %8s, NDaughts=%3i, Is Contained=%i, Time=%7.2f us\n",
            i, part_g4id[i], part_pdg[i], part_pathlen[i], part_KE[i], part_depEnergy[i],
	    part_mother[i], part_process->at(i).c_str(), part_nDaughters[i], part_isContained[i], part_startT[i]
          );
        }
        //map_trkid_pdg[trkid] = pdg; map_trkid_mother[trkid] = mom; map_trkid_index[trkid] = i;
        map_g4id_pdg[g4id] = pdg; map_g4id_mother[g4id] = mom; map_g4id_index[g4id] = i;

        TVector3 reco_nu_vtx( nu_reco_vtx_x, nu_reco_vtx_y, nu_reco_vtx_z );
	TVector3 part_start( part_startPointx[i], part_startPointy[i], part_startPointz[i] );
	float dist = (part_start-reco_nu_vtx).Mag();
	TVector3 part_P( part_Px[i], part_Py[i], part_Pz[i] );        
        TVector3 trace_vec( part_start-reco_nu_vtx );
        float ang_blip = part_P.Angle(trace_vec);
        float cos = TMath::Cos(ang_blip);
        // -------------------------------------------------------------
        // identify any gammas (PDG=22)
        if( part_pdg[i] == 22 ) {
	  totalGammas++;
	  //if( part_isPrimary ) prim_gams_g4trkid.push_back(edep_g4trkid);
	}
	// identify any electrons
	if( part_pdg[i] == 11 ) { total_electrons++; continue; }
      }//end loop over all MCParticles

      h_numTrks->Fill(ntrks);
      for( int i=0; i<ntrks; i++){ h_trackLength->Fill(trk_length[i]); }

      for(int i=0; i<nblips; i++){   // BLIP LOOP
	int edepid = blip_edepid[i]; //int g4trkid = edep_g4trkid[edepid];
	if( edepid < 0 || blip_g4id[i] < 0 ) { continue; }  // std::cout<<"\nedepid = "<<edepid<<", g4trkID = "<<g4trkid<<"\n"; 
	//int p = map_trkid_index[g4trkid], pdg = edep_pdg[edepid], mom_trkid = map_trkid_mother[g4trkid], mom_idx = map_trkid_index[mom_trkid];
        int p = map_g4id_index[blip_g4id[i]];    // blip_partid[i]; // this don't seem to work
	int pdg = part_pdg[p], g4id = part_g4id[p], mom_idx = map_g4id_index[part_mother[p]];
	//if( pdg != map_trkid_pdg[g4trkid] ) std::cout<<"\nNOPE, PDG = "<<pdg<<", map_PDG = "<<map_trkid_pdg[g4trkid];
	totalBlips++;
	//std::cout<<"\nBlip g4id: "<<g4id;   // g4trkid
	float energy = part_KE[p];   //edep_energy[edepid];
	totalBlipE+=energy;
        h_blipE->Fill(energy);

	float blipX=blip_x[i], blipY=blip_y[i], blipZ=blip_z[i];
	TVector3 blip_loc( blipX, blipY, blipZ );
	h_blipZY->Fill(blipZ, blipY);

        //int maxWires = std::round(blip_dw[i]/0.3); h_maxWires->Fill(maxWires);   // in cm so convert to Nwires by /0.3 cm	

	if( pdg == 11 && part_isPrimary[mom_idx] ) { //&& part_isPrimary[p] ) {
	    h_gammaBlipE->Fill(energy);
	    total_gamma_blips++;

	    TVector3 phot_start(part_startPointx[mom_idx], part_startPointy[mom_idx], part_startPointz[mom_idx]);
	    TVector3 elec_start(part_startPointx[p], part_startPointy[p], part_startPointz[p]);
	    TVector3 elec_end(part_endPointx[p], part_endPointy[p], part_endPointz[p]);
	    //std::cout<<"\nMagnitude of part_start - blip_loc is "<<(part_start-blip_loc).Mag();

	    TVector3 gam_dir( elec_start - phot_start );
	    TVector3 elec_dir( elec_end - elec_start );
 	    TVector3 elec_P( part_Px[p], part_Py[p], part_Pz[p] );
	    float ang_P = elec_P.Angle(gam_dir), ang_dir = elec_dir.Angle(gam_dir);
	    float cosP = TMath::Cos(ang_P), cosDIR = TMath::Cos(ang_dir);
            //float abs_cos = fabs(cos);
            h_cos_momentum->Fill(cosP);
	    h_cos_recoDirection->Fill(cosDIR);
	    h_cos_comparison->Fill(cosP - cosDIR);
	    //if( energy > 2 ) h_cosP_vs_energy->Fill(energy, cosP);
	    //if( energy > 2 ) h_cosDIR_vs_energy->Fill(energy, cosDIR);
	    if( fabs(part_KE[mom_idx] - 2.6) < 0.1 ) {   // round(part_KE[mom_idx] * 10) / 10.0 == 2.6
	      h_angP_vs_energy->Fill(energy, TMath::RadToDeg()*ang_P);
              h_angDIR_vs_energy->Fill(energy, TMath::RadToDeg()*ang_dir);
              if( energy > 2 ) h_angP->Fill(TMath::RadToDeg()*ang_P);
              if( energy > 2 ) h_angDIR->Fill(TMath::RadToDeg()*ang_dir);
	      //std::cout<<"2.6 MeV gamma->electron process: "<<part_process->at(p)<<std::endl;
	    }

	    std::vector<int> blip_plane; int Nhits=0;
	    for(auto hi : map_blipID_hitID[i] ) { //hit loop
	      Nhits++; blip_plane.push_back( hit_plane[hi] );
	    }
	    std::sort(blip_plane.begin(), blip_plane.end());
	    int uniqueCount = std::unique(blip_plane.begin(), blip_plane.end()) - blip_plane.begin();
	    if( Nhits>2 && uniqueCount!=1 && rnnBlipIDX < angles.size() ) {
	      h_res_vs_linearity->Fill( (cosP - cosDIR), angles[rnnBlipIDX] ); rnnBlipIDX++; 
	    }

            // specific neutron process
	    // all gamma->electron processes are "compt", "phot", "conv" (pp)
	    //if( part_process->at(p).compare("phot") ) { handler.handleProcess(part_process->at(mom_idx), cos, p_mom, energy); }
	}
      }//end blip loop

    }//end loop over events

    //handler.printEnergyStatistics();

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    h_res_vs_linearity->DrawCopy("hist");

    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    h_gammaBlipE->DrawCopy("hist");

    TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
    h_cos_momentum->DrawCopy("hist");

    TCanvas *c4 = new TCanvas("c4", "Canvas", 800, 600);
    h_cos_comparison->DrawCopy("hist");

    std::cout<<"\nNumber of events in sample: "<<numEvents<<
    " events\nTotal gammas in file: "<<totalGammas<<
    " gammas\nTotal electrons in file: "<<total_electrons<<
    " electrons\nTotal primary gamma blips: "<<total_gamma_blips<<
    " gammas\nTotal blips: "<<totalBlips<<" blips\n";

    Float_t a;
    TTree otree("evd_tree", "TTree with a single float branch for total evd"); // two branches
    TBranch* a_branch = otree.Branch("totalEvents", &a, "totalEvents/F");   //total_blips
    a = numEvents;
    otree.Fill();

    //===========================================
    // Write the output and close the file
    outfile->Write();
    outfile->Close();
    angleFile->Close();

  }
