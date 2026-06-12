  #include "variables.h"
  #include <iostream>
  #include "TFile.h"
  #include "TTree.h"
  #include "TClass.h"
  #include "TH1D.h"
  #include "TCanvas.h"
  #include <cmath>
  #include <fstream>
  #include <vector>
  using namespace std;

  // Max number of events to iterate (set negative = read all events)
  int maxEvents = -2000;

  // Input file name
  //std::string fFileName = "old_files/blipana_protons_0to200MeV_overlay_20250313.root";   // >200k protons
  //std::string fFileName = "blipana_genie_fullosc_overlay_run3_20250304.root";  // BNB nue
  //std::string fFileName = "blipana_eminus_0to10MeV_overlay_lowthresh_20250731.root";   // ELECTRONS
  std::string fFileName = "blipana_protons_0to200MeV_overlay_20250826.root";   // MORE protons (1mil?)
  //std::string fFileName = "blipana_MCC9_bnb_nu_overlay_mcc9_run3_20250904.root";  // BNB nu
  //std::string fFileName = "blipana_prodgenie_bnb_nu_overlay_mcc9_run3_numuCCinclusive_selection_20250926.root";  // numuCC inclusive

  void hit_info(){

    //part_process = new std::vector<std::string>();

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
    //tree->SetBranchAddress("part_process",    &part_process);
    tree->SetBranchAddress("part_Px",	      &part_Px);
    tree->SetBranchAddress("part_Py",	      &part_Py);
    tree->SetBranchAddress("part_Pz",	      &part_Pz);
    tree->SetBranchAddress("part_startPointx",&part_startPointx);
    tree->SetBranchAddress("part_startPointy",&part_startPointy);
    tree->SetBranchAddress("part_startPointz",&part_startPointz);
    tree->SetBranchAddress("part_endPointx",  &part_endPointx);
    tree->SetBranchAddress("part_endPointy",  &part_endPointy);
    tree->SetBranchAddress("part_endPointz",  &part_endPointz);
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
    tree->SetBranchAddress("blip_id",	      &blip_id);
    tree->SetBranchAddress("blip_x",          &blip_x);
    tree->SetBranchAddress("blip_y",          &blip_y);
    tree->SetBranchAddress("blip_z",          &blip_z);
    tree->SetBranchAddress("blip_size",       &blip_size);
    tree->SetBranchAddress("blip_dw",         &blip_dw);
    tree->SetBranchAddress("blip_energy",     &blip_energy);
    tree->SetBranchAddress("blip_g4id",       &blip_g4id);
    tree->SetBranchAddress("blip_edepid",     &blip_edepid);
    tree->SetBranchAddress("blip_nplanes",    &blip_nplanes);
    //tree->SetBranchAddress("blip_pl0_clustid",&blip_pl0_clustid);
    //tree->SetBranchAddress("blip_pl1_clustid",&blip_pl1_clustid);
    //tree->SetBranchAddress("blip_pl2_clustid",&blip_pl2_clustid);
    //tree->SetBranchAddress("clust_nwires",    &clust_nwires);
    //tree->SetBranchAddress("clust_timespan",  &clust_timespan);
    //tree->SetBranchAddress("clust_bydeadwire",&clust_bydeadwire);
    tree->SetBranchAddress("nhits",	      &nhits);
    tree->SetBranchAddress("hit_plane",	      &hit_plane);
    tree->SetBranchAddress("hit_wire",        &hit_wire);
    tree->SetBranchAddress("hit_driftT",      &hit_driftT); 
    tree->SetBranchAddress("hit_peakT",	      &hit_peakT);
    tree->SetBranchAddress("hit_amp",         &hit_amp);
    tree->SetBranchAddress("hit_rms",	      &hit_rms);
    tree->SetBranchAddress("hit_blipid",      &hit_blipid);
    //tree->SetBranchAddress("hit_charge",      &hit_charge);

    // Count the entries in this TTree file.
    int numEvents = tree->GetEntries();

    //if I want these to be global they have to be added to an outfile
    TH1D* h_NumProtons  = new TH1D("NumProtons","Number of Protons per Event", 150, 0, 15);
    TH1D* h_NumBlipHits	= new TH1D("NumBlipHits","Number of Hits per Blip; # Hits", 55, 0, 55);
    TH1D* h_HitsPerWire	= new TH1D("HitsPerWire","Number of Hits / Number of Wires", 100, 0, 5);   
    TH1D* h_BlipE	= new TH1D("BlipE","Blip Energy Spectrum", 75, 0, 150);
    TH1D* h_cos_comparison=new TH1D("cos_comparison","Cos b/w Initial Momentum and Particle Direction;Cos(ang);# blips",101,-1.01,1.01);
    TH1D* h_normPz	= new TH1D("normPz","Normalized Pz Used In RNN", 20, -2, 2);
    //std::vector<std::vector<float>> hit_info;
    //tree->Scan();   // prints the content on the screen

    ofstream outfile;
    outfile.open("info_hits.txt");
    outfile <<  "blip_g4id " << "hit_rms " << "hit_amp " << "hit_wire " << "hit_plane " << "hit_time " << "mom_x " << "mom_y " << "mom_z "<<"true_E\n";//"<<"real_g4id "<<"event\n";

    std::vector<int> blip_num;
    std::vector<float> rms;
    std::vector<float> amp;
    std::vector<float> wire;
    std::vector<int> plane;
    std::vector<double> time;
    std::vector<float> mom;
    std::vector<float> mom_x;
    std::vector<float> mom_y;
    std::vector<float> mom_z;
    std::vector<float> energy;
    //std::vector<int> g4id;
    //std::vector<int> event;
    std::vector<int> bad_num;

    std::map<int, int> map_g4id_pdg;
    std::map<int, int> map_g4id_index;
    std::map<int, int> blip_Nhits_map;
    std::map<int, std::map<int, std::set<float>>> blip_wire_per_plane_map;

    int blipN=0;
    int max_hits=0;
    int Ncosmics=0;
    int Nwires;
    int totalPrs=0;
    int over2=0;
    int over4=0;
    float minE=100;
    float maxE=0;

    //====================================================
    // Loop through the events and call the "GetEntry" command
    if( maxEvents < 0 ) maxEvents = numEvents;
    for(int iEvent = 0; iEvent < maxEvents; iEvent++ ){
      tree->GetEntry(iEvent);

      // set blipID map
      std::map<int,std::vector<int>> map_blipID_hitID;
      for(int ih=0; ih<nhits; ih++){ if (hit_blipid[ih] >= 0) map_blipID_hitID[hit_blipid[ih]].push_back(ih); }

      //std::vector<int> neutron_indices;
      std::vector<int> proton_indices;
      std::vector<int> pr_mother;
      std::vector<int> elec_idx;

      map_g4id_pdg.clear();
      map_g4id_index.clear();

      for(int i=0; i<nparticles; i++){

        map_g4id_pdg[part_g4id[i]] = part_pdg[i];
        map_g4id_index[part_g4id[i]] = i;

	// primary protons that do not leave the AV 
	if( part_pdg[i] == 2212 || part_isContained[i] || part_isPrimary[i] ) {   //&& part_process->at(i).compare("neutronInelastic") == 0 ) { 
	  totalPrs++; proton_indices.push_back(i);
	  //pr_mother.push_back(part_mother[i]);
	}
	// electrons (potential gamma products)
	if( part_pdg[i] == 11 && part_isPrimary[i] ) { elec_idx.push_back(i); }   // part_isContained[i]

      }//end loop over particles

      int proton_size = static_cast<int>(proton_indices.size()); 
      h_NumProtons->Fill(proton_size);
      int p;   // particle index
      //mother_idx = map_g4id_index[part_mother[p]];   // elec
      //if( part_pdg[mother_idx] == 22 && part_process->at(mother_idx).compare("neutronInelastic") == 0 && part_isContained[mother_idx])
      //============================== BLIP LOOP ==============================//
      for(int j=0; j<nblips; j++) {
	if( blip_edepid[j] < 0 ) { Ncosmics++; continue; }   // part_isContained[i] == 1 ) {
	if( map_g4id_index.count(blip_g4id[j]) > 0 ) p = map_g4id_index.at(blip_g4id[j]); else continue;
	//if ( p < 0 || p >= nparticles ) continue;
	//if( std::find( elec_idx.begin(), elec_idx.end(), p ) == elec_idx.end() ) continue;   // elec
        //if( std::find( proton_indices.begin(), proton_indices.end(), p ) == proton_indices.end() ) continue;   // proton
	if( !(part_pdg[p] == 2212) || !part_isContained[p] || !part_isPrimary[p] ) continue;
	// individual vectors for each blip
	std::vector<double> blip_time;
	std::vector<float> blip_wire0;
	std::vector<float> blip_wire1;
	std::vector<float> blip_wire2;
	std::vector<float> blip_amp;
        std::vector<float> blip_rms;
	std::vector<float> blip_plane;
	std::map<int, std::set<float>> wire_per_plane;
	float blipE = part_KE[p];
	if( blipE > maxE ) maxE = blipE;
	if( blipE < minE ) minE = blipE;
	h_BlipE->Fill(blipE);
	
	int Nhits=0;
	// for all the hit IDs in the blip ID container pertaining to nblip j
	for(auto hi : map_blipID_hitID[j] ) {
	  Nhits++;
	  energy.push_back(blipE);
	  //g4id.push_back(blip_g4id[j]);
	  //event.push_back(iEvent);
	  blip_num.push_back(blipN);
	  blip_amp.push_back(hit_amp[hi]);
	  blip_rms.push_back(hit_rms[hi]*(0.5*1.1)/3.0);   // convert to units of wire spacings
	  //blip_charge.push_back(hit_charge[hi]);
	  //wire.push_back(hit_wire[hi]);
	  if( hit_plane[hi] == 0 ) blip_wire0.push_back(hit_wire[hi]);
	  if( hit_plane[hi] == 1 ) blip_wire1.push_back(hit_wire[hi]);
	  if( hit_plane[hi] == 2 ) blip_wire2.push_back(hit_wire[hi]);
	  plane.push_back(hit_plane[hi]);
	  blip_plane.push_back(hit_plane[hi]);
	  blip_time.push_back(hit_driftT[hi]*(0.5*1.1)/3.0);   // convert to units of wire spacings
	  wire_per_plane[hit_plane[hi]].insert(hit_wire[hi]);

	  float Px = part_Px[p];
	  float Py = part_Py[p];
	  float Pz = part_Pz[p];
	  float P = std::sqrt(Px*Px + Py*Py + Pz*Pz);
	  float norm_Px = Px/P;
	  float norm_Py = Py/P;
	  float norm_Pz = Pz/P;
	  mom_x.push_back(norm_Px);   // don't standardize targets for regression
	  mom_y.push_back(norm_Py);
	  mom_z.push_back(norm_Pz);
	  mom.push_back(P);
	
	  TVector3 part_start(part_startPointx[p], part_startPointy[p], part_startPointz[p]);
	  TVector3 part_end(part_endPointx[p], part_endPointy[p], part_endPointz[p]);
	  TVector3 part_dir( part_end - part_start );
	  TVector3 part_P( part_Px[p], part_Py[p], part_Pz[p] );
	  float ang = part_P.Angle(part_dir);
	  float cos = TMath::Cos(ang);
	  h_cos_comparison->Fill(cos);
	
	}//end hit loop

	blip_wire_per_plane_map[blipN] = wire_per_plane;
	
	if( Nhits>2 ) over2++;
	if( Nhits>4 ) over4++;
	blip_Nhits_map[blipN] = Nhits;
	if (!h_NumBlipHits) {std::cerr << "ERROR: h_NumBlipHits is null\n";}
	h_NumBlipHits->Fill(Nhits);
	if( Nhits > max_hits ) max_hits = Nhits;
	
	// if all blip hits are (oddly) only on one plane
	std::sort(blip_plane.begin(), blip_plane.end());
	int uniqueCount = std::unique(blip_plane.begin(), blip_plane.end()) - blip_plane.begin(); 
	if(uniqueCount==1) bad_num.push_back(blipN);

	blipN++;
	//cout<<"\n\nblipN: "<<blipN<<" Nhits: "<<Nhits<<" blipID: "<<map_blipID_hitID[j];
	// not as precise for time as I would like
	//double time_mean = TMath::Mean(blip_time.begin(), blip_time.end());
	//double time_SD = TMath::StdDev(blip_time.begin(), blip_time.end());   // standard deviation
	//std::vector<float> time_z;	    // calculate z-score
	//for(auto& element : blip_time) time_z.push_back((element - time_mean)/time_SD);
	//time.insert(time.end(), time_z.begin(), time_z.end());

	// relative time within each blip (zero-shift scaling)
	float t_min = *std::min_element(blip_time.begin(), blip_time.end());//blip_time[0];
	for(auto& t : blip_time) {
	  time.push_back( t - t_min );
	  //float rel_t = t - t_min;
	  //time.push_back( abs(rel_t) < 1e-4 ? 0 : rel_t);
	}

        //for(auto& t : blip_time) { time.push_back( t ); }

	// max scale wire# per plane per blip
	std::vector<std::vector<float>> all_wires = { blip_wire0, blip_wire1, blip_wire2 };
	for( int pl=0; pl<3; pl++ ) {
	  std::vector<float> blip_wire = all_wires[pl]; if(blip_wire.empty()) continue;
	  auto wire_max = std::max_element( blip_wire.begin(), blip_wire.end());
	  float max_wire = *wire_max;
	  for(auto & wire_num : blip_wire) wire_num = max_wire-wire_num;
	  wire.insert(wire.end(), blip_wire.begin(), blip_wire.end());
	}

	// min-max scale rms and amp
	//float min_amp = *std::min_element(blip_amp.begin(), blip_amp.end());
	float max_amp = *std::max_element(blip_amp.begin(), blip_amp.end());
	for(auto& AMP : blip_amp) amp.push_back( AMP/max_amp );   //(AMP - min_amp)/(max_amp - min_amp) );   // can try charge again later
        //float min_rms = *std::min_element(blip_rms.begin(), blip_rms.end());
        //float max_rms = *std::max_element(blip_rms.begin(), blip_rms.end());
        for(auto& RMS : blip_rms) rms.push_back( RMS );   //(RMS - min_rms)/(max_rms - min_rms) );

	}//end blip loop
      //}//end particle (all or specific) loop
    }//end event loop

    int list_size = static_cast<int>(blip_num.size());
    cout<<"\nMax hits: "<<max_hits<<"\n"<<   // "\nNum hits: "<<list_size<<"\n\n";
	  "Num protons: "<<totalPrs<<"\n"<<
	  "Total sample blips: "<<blipN<<"\n"<<
	  "Blips with more than 2 hits: "<<over2<<"\n"
	  "Blips with more than 4 hits: "<<over4<<"\n";
    // loop through the hits
    for( int i=0; i<list_size; i++ ) {
      bool enough_wires_plane = false; auto wire_map = blip_wire_per_plane_map[blip_num[i]];
      for (const auto& plane_wires : wire_map) {
	if (plane_wires.second.size() >= 2)
	  { enough_wires_plane = true; break; }
      }   // only include blips with >=2 hits and spanning >1 plane
      if( blip_Nhits_map[blip_num[i]] > 2 && std::find(bad_num.begin(), bad_num.end(), blip_num[i]) == bad_num.end() && enough_wires_plane ) {
	h_normPz->Fill(mom_z[i]);
        outfile << blip_num[i] << " " <<
		   rms[i] << " " <<
                   //charge[i] << " " <<
		   amp[i] << " " <<
		   wire[i] << " " <<
		   plane[i] << " " <<
		   time[i] << " " <<
                   //mom[i] << "\n";
		   mom_x[i] << " " <<
		   mom_y[i] << " " <<
		   mom_z[i] << " " <<
		   energy[i]<< "\n";
		   //g4id[i]  << " " <<
		   //event[i] << "\n";
      }
    }

    std::cout << "\nTotal Cosmics: "<<Ncosmics<< std::endl;

    std::vector<std::pair<float, float>> energy_ranges = {{0, 10}, {10, 40}, {40, 50}, {50, 100}};
    std::cout << "\nEnergy calculation:" << std::endl;
    std::cout << "===================" << std::endl;
    std::cout << "Minimum blip energy: "<< minE << std::endl;
    for (const auto& range : energy_ranges) {
        float range_min = range.first;
        float range_max = range.second;
        
        std::map<int, int> blip_hits_count;
        
        // Count hits per blip ID in this energy range
        for (int i = 0; i < list_size; i++) {
          if (blip_Nhits_map[blip_num[i]] > 2 && std::find(bad_num.begin(), bad_num.end(), blip_num[i]) == bad_num.end() && energy[i] >= range_min && energy[i] < range_max) {
            blip_hits_count[blip_num[i]] = blip_Nhits_map[blip_num[i]];
          }
        } 
        if (!blip_hits_count.empty()) {
          float sum_hits = 0;
          for (const auto& pair : blip_hits_count) { sum_hits += pair.second; }
          float avg_hits = sum_hits / blip_hits_count.size();
          std::cout<<range_min<<"-"<<range_max<<" MeV: "<< std::fixed << std::setprecision(2) << avg_hits 
                   <<" avg hits ("<<blip_hits_count.size()<<" blips)\n";
        } else { std::cout<<range_min<<"-"<<range_max<<" MeV: "<<"0.00 avg hits (0 blips)\n"; }
    }std::cout<<"Maximum blip energy: "<<maxE<<std::endl;

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    h_BlipE->Draw("hist");

    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600);
    h_NumBlipHits->Draw("hist");

    TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600);
    h_cos_comparison->Draw("hist");

    TCanvas *c4 = new TCanvas("c4", "Canvas", 800, 600);
    h_normPz->Draw("hist");

    outfile.close();

  }




