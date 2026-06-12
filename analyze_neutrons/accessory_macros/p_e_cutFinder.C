
  #include "set_branches.cpp"
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
  #include "TArrow.h"

  int maxEvents = -2000;

  std::string fFileName1 = "../blipana_eminus_0to10MeV_overlay_lowthresh_20250731.root";
  std::string fFileName2 = "../blipana_protons_0to200MeV_overlay_20250826.root";

  //std::string fFileName1 = "../blipana_prodgenie_bnb_nu_overlay_mcc9_run3_numuCCinclusive_selection_20251028";

  TH2D*   h_elec_E_ds;
  TH2D*   h_pr_E_ds;

  void p_e_cutFinder(){

    //float min = std::min_element(blip_density.begin(), blip_density.end());
    //float max = std::max_element(blip_density.begin(), blip_density.end());
    //std::sort(blip_density.begin(), blip_density.end());
    //int guess = blip_density[size/2];

    // Read in the TFile and grab the TTree 
    TFile* file1 = new TFile(fFileName1.c_str(),"READ");
    TFile* file2 = new TFile(fFileName2.c_str(),"READ");
    TTree* elecTree = (TTree*)file1->Get("blipana/anatree");
    TTree* prTree = (TTree*)file2->Get("blipana/anatree");
    std::cout<<"Reading data files "<<fFileName1<<" and "<<fFileName2<<"\n";
    // make output file to store plots
    TFile* outfile = new TFile("fit_cut.root", "recreate");

    h_elec_E_ds = new TH2D("elec_E_ds","Electron Blip Energy Deposition Density;Reco blip E [MeVee]; Eblip/ds [MeV/cm]",500,0,10,250,0,5);
    h_pr_E_ds   = new TH2D("pr_E_ds","Proton Blip Energy Deposition Density;Reco blip E [MeVee]; Eblip/ds [MeV/cm]",500,0,10,400,0,8);

    int p; std::map<int, int> map_g4id_index; bool protonMode=false;
    float blipE; float density;
    std::map<float, float> map_E_den; std::vector<float> lowE_bins, highE_bins;
    float mean, std; int count; float z_score = 1.28;   // 90th percentile for normal distro
    float lowCutoff, highCutoff; std::vector<float> densities;
    std::vector<float> lowCutoffs, highCutoffs;

    int numEvents = elecTree->GetEntries(); get_branches(elecTree, false);
    //============================== ELECTRON EVENT LOOP ==============================//
    for( int iEvent = 0; iEvent < numEvents; iEvent++ ) {
      elecTree->GetEntry(iEvent);

      for( int i=0; i<nparticles; i++ ) { map_g4id_index[part_g4id[i]] = i; }

      for( int i=0; i<nblips; i++ ) {
	if( blip_edepid[i] < 0 ) continue;
	blipE = blip_energy[i]; density = blipE/blip_size[i];
	p = map_g4id_index[blip_g4id[i]];
	if( part_isPrimary[p] && part_pdg[p] == 11 ) { h_elec_E_ds->Fill(blipE, density); map_E_den[blipE] = density; }
      } // end blip loop

      map_g4id_index.clear();
    } // end event loop 

    int map_size=map_E_den.size(); float bin_min=3.0; float bin_max=3.5; float step_size=0.5;

    while( bin_max <=10 ){
      densities.clear();
      for (const auto& [blipE, density] : map_E_den) {
	if ( blipE > bin_min && blipE <= bin_max && std::isfinite(density) ) {
	  densities.push_back(density);
	}
      }
      mean = TMath::Mean(densities.begin(), densities.end());
      std = TMath::StdDev(densities.begin(), densities.end());
      // percentile point = mean + (z-score for percentile) * (standard deviation) = highCutoff = mean + ( z_score * std ); (NORMAL DISTRO)
      //std::cout<<"\nBin: "<<bin_min<<"\ncount = "<<count<<", mean = "<<mean<<", std = "<<std<<", cutoff = "<<cutoff<<std::endl;
      std::sort(densities.begin(), densities.end());
      count = densities.size(); highCutoff = densities[0.99*(count-1)];
      highCutoffs.push_back(highCutoff);

      bin_min = bin_max; bin_max += step_size;

    }
    map_E_den.clear();


    numEvents = prTree->GetEntries(); get_branches(prTree, false);
    //============================== PROTON EVENT LOOP ==============================//
    for( int iEvent = 0; iEvent < numEvents; iEvent++ ) {
      prTree->GetEntry(iEvent);

      for( int i=0; i<nparticles; i++ ) { map_g4id_index[part_g4id[i]] = i; }

      for( int i=0; i<nblips; i++ ) {
	if( blip_edepid[i] < 0 ) continue;
	blipE = blip_energy[i]; density = blipE/blip_size[i];
	p = map_g4id_index[blip_g4id[i]];
	if( part_isPrimary[p] && part_pdg[p] == 2212 ) { h_pr_E_ds->Fill(blipE, density); map_E_den[blipE] = density; }
      } // end blip loop

      map_g4id_index.clear();
    } // end event loop
 
    map_size = map_E_den.size(); bin_min = 0.0; bin_max = 0.5; step_size = 0.5;

    while( bin_max <=3 ){
      densities.clear();
      for (const auto& [blipE, density] : map_E_den) {
	if ( blipE > bin_min && blipE <= bin_max && std::isfinite(density) ) {
	  densities.push_back(density);
	}
      }
      mean = TMath::Mean(densities.begin(), densities.end());
      std = TMath::StdDev(densities.begin(), densities.end());

      lowCutoff = mean - ( z_score * std ); lowCutoffs.push_back(lowCutoff);

      bin_min = bin_max; bin_max += step_size;

    }

    for(float E = 3.25; E <= 9.75; E += 0.5) { highE_bins.push_back(E); }
    for(float E = 0.25; E <= 2.75; E += 0.5) { lowE_bins.push_back(E); }
    TGraph* highGr = new TGraph(highE_bins.size(), &highE_bins[0], &highCutoffs[0]);
    TGraph* lowGr = new TGraph(lowE_bins.size(), &lowE_bins[0], &lowCutoffs[0]);
    highGr->SetMarkerStyle(20); highGr->SetMarkerSize(1.1);
    lowGr->SetMarkerStyle(20); lowGr->SetMarkerStyle(20);

    highGr->SetMarkerColor(kAzure);
    TCanvas *c1 = new TCanvas("c1", "Canvas", 1100, 800);
    //gStyle->SetPalette();
    h_elec_E_ds->DrawCopy("hist");
    highGr->Draw("P SAME");
    lowGr->Draw("P SAME");

    lowGr->SetMarkerColor(kAzure);
    TCanvas *c2 = new TCanvas("c2", "Canvas", 1100, 800);
    gStyle->SetPalette(74);
    h_pr_E_ds->DrawCopy("hist");
    highGr->Draw("P SAME");
    lowGr->Draw("P SAME");

    TF1 *highE_func = new TF1("highE_func", "[0] * log( [1] * x ) + [2]*x + [3]", 3, 10);   // x range 
    TF1 *lowE_func = new TF1("lowE_func", "-[0]*x*x*x + [1]*x*x + [2]*x + [3]/x + [4]", 0.25, 3);   // x range
    // Set initial parameter estimates
    highE_func->SetParameters(1.0, 0.01, -0.01, 6.0); highE_func->SetLineColor(kAzure);//kRed);
    lowE_func->SetParameters(0.01, -0.1, 2.0, 0.1, -1.0); lowE_func->SetLineColor(kAzure);
    highE_func->SetLineWidth(6); lowE_func->SetLineWidth(6);
    highGr->Fit(highE_func, "", "", 3, 10); highGr->Write("highCutoff_graph");
    lowGr->Fit(lowE_func, "", "", 0.25, 3); lowGr->Write("lowCutoff_graph");
    TTree* param_tree = new TTree("param_tree", "Fit Parameters"); 
    float lowA = lowE_func->GetParameter(0); //double error_A = lowE_func->GetParError(0);
    float lowB = lowE_func->GetParameter(1);
    float lowC = lowE_func->GetParameter(2);
    float lowD = lowE_func->GetParameter(3);
    float lowE = lowE_func->GetParameter(4);
    float highA = highE_func->GetParameter(0);
    float highB = highE_func->GetParameter(1);
    float highC = highE_func->GetParameter(2);
    float highD = highE_func->GetParameter(3);
    param_tree->Branch("highParam_A", &highA, "highParam_A/F");
    param_tree->Branch("highParam_B", &highB, "highParam_B/F");
    param_tree->Branch("highParam_C", &highC, "highParam_C/F");
    param_tree->Branch("highParam_D", &highD, "highParam_D/F");
    param_tree->Branch("lowParam_A", &lowA, "lowParam_A/F");
    param_tree->Branch("lowParam_B", &lowB, "lowParam_B/F");
    param_tree->Branch("lowParam_C", &lowC, "lowParam_C/F");
    param_tree->Branch("lowParam_D", &lowD, "lowParam_D/F");
    param_tree->Branch("lowParam_E", &lowE, "lowParam_E/F");
    param_tree->Fill();

    // Close the histogram file  
    outfile->Write();
    outfile->Close();

  }

