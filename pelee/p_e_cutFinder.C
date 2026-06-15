
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

  int			     nblips_saved = 0;
  std::vector<float>       * blip_energy = 0;
  std::vector<float>       * blip_dx = 0;
  std::vector<float>       * blip_dw = 0;
  std::vector<float>       * blip_true_pdg = 0;
  std::vector<float>       * blip_true_g4id = 0;

  //std::string fFileName = "checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root";
  std::string fFileName = "MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260224.root";

  TH2D*   h_elec_E_ds;
  TH2D*   h_pr_E_ds;

  void p_e_cutFinder(){

    //float min = std::min_element(blip_density.begin(), blip_density.end());
    //float max = std::max_element(blip_density.begin(), blip_density.end());
    //std::sort(blip_density.begin(), blip_density.end());
    //int guess = blip_density[size/2];

    // Read in the TFile and grab the TTree
    TFile* file = new TFile(("files/"+fFileName).c_str(),"READ");
    TTree* tree = (TTree*)file->Get("nuselection/NeutrinoSelectionFilter");
    std::cout<<"Reading data file "<<fFileName<<"\n";

    tree->SetBranchAddress("nblips_saved",       &nblips_saved);
    tree->SetBranchAddress("blip_energy",	 &blip_energy);
    tree->SetBranchAddress("blip_dx",		 &blip_dx);
    tree->SetBranchAddress("blip_dw",            &blip_dw);
    tree->SetBranchAddress("blip_true_pdg",	 &blip_true_pdg);
    tree->SetBranchAddress("blip_true_g4id",	 &blip_true_g4id);

    // make output file to store plots
    TFile* outfile = new TFile("fit_cut.root", "recreate");

    h_elec_E_ds = new TH2D("elec_E_ds","Electron Blip Energy Deposition Density;Reco blip E [MeVee]; Eblip/ds [MeV/cm]",500,0,10,400,0,8);
    h_pr_E_ds   = new TH2D("pr_E_ds","Proton Blip Energy Deposition Density;Reco blip E [MeVee]; Eblip/ds [MeV/cm]",500,0,10,400,0,8);

    int p; std::map<int, int> map_g4id_index; bool protonMode=false;
    float blipE; float density;
    std::map<float, float> map_E_den; std::vector<float> E_bins;
    float mean, std; int count; float z_score = 1.28;   // 90th percentile for normal distro
    float cutoff; std::vector<float> densities;
    std::vector<float> cutoffs;

    int numEvents = tree->GetEntries(); //get_branches(tree, false);
    //============================== EVENT LOOP ==============================//
    for( int iEvent = 0; iEvent < numEvents; iEvent++ ) {
      tree->GetEntry(iEvent);

      for( int i=0; i<nblips_saved; i++ ) {
	if( blip_true_g4id->at(i) < 0 ) continue;
	blipE = blip_energy->at(i);
	float dx = blip_dx->at(i); float dw = blip_dw->at(i);
	float blip_size = sqrt( dx*dx + dw*dw );
	density = blipE/blip_size;
	int blip_pdg = blip_true_pdg->at(i);
	if( blip_pdg == 11 ) {
	  h_elec_E_ds->Fill(blipE, density);
	  map_E_den[blipE] = density;
	}
	else if( blip_pdg == 2212 ) { h_pr_E_ds->Fill(blipE, density); }
      } // end blip loop

      map_g4id_index.clear();
    } // end event loop 

    int map_size=map_E_den.size(); float bin_min=3.0; float bin_max=3.25; float step_size=0.25;

    while( bin_max <=10 ){
      densities.clear();
      for (const auto& [blipE, density] : map_E_den) {
	if ( blipE > bin_min && blipE <= bin_max && std::isfinite(density) ) {
	  densities.push_back(density);
	}
      }
      //mean = TMath::Mean(densities.begin(), densities.end());
      //std = TMath::StdDev(densities.begin(), densities.end());
      // percentile point = mean + (z-score for percentile) * (standard deviation) = cutoff = mean + ( z_score * std ); (NORMAL DISTRO)
      //std::cout<<"\nBin: "<<bin_min<<"\ncount = "<<count<<", mean = "<<mean<<", std = "<<std<<", cutoff = "<<cutoff<<std::endl;
      std::sort(densities.begin(), densities.end());
      count = densities.size(); cutoff = densities[0.99*(count-1)];
      cutoffs.push_back(cutoff);

      bin_min = bin_max; bin_max += step_size;

    }

    for(float E = 3.125; E <= 9.875; E += 0.25) { E_bins.push_back(E); }
    TGraph* gr = new TGraph(E_bins.size(), &E_bins[0], &cutoffs[0]);
    gr->SetMarkerStyle(20); gr->SetMarkerSize(1.1);

    gr->SetMarkerColor(kAzure);
    TCanvas *c1 = new TCanvas("c1", "Canvas", 1100, 800);
    //gStyle->SetPalette();
    h_elec_E_ds->DrawCopy("hist");
    gr->Draw("PL SAME");

    TCanvas *c2 = new TCanvas("c2", "Canvas", 1100, 800);
    gStyle->SetPalette(74);
    h_pr_E_ds->DrawCopy("hist");
    gr->Draw("PL SAME");

    TF1 *E_func = new TF1("E_func", "[0] * log( [1] * x ) + [2]*x + [3]", 3, 10);   // x range 
    // Set initial parameter estimates
    E_func->SetParameters(0.1, 0.01, -0.3, 4.0); E_func->SetLineColor(kAzure);//kRed);
    E_func->SetLineWidth(6); gr->Fit(E_func, "", "", 3, 10);
    E_func->Write("cutoff_graph");
    TTree* param_tree = new TTree("param_tree", "Fit Parameters"); 
    float A = E_func->GetParameter(0);
    float B = E_func->GetParameter(1);
    float C = E_func->GetParameter(2);
    float D = E_func->GetParameter(3);
    param_tree->Branch("Param_A", &A, "Param_A/F");
    param_tree->Branch("Param_B", &B, "Param_B/F");
    param_tree->Branch("Param_C", &C, "Param_C/F");
    param_tree->Branch("Param_D", &D, "Param_D/F");
    param_tree->Fill();

    // Close the histogram file  
    outfile->Write();
    outfile->Close();

  }

