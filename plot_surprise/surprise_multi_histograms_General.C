#include "TDirectory.h"
#include "TText.h"
#include "THStack.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TVector.h"
#include <string>
#include <iostream>
#include "TGraph.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLine.h"
#include<stdio.h>
#include "TAttFill.h"
#include "TBox.h"
#include "THStack.h"
#include "TTree.h"
#include <TLatex.h>
#include "TLegend.h"
#include <iomanip>
#include <fstream>
#include <sstream>

void surprise_multi_histograms_General(const char* file1, const char* file2, const char* file3, const char* file4, const char* config_file = "plot_configs.txt") {

    gStyle->SetOptStat(0);
    const int numHists = 8;
  
    // Structure to hold variable set configuration
    struct VariableSet {
      std::string name;
      std::string histTitle;
      std::string tag;
      std::string hYaxisName;
      std::string hXaxisName;
      float SRUx1;
      float SRUx2;
      bool printBox;
      std::vector<std::string> hist_names;
    };
    
    // Read configuration from file
    std::vector<VariableSet> variable_sets;
    std::ifstream config(config_file);
    if (!config.is_open()) {
      std::cerr << "Error: Could not open configuration file: " << config_file << std::endl;
      return;
    }
    
    std::string line;
    while (std::getline(config, line)) {
      // Skip empty lines and comments
      if (line.empty() || line[0] == '#') continue;
      
      std::istringstream iss(line);
      std::vector<std::string> tokens;
      std::string token;
      
      // Split line into tokens
      while (iss >> token) {
          tokens.push_back(token);
      }
      
      // Check if we have enough tokens (8 basic fields + 8 hist names = 16)
      if (tokens.size() < 16) {
          std::cerr << "Error: Expected 16 tokens but got " << tokens.size() << " in config line: " << line << std::endl;
          continue;
      }
      
      VariableSet vs;
      int idx = 0;
      
      // Read basic configuration
      vs.name = tokens[idx++];
      vs.histTitle = tokens[idx++];
      vs.tag = tokens[idx++];
      vs.hYaxisName = tokens[idx++];
      vs.hXaxisName = tokens[idx++];
      vs.SRUx1 = std::stof(tokens[idx++]);
      vs.SRUx2 = std::stof(tokens[idx++]);
      vs.printBox = (tokens[idx++] == "true");
      
      // Read histogram names
      for (int i = 0; i < numHists; i++) { vs.hist_names.push_back(tokens[idx++]); }

      variable_sets.push_back(vs);

      std::cout << "Successfully parsed variable set: " << vs.name << " with " << vs.hist_names.size() << " histograms" << std::endl;
    }
    config.close();

    if (variable_sets.empty()) {
      std::cerr << "Error: No valid variable sets found in configuration file." << std::endl;
      return;
    }

    const int NUM_VARIABLE_SETS = variable_sets.size();

    // Define colors array (same for all variable sets)
    Int_t colors[numHists] = {
        kBlack,      // h1 - line only
        kGreen-2,    // h2
        kOrange+2,   // h3  
        kPink-5,     // h4
        kBlue-5,     // h5
        kYellow-5,   // h6
        //kRed-5,      // h7
        //kMagenta-5,  // h8
	kAzure-5,    // h9
	kViolet-5    // h10
    };

    // Open or create output ROOT file (recreate mode to include all plots)
    TFile* outFile = new TFile("plots_output.root", "recreate");
    TDirectory* linear_plots           = outFile->mkdir("linear_plots");
    TDirectory* log_plots              = outFile->mkdir("log_plots");
    TDirectory* overlay_plots          = outFile->mkdir("overlay_plots");
    TDirectory* overlay_linear_plots          = outFile->mkdir("overlay_linear_plots");

    // Open the ROOT files
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    TFile *f3 = TFile::Open(file3);
    TFile *f4 = TFile::Open(file4);

    if (!f1 || !f2 || !f3 || !f4) {
        std::cerr << "Error: Could not open one of the files." << std::endl;
        outFile->Close();
        return;
    } 

    // Get total events from tree (do this once)
    //TTree* treedata_f1 = dynamic_cast<TTree*>(f1->Get("evd_tree"));
    //TBranch* eventNumberBranch_f1 = treedata_f1->GetBranch("totalEvents");
    //float TotalEventDisplays_f1;
    //eventNumberBranch_f1->SetAddress(&TotalEventDisplays_f1);
    //treedata_f1->GetEntry(0);
    //cout<<"==============================================================="<<endl;
    //std::cout << "Total Event Displays from blipana outEVD f1: " << TotalEventDisplays_f1 << std::endl;
    //cout<<"==============================================================="<<endl;

    // Loop through all variable sets
    for(int var_set = 0; var_set < NUM_VARIABLE_SETS; var_set++) {
        
        std::cout << "\n=== Processing variable set: " << variable_sets[var_set].name << " ===" << std::endl;
        
        // Declare histogram array for this variable set
        TH1F *hist_var[numHists];
        
        // Retrieve histograms using loop
        bool all_histograms_found = true;
	for(int i = 0; i < numHists; i++) {
	    std::string hist_name = variable_sets[var_set].hist_names[i];
	    std::string dir_path = variable_sets[var_set].name + "_plots";
	    std::string full_path = dir_path + "/" + hist_name;

	    hist_var[i] = (TH1F*)f1->Get(full_path.c_str());
	    if (!hist_var[i]) {
	        // Try as TH1D if TH1F fails
	        TH1D* temp_hist = (TH1D*)f1->Get(full_path.c_str());
	        if (temp_hist) {
	            hist_var[i] = (TH1F*)temp_hist;
	        } else {
	            std::cerr << "Warning: Could not retrieve histogram " << hist_name 
	                     << " from file. Skipping variable set " << variable_sets[var_set].name << std::endl;
	            all_histograms_found = false;
	            break;
	        }
	    }
            std::cout << "Successfully loaded histogram: " << full_path << std::endl;
        }
        
        if (!all_histograms_found) {
            continue; // Skip this variable set
        }

        // Set histogram properties using loop
        for(int i = 0; i < numHists; i++) {
            hist_var[i]->SetTitle(Form("%s", variable_sets[var_set].histTitle.c_str()));
            hist_var[i]->GetXaxis()->SetTitle(Form("%s", variable_sets[var_set].hXaxisName.c_str()));
            hist_var[i]->GetYaxis()->SetTitle(Form("%s", variable_sets[var_set].hYaxisName.c_str()));
            
            // Normalize histograms
            //if (hist_var[i]->GetEntries() > 0) {
            //    hist_var[i]->Scale(1.0 / TotalEventDisplays_f1);
            //}
            
            // Set colors and styles
            hist_var[i]->SetLineColor(colors[i]);
            hist_var[i]->SetLineWidth(3);
            
            // Special styling for first histogram (line only)
            if(i == 0) {
                hist_var[i]->SetLineStyle(kDashed);
                hist_var[i]->SetLineWidth(4);
                hist_var[i]->SetMarkerColor(colors[i]);
            } else {
                // Fill color for stacked histograms
                hist_var[i]->SetFillColor(colors[i]);
            }
        }

        // Create canvas and draw histograms
        TCanvas *c = new TCanvas(Form("canvas_%s_%s", variable_sets[var_set].name.c_str(), variable_sets[var_set].tag.c_str()), 
                               "Stacked Histograms", 1200, 900);

        // Create stack (skip first histogram since it's the line overlay)
        THStack *stack = new THStack("stack", Form("%s;%s;%s", variable_sets[var_set].histTitle.c_str(), 
                                                variable_sets[var_set].hXaxisName.c_str(), variable_sets[var_set].hYaxisName.c_str()));

        // Add histograms to stack in specific order (indices 2,3,7,5,4,1,6)
        //int stack_order[7] = {1, 2, 3, 4, 5, 6, 7};
        for(int i = 1; i < numHists; i++) {
          stack->Add(hist_var[i]);
        }

	double maxValue = hist_var[0]->GetMaximum();
	stack->SetMaximum(maxValue * 1.1);
        
        stack->Draw("HIST");
        hist_var[0]->Draw("same HIST"); // Draw the line histogram on top

        // Get percentage makeup
	int Nbins = hist_var[0]->GetNbinsX();
	double all = hist_var[0]->Integral(0,Nbins+1); double integrals[numHists];
        for(int i = 0; i < numHists; i++) {
	  Nbins = hist_var[i]->GetNbinsX();
	  integrals[i] = hist_var[i]->Integral(0,Nbins+1);
	}

	double sigTrue = hist_var[1]->GetBinContent(2) + hist_var[2]->GetBinContent(2);
	double sigFalse = hist_var[1]->GetBinContent(1) + hist_var[2]->GetBinContent(1);

	double bkgTrue = hist_var[3]->GetBinContent(2) + hist_var[4]->GetBinContent(2) + hist_var[5]->GetBinContent(2) + hist_var[6]->GetBinContent(2) + hist_var[7]->GetBinContent(2);
	double bkgFalse = hist_var[3]->GetBinContent(1) + hist_var[4]->GetBinContent(1) + hist_var[5]->GetBinContent(1) + hist_var[6]->GetBinContent(1) + hist_var[7]->GetBinContent(1);

	std::cout<<"True signal: "<<sigTrue<<" False signal: "<<sigFalse<<"\nTrue background: "<<bkgTrue<<" False background: "<<bkgFalse<<std::endl;

        std::cout<<"signal: "<<integrals[2]<<", non-prim mom: "<<integrals[3]<<std::endl;

        // Add legend
        TLegend *legend = new TLegend(0.62, 0.65, 0.87, 0.89);
        legend->SetFillColorAlpha(0,0.0);
        legend->SetLineColorAlpha(0,0.0);
        
        // Add legend entries with hardcoded labels
        const char* legend_labels[numHists] = {
            "all blips",
            "(n,1p) primary neutrons",
	    "(n, Np) primary neutrons",
            "(n,p) non-primary neutrons",
            "primary protons",
            "photons, e+, and e-",
            "other interactions",
	    //"shower products",
            "cosmics & radiogenics"
        };

        legend->AddEntry(hist_var[0], TString::Format("%s   (%.0f)", legend_labels[0], all), "l");
        for(int i=1; i<numHists; i++) { legend->AddEntry(hist_var[i], TString::Format("%s = %.1f%%", legend_labels[i], (integrals[i]/all)*100), "f"); } //legend->AddEntry(hist_var[i], TString::Format("%s = %.3f%%   (%.0f)", legend_labels[i], (integrals[i]/all)*100, integrals[i]), "f"); }    //, TString::Format("%s", legend_labels[i]), "f"); }       
        legend->Draw();

        c->Update(); c->Draw();
        
        // Write to output ROOT file
        outFile->cd(); linear_plots->cd();
        
        // Save the overlayed linear plot
        c->SetName(Form("linear_%s_%s", variable_sets[var_set].name.c_str(), variable_sets[var_set].tag.c_str()));
        c->SetTitle(Form("Linear Scale: %s", variable_sets[var_set].histTitle.c_str()));
        c->Write();
        
        // new canvas for log scale
        TCanvas *c_log = new TCanvas(Form("log_%s_%s", variable_sets[var_set].name.c_str(), variable_sets[var_set].tag.c_str()), 
                                   Form("Log Scale: %s", variable_sets[var_set].histTitle.c_str()), 1200, 900); 
        
        // new stack for log plot
        THStack *stack_log = new THStack("stack_log", Form("%s;%s;%s", variable_sets[var_set].histTitle.c_str(), 
                                                          variable_sets[var_set].hXaxisName.c_str(), variable_sets[var_set].hYaxisName.c_str()));
        for(int i = 1; i < numHists; i++) { stack_log->Add(hist_var[i]); }

        stack_log->SetMinimum(0.1);  // or your desired minimum
	stack_log->SetMaximum(maxValue * 2);  // 2x headroom for log scale

	std::string histName = variable_sets[var_set].name;
	if( histName == "recoVtxDist" ) { stack_log->SetMinimum(1.0); } else if( histName == "cosAlignment" ) { stack_log->SetMinimum(1.0); }
        else if( histName == "alignmentAngle" ) { stack_log->SetMinimum(1.0); } else if( histName == "resolution" ) { stack_log->SetMinimum(1.0); }
	else if( histName == "byDeadWire" ) { stack_log->SetMinimum(10.0); stack_log->SetMaximum(pow(10,7)); } //stack_log->GetXaxis()->SetNdivisions(2, kFALSE); }
	else if( histName == "dw" ) { stack_log->SetMinimum(10); } else if( histName == "size" ) { stack_log->SetMinimum(3); }
	else if( histName == "inShower" ) { stack_log->SetMaximum(pow(10,5)); } //stack_log->GetXaxis()->SetNdivisions(2, kFALSE); }

	stack_log->Draw("HIST"); hist_var[0]->Draw("same HIST");
        legend->Draw();
        
        // log scale on y-axis
        c_log->SetLogy(1); c_log->Update();
        
        // Write log version to file
        log_plots->cd(); c_log->Write();

        // new canvas for log overlay
        TCanvas *c_over = new TCanvas(Form("overlay_%s_%s", variable_sets[var_set].name.c_str(), variable_sets[var_set].tag.c_str()),
					     "Overlayed Histograms", 1200, 900);

        hist_var[1]->SetMinimum(0.1);
	hist_var[1]->SetMaximum(maxValue * 2);  // 2x headroom for log scale

	int styles[] = {1,2,9,10,7,1,2,9};
	for (int i = 1; i < numHists; i++) {

	  hist_var[i]->SetFillStyle(0);
	  hist_var[i]->SetLineWidth(6);

	  hist_var[i]->SetLineStyle(styles[i-1]);

	  hist_var[i]->Draw("same HIST");
	}

	legend->Draw();
        c_over->SetLogy(1); c_over->Update();
	overlay_plots->cd(); c_over->Write();

        // new canvas for overlay
        TCanvas *c_overLin = new TCanvas(Form("overlay_linear_%s_%s", variable_sets[var_set].name.c_str(), variable_sets[var_set].tag.c_str()),
                                             "Overlayed Linear Histograms", 1200, 900);

        //hist_var[1]->SetMinimum(0.1);
        //hist_var[1]->SetMaximum(maxValue * 2);  // 2x headroom for log scale

        for (int i = 1; i < numHists; i++) {

          hist_var[i]->SetFillStyle(0);
          hist_var[i]->SetLineWidth(6);

          hist_var[i]->SetLineStyle(styles[i-1]);

          hist_var[i]->Draw("same HIST");
        }

        legend->Draw();
        c_overLin->Update();
        overlay_linear_plots->cd(); c_overLin->Write();


        // Clean up canvases for this iteration
        delete c;
        delete c_log;
	delete c_over;
	delete c_overLin;

        std::cout << "Completed processing: " << variable_sets[var_set].name << std::endl;
    }
    
    // Close and save the output file
    outFile->Close();
    std::cout << "\n=== All plots saved to plots_output.root ===" << std::endl;
    
    // Clean up files
    f1->Close();
    f2->Close(); 
    f3->Close();
    f4->Close();
    delete outFile;
}
