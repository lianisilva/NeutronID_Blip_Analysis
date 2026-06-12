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
#include <iomanip>      // std::setprecision
#include <fstream>

    void multi_histograms_General(const char* file1, const char* file2, const char* file3, const char* file4,const char* histTitle,const char* histName, const char* tag, const char* hYaxisName, const char* hXaxisName, float SRUx1, float SRUx2, bool printBox) {
	gStyle->SetOptStat(0);
    
	// Open or create output ROOT file (append mode)
	TFile* outFile = new TFile("all_plots_output.root", "recreate");

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

        // Retrieve histograms
        //TH1F *h1 = (TH1F*)f1->Get("blipE");
        //TH1F *h2 = (TH1F*)f1->Get("primPrBlipE");
        //TH1F *h3 = (TH1F*)f1->Get("prBlipE");
        //TH1F *h4 = (TH1F*)f1->Get("prNonPrimBlipE");
        //TH1F *h5 = (TH1F*)f1->Get("gammaBlipE");
        //TH1F *h6 = (TH1F*)f1->Get("cosmicBlipE");
        //TH1F *h7 = (TH1F*)f1->Get("neutronBlipE");
        //TH1F *h8 = (TH1F*)f1->Get("nonNeutBlipE");

        //TH1F *h1 = (TH1F*)f1->Get("trueBlipE");
        //TH1F *h2 = (TH1F*)f1->Get("neutronTrueE");
        //TH1F *h5 = (TH1F*)f1->Get("gammaTrueE");
        //TH1F *h3 = (TH1F*)f1->Get("prTrueE");
        //TH1F *h4 = (TH1F*)f1->Get("allPrTrueE");
        //TH1F *h6 = (TH1F*)f1->Get("trueBlipE");

        TH1F *h1 = (TH1F*)f1->Get("nblips_vs_vtxDist");
        TH1F *h2 = (TH1F*)f1->Get("primPrBlips_vs_vtxDist"); 
        TH1F *h3 = (TH1F*)f1->Get("prBlips_vs_vtxDist");
	TH1F *h4 = (TH1F*)f1->Get("prNonPrimBlips_vs_vtxDist");
        TH1F *h5 = (TH1F*)f1->Get("gamBlips_vs_vtxDist");
        TH1F *h6 = (TH1F*)f1->Get("cosmicBlips_vs_vtxDist");
        TH1F *h7 = (TH1F*)f1->Get("N_blips_vs_vtxDist");
        TH1F *h8 = (TH1F*)f1->Get("nonN_blips_vs_vtxDist");
 
        //TH1F *h1 = (TH1F*)f1->Get("cos_all");
        //TH1F *h2 = (TH1F*)f1->Get("cos_primPr");
        //TH1F *h3 = (TH1F*)f1->Get("cos_pr");
        //TH1F *h4 = (TH1F*)f1->Get("cos_prNonPrim");
        //TH1F *h5 = (TH1F*)f1->Get("cos_gamma");
        //TH1F *h6 = (TH1F*)f1->Get("cos_cosmic");
        //TH1F *h7 = (TH1F*)f1->Get("cos_neutron");
        //TH1F *h8 = (TH1F*)f1->Get("cos_other");

        //TH1F *h1 = (TH1F*)f1->Get("allRes");
        //TH1F *h2 = (TH1F*)f1->Get("primPrRes");
        //TH1F *h3 = (TH1F*)f1->Get("prRes");
        //TH1F *h4 = (TH1F*)f1->Get("prNonPrimRes");
	//TH1F *h5 = (TH1F*)f1->Get("gammaRes");
        //TH1F *h6 = (TH1F*)f1->Get("cosmicRes");
        //TH1F *h7 = (TH1F*)f1->Get("neutronRes");
        //TH1F *h8 = (TH1F*)f1->Get("nonNeutRes");

        //TH1F *h1 = (TH1F*)f1->Get("allDensity");
        //TH1F *h2 = (TH1F*)f1->Get("primPrDensity");
        //TH1F *h3 = (TH1F*)f1->Get("prDensity");
        //TH1F *h4 = (TH1F*)f1->Get("prNonPrimDensity");
        //TH1F *h5 = (TH1F*)f1->Get("gammaDensity");
        //TH1F *h6 = (TH1F*)f1->Get("cosmicDensity");
        //TH1F *h7 = (TH1F*)f1->Get("neutronDensity");
        //TH1F *h8 = (TH1F*)f1->Get("nonNeutDensity");

        if (!h1 || !h2 || !h3 || !h4 || !h5 || !h6 || !h7 || !h8) {
            std::cerr << "Error: Could not retrieve histogram from one of the files." << std::endl;
            outFile->Close();
            return;
        }

	//h1->Sumw2();
	//h2->Sumw2();
	//h3->Sumw2();
	//h4->Sumw2();
	//h5->Sumw2();
	//h6->Sumw2();
        //h7->Sumw2();
        //h8->Sumw2();

	TTree* treedata_f1 = dynamic_cast<TTree*>(f1->Get("evd_tree"));
	TBranch* eventNumberBranch_f1 = treedata_f1->GetBranch("totalEvents");   //total_blips
	float TotalEventDisplays_f1;
	eventNumberBranch_f1->SetAddress(&TotalEventDisplays_f1);
	treedata_f1->GetEntry(0);
	cout<<"==============================================================="<<endl;
	std::cout << "Total Event Displays from blipana outEVD f1: " << TotalEventDisplays_f1 << std::endl;
	cout<<"==============================================================="<<endl;

        h1->SetTitle(Form("%s", histTitle));
	h1->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
	h1->GetYaxis()->SetTitle(Form("%s", hYaxisName ));
	h2->SetTitle(Form("%s", histTitle));
        h2->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
        h2->GetYaxis()->SetTitle(Form("%s", hYaxisName ));
        h3->SetTitle(Form("%s", histTitle));
	h3->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
        h3->GetYaxis()->SetTitle(Form("%s", hYaxisName ));
        h4->SetTitle(Form("%s", histTitle));
        h4->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
        h4->GetYaxis()->SetTitle(Form("%s", hYaxisName ));
        h5->SetTitle(Form("%s", histTitle));
        h5->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
        h5->GetYaxis()->SetTitle(Form("%s", hYaxisName ));
        h6->SetTitle(Form("%s", histTitle));
        h6->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
        h6->GetYaxis()->SetTitle(Form("%s", hYaxisName ));
        h7->SetTitle(Form("%s", histTitle));
        h7->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
        h7->GetYaxis()->SetTitle(Form("%s", hYaxisName ));
        h8->SetTitle(Form("%s", histTitle));
        h8->GetXaxis()->SetTitle(Form("%s", hXaxisName ));
        h8->GetYaxis()->SetTitle(Form("%s", hYaxisName ));

        // Normalize histograms by the number of entries
        if (h1->GetEntries() > 0) h1->Scale(1.0 / TotalEventDisplays_f1);
        if (h2->GetEntries() > 0) h2->Scale(1.0 / TotalEventDisplays_f1); // normalizing by f1 too 
        if (h3->GetEntries() > 0) h3->Scale(1.0 / TotalEventDisplays_f1);
	if (h4->GetEntries() > 0) h4->Scale(1.0 / TotalEventDisplays_f1);
        if (h5->GetEntries() > 0) h5->Scale(1.0 / TotalEventDisplays_f1);
        if (h6->GetEntries() > 0) h6->Scale(1.0 / TotalEventDisplays_f1);
        if (h7->GetEntries() > 0) h7->Scale(1.0 / TotalEventDisplays_f1);
        if (h8->GetEntries() > 0) h8->Scale(1.0 / TotalEventDisplays_f1);

	// Set X-axis range if specified
	//if(SRUx1 != SRUx2) {
	//  if(SRUx1 != SRUx2) h1->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//  if(SRUx1 != SRUx2) h2->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//  if(SRUx1 != SRUx2) h3->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//  if(SRUx1 != SRUx2) h4->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//  if(SRUx1 != SRUx2) h5->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//  if(SRUx1 != SRUx2) h6->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//}

        // Create canvas and draw histograms
        TCanvas *c = new TCanvas(Form("canvas_%s_%s_%s", histTitle, histName, tag), "Stacked Histograms", 1200, 900);

	h1->SetLineColor(kBlack);
	h1->SetMarkerColor(kBlack);
	h1->SetLineStyle(kDashed);
	h1->SetLineWidth(4);

	h2->SetLineColor(kGreen-2);
	h2->SetFillColor(kGreen-2);
        h2->SetLineWidth(3);

	h3->SetLineColor(kOrange+2);
	h3->SetFillColor(kOrange+2);
	//h3->SetFillStyle(3005); // Cross-hatch fill pattern
        h3->SetLineWidth(3);

        h4->SetLineColor(kPink-5);
        h4->SetFillColor(kPink-5);
        h4->SetLineWidth(3);

        h5->SetLineColor(kBlue-5);
        h5->SetFillColor(kBlue-5);
        h5->SetLineWidth(3);

        h6->SetLineColor(kYellow-5);
        h6->SetFillColor(kYellow-5) ;
        h6->SetLineWidth(3);

        h7->SetLineColor(kRed-5);
        h7->SetFillColor(kRed-5) ;
        h7->SetLineWidth(3);

        h8->SetLineColor(kMagenta-5);
        h8->SetFillColor(kMagenta-5) ;
        h8->SetLineWidth(3);

        //h1->Draw("HIST");//h5->Draw("HIST SAME");//h4->Draw("HIST SAME");//h3->Draw("HIST SAME");

        THStack *stack = new THStack("stack", Form("%s;%s;%s", histTitle, hXaxisName, hYaxisName));
	stack->Add(h3);
 	stack->Add(h4);
        stack->Add(h8);
        stack->Add(h6);
        stack->Add(h5);
        stack->Add(h2);
        stack->Add(h7);
	stack->Draw("HIST");
	h1->Draw("same HIST");

	// Get percentage makeup
	double all = h1->Integral();   // Return the array of bins integral.
	double prims = h2->Integral();
        double momPrim = h3->Integral();
        double momNoPrim = h4->Integral();
        double gams = h5->Integral();
        double cosmics = h6->Integral();
	double otherN = h7->Integral();
        double nonN = h8->Integral();

	std::cout<<"\nsignal: "<<momPrim<<"\nnon-prim mom: "<<momNoPrim<<"\n";

        // Add legend
        TLegend *legend = new TLegend(0.62, 0.65, 0.87, 0.89);   // Xlow, Ylow, Xup, Yup
        legend->SetFillColorAlpha(0,0.0);
        legend->SetLineColorAlpha(0,0.0);
 	legend->AddEntry(h1,"all blips in event", "l");
 	legend->AddEntry(h3, TString::Format("(n,p) - primary neutrons = %.1f%%", (momPrim/all)*100 ), "f");
 	legend->AddEntry(h4, TString::Format("(n,p) - NON-primary neutrons = %.1f%%", (momNoPrim/all)*100 ), "f");
	legend->AddEntry(h2, TString::Format("primary protons = %.1f%%", (prims/all)*100 ), "f");
	legend->AddEntry(h5, TString::Format("(n,gamma) - all neutrons = %.1f%%", (gams/all)*100 ), "f");
        legend->AddEntry(h7, TString::Format("other neutrons = %.1f%%", (otherN/all)*100 ), "f");
        legend->AddEntry(h8, TString::Format("other interactions = %.1f%%", (nonN/all)*100 ), "f");
        legend->AddEntry(h6, TString::Format("cosmics = %.1f%%", (cosmics/all)*100 ), "f");
        legend->Draw();

	// Update canvas
	c->Update();
	c->Draw();
	
	// Write to output ROOT file
	outFile->cd();
	
	// Save the overlayed linear plot
	c->SetName(Form("canvas_linear_%s_%s_%s", histTitle, histName, tag));
	c->SetTitle(Form("Linear Scale: %s", histTitle));
	c->Write();
	
	// Create a new canvas for log scale to preserve the linear one
	TCanvas *c_log = new TCanvas(Form("canvas_log_%s_%s_%s", histTitle, histName, tag), 
	                           Form("Log Scale: %s", histTitle), 1200, 900); 
	
	// Create new stack for log plot
	THStack *stack_log = new THStack("stack_log", Form("%s;%s;%s", histTitle, hXaxisName, hYaxisName));
	stack_log->Add(h3);
	stack_log->Add(h4);
	stack_log->Add(h2);
	stack_log->Add(h5);
	stack_log->Add(h7);
	stack_log->Add(h8);
	stack_log->Add(h6);
	stack_log->Draw("HIST");
	//stack_log->GetYaxis()->SetRangeUser(0.1, 1e3.2);//1e-4, 0.04);
	stack_log->SetMinimum(TMath::Power(10, -4.3));
	stack_log->SetMaximum(0.04);
	h1->Draw("same HIST");
	
	legend->Draw();
	
	// Set log scale on y-axis
	c_log->SetLogy(1); 
	c_log->Update();
	
	// Write log version to file
	c_log->Write();
	
	
	// Close and save the output file
	outFile->Close();
	
	// Clean up
	delete outFile;
      }	
