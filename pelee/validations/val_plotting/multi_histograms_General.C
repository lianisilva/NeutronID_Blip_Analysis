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
    
	// Open the ROOT files
        TFile *f1 = TFile::Open(file1);
        TFile *f2 = TFile::Open(file2);
        TFile *f3 = TFile::Open(file3);
	TFile *f4 = TFile::Open(file4);
        
        if (!f1 || !f2 || !f3 || !f4) {
            std::cerr << "Error: Could not open one of the files." << std::endl;
            return;
        }
        
        // Retrieve histograms
        TH1F *h1 = (TH1F*)f1->Get(histName);
        TH1F *h2 = (TH1F*)f2->Get(histName);
        TH1F *h3 = (TH1F*)f3->Get(histName);
	TH1F *h4 = (TH1F*)f4->Get(histName);
        
        if (!h1 || !h2 || !h3 || !h4) {
            std::cerr << "Error: Could not retrieve histogram from one of the files." << std::endl;
            return;
        }


	//h1->Sumw2();
	//h2->Sumw2();
	//h3->Sumw2();
	//h4->Sumw2();


	TTree* treedata_f1 = dynamic_cast<TTree*>(f1->Get("evd_tree"));
	TBranch* eventNumberBranch_f1 = treedata_f1->GetBranch("total_event_displays");
	float TotalEventDisplays_f1;
	eventNumberBranch_f1->SetAddress(&TotalEventDisplays_f1);
	treedata_f1->GetEntry(0);
	cout<<"==============================================================="<<endl;
	std::cout << "Total Event Displays from blipana outEVD f1: " << TotalEventDisplays_f1 << std::endl;
	cout<<"==============================================================="<<endl;
	
	TTree* treedata_f2 = dynamic_cast<TTree*>(f2->Get("evd_tree"));
	TBranch* eventNumberBranch_f2 = treedata_f2->GetBranch("total_event_displays");
	float TotalEventDisplays_f2;
	eventNumberBranch_f2->SetAddress(&TotalEventDisplays_f2);
	treedata_f2->GetEntry(0);
	cout<<"==============================================================="<<endl;
	std::cout << "Total Event Displays from blipana outEVD f2: " << TotalEventDisplays_f2 << std::endl;
	cout<<"==============================================================="<<endl;
	
	
	TTree* treedata_f3 = dynamic_cast<TTree*>(f3->Get("evd_tree"));
	TBranch* eventNumberBranch_f3 = treedata_f3->GetBranch("total_event_displays");
	float TotalEventDisplays_f3;
	eventNumberBranch_f3->SetAddress(&TotalEventDisplays_f3);
	treedata_f3->GetEntry(0);
	cout<<"==============================================================="<<endl;
	std::cout << "Total Event Displays from blipana outEVD f3: " << TotalEventDisplays_f3 << std::endl;
	cout<<"==============================================================="<<endl;
	
	TTree* treedata_f4 = dynamic_cast<TTree*>(f4->Get("evd_tree"));
	TBranch* eventNumberBranch_f4 = treedata_f4->GetBranch("total_event_displays");
	float TotalEventDisplays_f4;
	eventNumberBranch_f4->SetAddress(&TotalEventDisplays_f4);
	treedata_f4->GetEntry(0);
	cout<<"==============================================================="<<endl;
	std::cout << "Total Event Displays from blipana outEVD f4: " << TotalEventDisplays_f4 << std::endl;
	cout<<"==============================================================="<<endl;





        // Normalize histograms by the number of entries
        //if (h1->GetEntries() > 0) h1->Scale(1.0 * TotalEventDisplays_f1);
        //if (h2->GetEntries() > 0) h2->Scale(1.0 / TotalEventDisplays_f2);
        //if (h3->GetEntries() > 0) h3->Scale(1.0 / TotalEventDisplays_f3);
	//if (h4->GetEntries() > 0) h4->Scale(1.0 / TotalEventDisplays_f4);

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

	//if(SRUx1 != SRUx2) h1->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//if(SRUx1 != SRUx2) h2->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//if(SRUx1 != SRUx2) h3->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );
	//if(SRUx1 != SRUx2) h4->GetXaxis()->SetRangeUser( SRUx1, SRUx2 );

        h1->GetYaxis()->SetRangeUser( 1.e-7, 0.1);
        h2->GetYaxis()->SetRangeUser( 1.e-7, 0.1);
        //h1->GetYaxis()->SetRangeUser( 0.0, 0.04);   // linear
        //h2->GetYaxis()->SetRangeUser( 0.0, 0.04);

        // Print mean and standard deviation
//	std::cout << std::fixed << std::setprecision(2);
//        std::cout << /*"Histogram from " << file1 <<*/ "  Mean: " << h1->GetMean() << " ± " << h1->GetStdDev() << std::endl;
//        std::cout << /*"Histogram from " << file2 <<*/ "  Mean: " << h2->GetMean() << " ± " << h2->GetStdDev() << std::endl;
//        std::cout << /*"Histogram from " << file3 <<*/ "  Mean: " << h3->GetMean() << " ± " << h3->GetStdDev() << std::endl;
       
//        std::cout << std::fixed << std::setprecision(2);
//        std::cout << /*"Histogram from " << file1 <<*/ "  Mean: " << h1->GetMean() << " ± " << h1->GetEntries() << std::endl;
//        std::cout << /*"Histogram from " << file2 <<*/ "  Mean: " << h2->GetMean() << " ± " << h2->GetEntries() << std::endl;
//        std::cout << /*"Histogram from " << file3 <<*/ "  Mean: " << h3->GetMean() << " ± " << h3->GetEntries() << std::endl;

        std::cout << /*"Histogram from " << file1 <<*/ "  Mean ± [Std / sqrt(Entries)]  " << std::endl;
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "EXT, 0p, Np, 0p+Np : " << std::endl;	
        std::cout << "  Mean: " << h1->GetMean() << " ± " << h1->GetStdDev() / sqrt(TotalEventDisplays_f1) << std::endl;
        std::cout << "  Mean: " << h3->GetMean() << " ± " << h3->GetStdDev() / sqrt(TotalEventDisplays_f2) << std::endl;
        std::cout << "  Mean: " << h4->GetMean() << " ± " << h4->GetStdDev() / sqrt(TotalEventDisplays_f3) << std::endl;
        std::cout << "  Mean: " << h2->GetMean() << " ± " << h2->GetStdDev() / sqrt(TotalEventDisplays_f4) << std::endl;
	std::cout << "\n" << std::endl;



        // Create canvas and draw histograms
        TCanvas *c = new TCanvas("c", "Normalized Histograms", 1200, 900);

	h1->SetLineColor(kRed-3) ;
	h1->SetMarkerColor(kRed-3) ;
	h1->SetMarkerStyle(20) ;
	h1->SetLineWidth(3);

	h2->SetLineColor(kBlue-6);
	h2->SetMarkerColor(kBlue-6) ;
	h2->SetMarkerStyle(20) ;
        h2->SetLineWidth(3);

	h3->SetLineColor(8);
	h3->SetMarkerColor(8) ;
	h3->SetMarkerStyle(20) ;

        h4->SetLineColor(4);
        h4->SetMarkerColor(4) ;
        h4->SetMarkerStyle(20) ;


        h1->Draw("E");
	//h1->Draw("HIST");
	//h4->Draw("E SAME");
	//h3->Draw("E SAME");
	h2->Draw("E SAME");
	//h2->Draw("HIST SAME");
        
        // Add legend
        TLegend *legend = new TLegend(0.5, 0.6, 0.9, 0.9);//0.25,0.2,0.55,0.5);//
        legend->SetFillColorAlpha(0,0.0);
        legend->SetLineColorAlpha(0,0.0);
 	legend->AddEntry(h1, "Run 4b BNB nu overlay", "lpe");   // (Blips Saved in Reco2)
//      legend->AddEntry(h3, "WC 1g Incl. 0p", "lpe");
//	legend->AddEntry(h4, "WC 1g Incl. Np", "lpe");
	legend->AddEntry(h2, "Run 3 BNB nu overlay", "lpe");
        legend->Draw();


	// Get mean and standard deviation
    	double mean1 = h1->GetMean();
	double stdE1 = h1->GetStdDev() / sqrt(h1->GetEntries()) ; // StdE = Std / sqrt (Entries)
    
    	double mean2 = h2->GetMean();
    	double stdE2 = h2->GetStdDev() / sqrt(h2->GetEntries()) ;
    
    	double mean3 = h3->GetMean();
    	double stdE3 = h3->GetStdDev() / sqrt(h3->GetEntries());
    
    	double mean4 = h4->GetMean();
    	double stdE4 = h4->GetStdDev() / sqrt(h4->GetEntries()) ;
    
    // Create TLatex text box for each histogram
    TLatex text;
    text.SetTextSize(0.04);
    text.SetTextFont(42);    // Set bold Helvetica font
    if(printBox){ 
    text.SetTextColor(kRed);
    text.DrawLatexNDC(0.6, 0.55, Form("Mean: %.3f #pm %.3f", mean1, stdE1));

    text.SetTextColor(8);
    text.DrawLatexNDC(0.6, 0.50, Form("Mean: %.2f #pm %.2f", mean3, stdE3));

    text.SetTextColor(kBlue);
    text.DrawLatexNDC(0.6, 0.45, Form("Mean: %.2f #pm %.2f", mean4, stdE4));

//    text.SetTextColor(kMagenta);
//    text.DrawLatexNDC(0.6, 0.40, Form("Mean: %.2f #pm %.2f", mean2, stdE2));

    }











        // Show canvas
        c->Update();
        c->Print(Form("Val_Gen_Lin/png/%s_%s_%s.png", histTitle, histName, tag )) ;
        //c->Print(Form("Val_Gen_Lin/pdf/%s_%s_%s.pdf", histTitle, histName, tag )) ;
 
       	// Set log scale on y-axis
        gPad->SetLogy(1); 
 	c->Update();
        c->Print(Form("Val_Gen_Log/png/%s_%s_%s.png", histTitle, histName, tag )) ;
        //c->Print(Form("Val_Gen_Log/pdf/%s_%s_%s.pdf", histTitle, histName, tag )) ;



    
    }
