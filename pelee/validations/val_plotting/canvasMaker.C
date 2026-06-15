

void canvasMaker(const char* filename = "nu_v_out.root") {   // change this to run on different files
    // open input file
    TFile* infile = new TFile(filename, "READ");
    if (!infile || infile->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << filename << std::endl;
        return;
    }

    infile->ls();

infile->GetListOfKeys()->Print();
    
    // output file
    TFile* outfile = new TFile("eight_histograms.root", "RECREATE");
    
    // canvas divided into 4x2 grid (4 columns, 2 rows) for 8 histograms
    TCanvas* canvas = new TCanvas("multiHist1", "Muliple Histogram Canvas", 1600, 800);
    canvas->Divide(4, 2);

    TH1* hist1[8]; TH1* hist2[8];
    
    TString hist1_names[8] = {"BlipsPerEvent", "blipX", "blipX", "blipY", "blipZ", "blipZY", "blipdx", "blipdw"};
    TString hist2_names[8] = {"blipE", "charge", "nPlanes", "proxTrkDist", "proxTrkID", "touchTrkID", "trueE", "blipZYdistToVtx"};

    // Retrieve histograms from input file
    for(int i = 0; i < 8; i++) {
        hist1[i] =(TH1*)infile->Get(hist1_names[i]);
        if (!hist1[i]) {
            std::cerr << "Error: Cannot find histogram " << hist1_names[i] << " in input file" << std::endl;
            infile->Close();
            outfile->Close();
            return;
        }
        
        // Set different colors for each histogram (optional)
        //hist[i]->SetFillColor(i + 2);
        //hist[i]->SetLineColor(i + 2);
    }

   for(int i = 0; i < 8; i++) {
        hist2[i] = (TH1*)infile->Get(hist2_names[i]);
        if (!hist2[i]) {
            std::cerr << "Error: Cannot find histogram " << hist2_names[i] << " in input file" << std::endl;
            infile->Close();
            outfile->Close();
            return;
        }
    }
    
    for(int i = 0; i < 8; i++) {
        canvas->cd(i + 1); // ROOT canvas pads are numbered starting from 1
        hist1[i]->Draw("HIST");
	//gPad->Update();
    }
    
    // update canvas to ensure all plots are drawn
    canvas->Update();
    for(int i = 0; i < 8; i++) { if (hist1[i]) { hist1[i]->Write("HIST"); } }
    // Save canvas as image files
    //canvas->SaveAs("eight_histograms.png");
    //canvas->SaveAs("eight_histograms.pdf"); 
    canvas->Write();

    TCanvas* canvas2 = new TCanvas("multiHist2", "Muliple Histogram Canvas", 1600, 800);
    canvas2->Divide(4, 2);
    
    for(int i = 0; i < 8; i++) {
        canvas2->cd(i + 1);
        hist2[i]->Draw("HIST");
	//gPad->Update();
        //hist2[i]->Write("HIST");
    }

    canvas2->Update();
    for(int i = 0; i < 8; i++) { if (hist2[i]) { hist2[i]->Write(); } }
    canvas2->Write();
 
    // cleanup
    infile->Close();
    outfile->Close();
    delete infile;
    delete outfile;
    
}
