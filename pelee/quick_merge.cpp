#include "merge_plot2.C"

void quick_merge(){
    TFile* f = new TFile("test_out.root", "READ");
    TH1D* h_data = (TH1D*)f->Get("shower_plots/shrBlipE_data");
    TH1D* h_MC   = (TH1D*)f->Get("shower_plots/shrBlipE_MC");
    
    merge_plot2(h_data, "Data", h_MC, "MC", "Energy Spectrum of Blips Inside Showers");
}
