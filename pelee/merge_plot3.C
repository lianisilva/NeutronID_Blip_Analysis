  #include "TROOT.h"
  #include "TLine.h"
  #include "TMarker.h"

  void merge_plot3(TH1D* hist1, const char* leg1, TH1D* hist2, const char* leg2, TH1D* hist3, const char* leg3) {

    //hist1->Scale( 1./hist1->Integral());
    //hist2->Scale( 1./hist2->Integral());

    double minY, maxY, minX, maxX;
    minY = std::min(hist1->GetMinimum(), hist2->GetMinimum());
    maxY = std::max(hist1->GetMaximum(), hist2->GetMaximum());
    
    std::cout<<"\nMinYaxis: "<<minY<<"\nMaxYaxis: "<<maxY<<"\n";
    
    // Set the Y-axis range on the first histogram before drawing
    hist1->SetMinimum(0);
    hist1->SetMaximum(maxY+10);
    hist3->SetMinimum(0);
    hist3->SetMaximum(maxY+10);
    
    hist1->Sumw2();
    hist2->Sumw2();
    hist3->Sumw2();

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    auto legend = new TLegend(0.42,0.2,0.58,0.4);
    //hist1->SetLineColor(kRed);
    //hist1->SetFillColorAlpha(kRed, 0.55);
    //hist1->SetMarkerColor(kRed);
    hist1->SetLineColor(kGreen+2);
    hist1->SetFillColor(0);
    //hist1->SetMarkerColor(kGreen+2);
    hist1->SetLineWidth(3);
    hist1->DrawCopy("HIST");
    legend->AddEntry(hist1, leg1, "lf");

    hist2->SetLineColor(kMagenta+2);
    hist2->SetFillColor(0);
    //hist2->SetMarkerColor(kMagenta+2);
    //hist2->SetLineColor(kBlue);
    //hist2->SetFillColorAlpha(kBlue, 0.55);
    //hist2->SetMarkerColor(kBlue);
    hist2->SetLineWidth(3);
    hist2->DrawCopy("HIST SAME");
    legend->AddEntry(hist2, leg2, "lf");

    hist3->SetLineColor(kOrange+2);
    hist3->SetFillColor(0);
    hist3->SetLineWidth(3);
    hist3->DrawCopy("HIST SAME");
    legend->AddEntry(hist3, leg3, "lf");
    legend->Draw();

  }



