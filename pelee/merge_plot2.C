  #include "TROOT.h"
  #include "TLine.h"
  #include "TMarker.h"

  void merge_plot2(TH1D* hist1, const char* leg1, TH1D* hist2, const char* leg2, const char* plotTitle) {

    //hist1->Scale( 1./hist1->Integral());
    //hist2->Scale( 1./hist2->Integral());

    double minY, maxY, minX, maxX;
    minY = std::min(hist1->GetMinimum(), hist2->GetMinimum());
    maxY = std::max(hist1->GetMaximum(), hist2->GetMaximum());

    std::cout<<"\nMinYaxis: "<<minY<<"\nMaxYaxis: "<<maxY<<"\n";

    // Set the Y-axis range on the first histogram before drawing
    hist1->SetMinimum(minY);
    hist1->SetMaximum(maxY+(maxY/10));

    //hist1->Sumw2(); hist2->Sumw2();

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    TH1D* h1_drawn = (TH1D*)hist1->DrawCopy("hist");
    h1_drawn->SetTitle(plotTitle);
    //h1_drawn->SetLineColor(kRed);
    //h1_drawn->SetFillColorAlpha(kRed, 0.55);
    //h1_drawn->SetMarkerColor(kRed);
    h1_drawn->SetLineColor(kBlue-4);
    h1_drawn->SetFillColor(0);
    //h1_drawn->SetMarkerColor(kGreen+2);
    h1_drawn->SetLineWidth(3);

    TH1D* h2_drawn = (TH1D*)hist2->DrawCopy("hist same");
    h2_drawn->SetLineColor(kRed+2);
    h2_drawn->SetFillColor(0);
    //h2_drawn->SetMarkerColor(kMagenta+2);
    //h2_drawn->SetLineColor(kBlue);
    //h2_drawn->SetFillColorAlpha(kBlue, 0.55);
    //h2_drawn->SetMarkerColor(kBlue);
    h2_drawn->SetLineWidth(3);

    auto legend = new TLegend(0.42,0.2,0.58,0.4);
    legend->AddEntry(h1_drawn, leg1, "lf");
    legend->AddEntry(h2_drawn, leg2, "lf");
    legend->Draw();

  }



