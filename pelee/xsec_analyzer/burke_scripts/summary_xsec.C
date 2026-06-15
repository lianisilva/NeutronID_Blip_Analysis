void summary_xsec()
{
//=========Macro generated from canvas: c/Cross Section Plot
//=========  (Wed Feb  4 15:12:38 2026) by ROOT version 6.12/06
   TCanvas *c = new TCanvas("c", "Cross Section Plot",0,0,800,600);
   gStyle->SetOptStat(0);
   c->SetHighLightColor(2);
   c->Range(-0.125,1.475171e-38,1.125,4.468968e-38);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetFrameBorderMode(0);
   c->SetFrameBorderMode(0);
   
   TH1F *h_xsec__1 = new TH1F("h_xsec__1","",1,0,1);
   h_xsec__1->SetBinContent(1,3.33567e-38);
   h_xsec__1->SetBinError(1,2.675854e-39);
   h_xsec__1->SetMinimum(1.774551e-38);
   h_xsec__1->SetMaximum(4.169588e-38);
   h_xsec__1->SetEntries(1);
   h_xsec__1->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   h_xsec__1->SetLineColor(ci);
   h_xsec__1->SetLineWidth(2);
   h_xsec__1->GetXaxis()->SetTitle("CC1#muNp0#piNn total #sigma");
   h_xsec__1->GetXaxis()->SetLabelFont(42);
   h_xsec__1->GetXaxis()->SetLabelSize(0);
   h_xsec__1->GetXaxis()->SetTitleSize(0.035);
   h_xsec__1->GetXaxis()->SetTickLength(0);
   h_xsec__1->GetXaxis()->SetTitleFont(42);
   h_xsec__1->GetYaxis()->SetTitle("#sigma [cm^{2}]");
   h_xsec__1->GetYaxis()->SetLabelFont(42);
   h_xsec__1->GetYaxis()->SetLabelSize(0.035);
   h_xsec__1->GetYaxis()->SetTitleSize(0.035);
   h_xsec__1->GetYaxis()->SetTitleOffset(0);
   h_xsec__1->GetYaxis()->SetTitleFont(42);
   h_xsec__1->GetZaxis()->SetLabelFont(42);
   h_xsec__1->GetZaxis()->SetLabelSize(0.035);
   h_xsec__1->GetZaxis()->SetTitleSize(0.035);
   h_xsec__1->GetZaxis()->SetTitleFont(42);
   h_xsec__1->Draw("E1");
   
   TH1F *data_h_xsec__2 = new TH1F("data_h_xsec__2","",1,0,1);
   data_h_xsec__2->SetBinContent(1,2.366068e-38);
   data_h_xsec__2->SetBinError(1,3.930207e-39);
   data_h_xsec__2->SetEntries(1);
   data_h_xsec__2->SetStats(0);
   data_h_xsec__2->SetLineWidth(2);
   data_h_xsec__2->GetXaxis()->SetLabelFont(42);
   data_h_xsec__2->GetXaxis()->SetLabelSize(0.035);
   data_h_xsec__2->GetXaxis()->SetTitleSize(0.035);
   data_h_xsec__2->GetXaxis()->SetTitleFont(42);
   data_h_xsec__2->GetYaxis()->SetLabelFont(42);
   data_h_xsec__2->GetYaxis()->SetLabelSize(0.035);
   data_h_xsec__2->GetYaxis()->SetTitleSize(0.035);
   data_h_xsec__2->GetYaxis()->SetTitleOffset(0);
   data_h_xsec__2->GetYaxis()->SetTitleFont(42);
   data_h_xsec__2->GetZaxis()->SetLabelFont(42);
   data_h_xsec__2->GetZaxis()->SetLabelSize(0.035);
   data_h_xsec__2->GetZaxis()->SetTitleSize(0.035);
   data_h_xsec__2->GetZaxis()->SetTitleFont(42);
   data_h_xsec__2->Draw("E1 SAME");
   
   TLegend *leg = new TLegend(0.6,0.7,0.88,0.85,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("h_xsec","MicroBooNE Prediction","l");

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("data_h_xsec","NuWro Fake Data","l");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TText *text = new TText(0.4,0.91,"MicroBooNE Simulation");
   text->SetNDC();
   text->SetTextFont(42);
   text->SetTextSize(0.04);
   text->Draw();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
