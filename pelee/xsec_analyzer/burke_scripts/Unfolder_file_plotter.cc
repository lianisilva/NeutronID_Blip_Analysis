#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMath.h"
#include <iostream>
#include <vector>
#include <map>

void extract_and_style_histogram() {
  //TFile* file = TFile::Open("/exp/uboone/app/users/birwin/LiangRepo/xsec_analyzer/CC1muNp0piNn_unfolder_output/Thesis_unfold.root", "READ");
  //TFile* file = TFile::Open("/exp/uboone/app/users/birwin/LiangRepo/xsec_analyzer/CC1muNp0piNn_unfolder_output/Thesis_singlebin_BNB_Data.root", "READ");
  TFile* file = TFile::Open("/exp/uboone/app/users/mhernan/uboonexsec_pv/dev/xsec_analyzer/CC1muNp0piNn_unfolder_output_burkelast_1bin/CC1muNp0piNn_unfold.root", "READ");
  if (!file || file->IsZombie()) return;

  file->cd("EventCountUnits/EventCategory/");
  TH1D* raw_event_hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_total"));
  if (raw_event_hist) {
    double raw_val = raw_event_hist->GetBinContent(1);
    double raw_err = raw_event_hist->GetBinError(1);
    std::cout << "[EVENT COUNT (UNFOLDED)] " << raw_val << " plus or minus " << raw_err << std::endl;
  } else {
    std::cerr << "[ERROR] Could not load EventCountUnits/EventCategory/EventCategory_total" << std::endl;
  }

  file->cd("XsecUnits/EventCategory/");
  TH1D* hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_total"));
  TH1D* mcStats_xsec_units_hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_MCstats"));
  if (!hist || !mcStats_xsec_units_hist) { file->Close(); return; }
  hist->SetDirectory(nullptr);
  mcStats_xsec_units_hist->SetDirectory(nullptr);
  file->Close();

  TFile* xsecFile = TFile::Open("microboone_xsec_output.root", "READ");
  if (!xsecFile || xsecFile->IsZombie()) return;
  TH1F* h_xsec = dynamic_cast<TH1F*>(xsecFile->Get("h_xsec"));
  if (!h_xsec) { xsecFile->Close(); return; }
  h_xsec->SetDirectory(nullptr);
  xsecFile->Close();

  double xsec_value = h_xsec->GetBinContent(1);
  int nBins = hist->GetNbinsX();
  double xMin = hist->GetXaxis()->GetXmin();
  double xMax = hist->GetXaxis()->GetXmax();
  TH1D* h_line = new TH1D("h_line", "", nBins, xMin, xMax);
  for (int i = 1; i <= nBins; ++i) h_line->SetBinContent(i, xsec_value);
  std::cout << "MicroBooNE Tune: " << xsec_value << std::endl;
  h_line->SetLineColor(kRed + 1);
  h_line->SetLineWidth(2);
  h_line->SetLineStyle(2);

  TFile* signalFile = TFile::Open("ClosureTest_signal_counts.root", "READ");
  if (!signalFile || signalFile->IsZombie()) return;

  std::vector<std::string> generators = {
    "neut_5_4_0_1",
    "AR23_20i_00_000",
    //"gibuu_2025",
    "NuWro_numu_21_09_2",
    "genie_v3_00_06",
  };

  std::vector<std::string> titles = {
    "NEUT 5.4.0.1",
    "GENIE AR23.20i.00.000",
    //"GiBUU 2025",
    "NuWro 21.09.2",
    "GENIE G18 v3.00.06",
  };

  std::vector<int> colors = {
    3,
    4,
    //6,
    7,
    51,
  };

  std::vector<int> styles = {
    2,
    4,
    //7,
    9,
    10,
  };

  std::vector<TH1D*> model_lines;
  std::vector<double> pvalues(generators.size(), -1.0);
  double pvalue_mb = -1.0;   // MicroBooNE prediction p-value

  for (size_t i = 0; i < generators.size(); ++i) {
    const auto& gen = generators[i];
    std::string hist_name = "h_" + gen;
    TH1D* h_model = dynamic_cast<TH1D*>(signalFile->Get(hist_name.c_str()));
    if (h_model) {
      if (gen == "gibuu_2025") {
        h_model->Scale(1.0 / 250.0);
      }
      h_model->SetDirectory(nullptr);
      h_model->SetLineColor(colors[i]);
      h_model->SetLineStyle(styles[i]);
      h_model->SetLineWidth(2);
      model_lines.push_back(h_model);
    } else {
      model_lines.push_back(nullptr);
      std::cerr << "[WARN] Could not find histogram " << hist_name << " in signal_counts.root" << std::endl;
    }
  }
  signalFile->Close();

  hist->GetXaxis()->SetLabelSize(0);
  hist->GetXaxis()->SetTickLength(0);
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(kFullCircle);
  hist->SetMarkerSize(1.0);
  hist->SetTitle("");
  hist->GetXaxis()->SetTitle("");
  hist->GetYaxis()->SetTitle("#sigma (cm^{2}/nucleus)");

  double unfolded_val = hist->GetBinContent(1);
  std::cout << "Unfolded value: " << unfolded_val << std::endl;
  double unfolded_err = hist->GetBinError(1);
  std::cout << "Unfolded err: " << unfolded_err << std::endl;
  double max_val = unfolded_val + unfolded_err;
  double min_val = unfolded_val - unfolded_err;

  for (TH1D* h : model_lines) {
    if (!h) continue;
    double val = h->GetBinContent(1);
    if (val > max_val) max_val = val;
    if (val < min_val) min_val = val;
  }
  hist->SetMaximum(1.15 * max_val);
  hist->SetMinimum(0.0);

  TCanvas* c = new TCanvas("c", "Styled Histogram", 800, 600);
  hist->Draw("E1");

  /*TH1D* h_mcstat_band = (TH1D*)mcStats_xsec_units_hist->Clone("h_mcstat_band");
  h_mcstat_band->SetDirectory(nullptr);
  h_mcstat_band->SetFillColor(kGray+3);
  h_mcstat_band->SetFillStyle(3004);
  h_mcstat_band->SetLineColor(kGray+2);
  h_mcstat_band->SetLineWidth(1);
  h_mcstat_band->Draw("E2 same");*/

  for (size_t i = 0; i < generators.size(); ++i) {
    TH1D* model_hist = model_lines[i];
    if (!model_hist) continue;

    double val = model_hist->GetBinContent(1);
    std::cout << generators[i] << " cross section: " << val << " (10^-39 cm^2/nucleus)" << std::endl;

    double p_val = -1.0;
    if (unfolded_err > 0.0) {
      double z = (val - unfolded_val) / unfolded_err;
      double z_abs = std::fabs(z);
      p_val = TMath::Erfc(z_abs / std::sqrt(2.0));  // two-sided p-value
    }
    pvalues[i] = p_val;

    std::cout << "  -> two-sided p-value vs unfolded xsec: " << p_val << std::endl;

    TH1D* flat_line = new TH1D(("line_" + generators[i]).c_str(), "", nBins, xMin, xMax);
    for (int j = 1; j <= nBins; ++j) flat_line->SetBinContent(j, val);
    flat_line->SetLineColor(colors[i]);
    flat_line->SetLineStyle(styles[i]);
    flat_line->SetLineWidth(2);
    flat_line->Draw("same HIST");
  }

  // p-value for MicroBooNE prediction vs unfolded value
  if (unfolded_err > 0.0) {
    double z_mb = (xsec_value - unfolded_val) / unfolded_err;
    double z_mb_abs = std::fabs(z_mb);
    double p_mb = TMath::Erfc(z_mb_abs / std::sqrt(2.0));
    pvalue_mb = p_mb;
    std::cout << "MicroBooNE prediction p-value vs unfolded xsec: " << p_mb << std::endl;
  }


  h_line->Draw("same HIST");

  TLegend* legend = new TLegend(0.58, 0.65, 0.85, 0.85);
  legend->AddEntry(hist, "Unfolded BNB Data (1.15e21 POT)", "lpe");
  //legend->AddEntry(h_mcstat_band, "MC Stat. Uncertainty", "f");
  {
    TString label_mb = Form("MicroBooNE Prediction (p = %.3f)", pvalue_mb);
    legend->AddEntry(h_line, label_mb, "l");
  }
  for (size_t i = 0; i < model_lines.size(); ++i) {
    if (!model_lines[i]) continue;
    TString label = Form("%s (p = %.3f)", titles[i].c_str(), pvalues[i]);
    legend->AddEntry(model_lines[i], label, "l");
  }
  legend->SetBorderSize(1);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(1001);
  legend->Draw();

  c->SaveAs("Thesis_singlebin_xsec_no_gibuu.pdf");
}

int main() {
  extract_and_style_histogram();
  return 0;
}

