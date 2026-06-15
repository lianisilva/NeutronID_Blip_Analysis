#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TAxis.h"
#include <iostream>
#include <vector>

void extract_and_style_histogram() {
    //TFile* file = TFile::Open("/exp/uboone/app/users/birwin/LiangRepo/xsec_analyzer/CC1muNp0piNn_unfolder_output/CC1muNp0piNn_unfold_nAr_weight.root", "READ");
    //TFile* file = TFile::Open("/exp/uboone/app/users/birwin/LiangRepo/xsec_analyzer/CC1muNp0piNn_unfolder_output/MicroBooNE_Prediction_Fake_data_closure_single_bin_Unfolder.root", "READ");
    TFile* file = TFile::Open("/exp/uboone/app/users/mhernan/uboonexsec_pv/dev/xsec_analyzer/CC1muNp0piNn_unfolder_output_burkelast_1bin/CC1muNp0piNn_unfold.root", "READ");

    if (!file || file->IsZombie()) return;

    file->cd("EventCountUnits/EventCategory/");
    TH1D* DataStats_hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_DataStats"));
    TH1D* xsec_total_hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_xsec_total"));
    if (DataStats_hist && xsec_total_hist) {
      double raw_val = xsec_total_hist->GetBinContent(1);
      double xsec_err = xsec_total_hist->GetBinError(1);
      double stat_err = DataStats_hist->GetBinError(1);
      double raw_err = std::sqrt(std::pow(xsec_err, 2) + std::pow(stat_err, 2));
      std::cout << "[EVENT COUNT (UNFOLDED)] " << raw_val << " plus or minus " << raw_err << std::endl;
    } else {
      std::cerr << "[ERROR] Could not load EventCountUnits/EventCategory/EventCategory_total" << std::endl;
    }

    file->cd("XsecUnits/EventCategory/");
    TH1D* dataStats_xsec_units_hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_DataStats"));
    TH1D* xsec_total_xsec_units_hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_xsec_total"));
    TH1D* mcStats_xsec_units_hist = dynamic_cast<TH1D*>(gDirectory->Get("EventCategory_MCstats"));

    if (!dataStats_xsec_units_hist || !xsec_total_xsec_units_hist || !mcStats_xsec_units_hist) { file->Close(); return; }
    dataStats_xsec_units_hist->SetDirectory(nullptr);
    xsec_total_xsec_units_hist->SetDirectory(nullptr);
    mcStats_xsec_units_hist->SetDirectory(nullptr);

    file->Close();

    //TFile* xsecFile = TFile::Open("microboone_xsec_output.root", "READ");
    TFile* xsecFile = TFile::Open("microboone_xsec_closure_test.root", "READ");
    if (!xsecFile || xsecFile->IsZombie()) return;
    TH1F* h_xsec = dynamic_cast<TH1F*>(xsecFile->Get("h_xsec"));
    if (!h_xsec) { xsecFile->Close(); return; }
    h_xsec->SetDirectory(nullptr);
    xsecFile->Close();

    double xsec_value = h_xsec->GetBinContent(1);
    int nBins = dataStats_xsec_units_hist->GetNbinsX();
    double xMin = dataStats_xsec_units_hist->GetXaxis()->GetXmin();
    double xMax = dataStats_xsec_units_hist->GetXaxis()->GetXmax();
    TH1D* h_line = new TH1D("h_line", "", nBins, xMin, xMax);
    for (int i = 1; i <= nBins; ++i) h_line->SetBinContent(i, xsec_value);
    std::cout << "MicroBooNE Tune: " << xsec_value << std::endl;
    h_line->SetLineColor(kRed + 1);
    h_line->SetLineWidth(2);
    h_line->SetLineStyle(2);

    TFile* signalFile = TFile::Open("ClosureTest_signal_counts.root", "READ");
    if (!signalFile || signalFile->IsZombie()) return;
    std::vector<std::string> generators = {
      //"genie_v2_12_10",
      //"genie_v2_12_10.mec",
      //"genie_v3_00_06",
      "neut_5_4_0_1",
      //"nuwro_19_02_1"
      "AR23_20i_00_000",
      "gibuu_2025",
      "NuWro_numu_21_09_2",
      "genie_v3_00_06",
    };

    std::vector<std::string> titles = {
	"NEUT 5.4.0.1",
	"GENIE AR23.20i.00.000",
	"GiBUU 2025",
	"NuWro 21.09.2",
	"GENIE G18 v3.00.06",
    };

    std::vector<int> colors = {
	3,
	4,
	6,
	7,
	51,
    };

    std::vector<int> styles = {
	2,
	4,
	7,
	9,
	10,
    };

    std::vector<TH1D*> model_lines;
    int color_it = 0;
    for (const auto& gen : generators) {
      std::string hist_name = "h_" + gen;
      TH1D* h_model = dynamic_cast<TH1D*>(signalFile->Get(hist_name.c_str()));
      if (h_model) {
	std::cout << "model line" << std::endl;
        if ( gen == "gibuu_2025" ) { h_model->Scale(1.0/250.0); }
        h_model->SetDirectory(nullptr);
        h_model->SetLineColor(colors[color_it]);
        h_model->SetLineStyle(styles[color_it]);
        h_model->SetLineWidth(2);
        model_lines.push_back(h_model);
        color_it++;
      }
    }
    signalFile->Close();

    dataStats_xsec_units_hist->GetXaxis()->SetLabelSize(0);
    dataStats_xsec_units_hist->GetXaxis()->SetTickLength(0);
    dataStats_xsec_units_hist->SetLineColor(kBlue + 1);
    dataStats_xsec_units_hist->SetLineWidth(2);
    dataStats_xsec_units_hist->SetMarkerStyle(kFullCircle);
    dataStats_xsec_units_hist->SetMarkerSize(1.0);
    dataStats_xsec_units_hist->SetTitle("");
    dataStats_xsec_units_hist->GetXaxis()->SetTitle("");
    dataStats_xsec_units_hist->GetYaxis()->SetTitle("#sigma (cm^{2}/nucleus)");

    double unfolded_val = dataStats_xsec_units_hist->GetBinContent(1);
    std::cout << "Unfolded value: " << unfolded_val << std::endl;
    double dataStats_xsec_units_err = dataStats_xsec_units_hist->GetBinError(1);
    double xsec_total_xsec_units_err = xsec_total_xsec_units_hist->GetBinError(1);
    double unfolded_err = std::sqrt(std::pow(xsec_total_xsec_units_err, 2) + std::pow(dataStats_xsec_units_err, 2));
    std::cout << "Unfolded err: " << unfolded_err << std::endl;
    double max_val = unfolded_val + unfolded_err;
    double min_val = unfolded_val - unfolded_err;
    for (TH1D* h : model_lines) {
      double val = h->GetBinContent(1);
      if (val > max_val) max_val = val;
      if (val < min_val) min_val = val;
    }
    dataStats_xsec_units_hist->SetMaximum(1.15 * max_val);
    dataStats_xsec_units_hist->SetMinimum(0.0);
    dataStats_xsec_units_hist->SetBinError(1, unfolded_err);

    TCanvas* c = new TCanvas("c", "Styled Histogram", 800, 600);
    dataStats_xsec_units_hist->Draw("E1");

    std::cout << "before mc stat band" << std::endl;

    TH1D* h_mcstat_band = (TH1D*)mcStats_xsec_units_hist->Clone("h_mcstat_band");
    h_mcstat_band->SetDirectory(nullptr);
    h_mcstat_band->SetFillColor(kGray+3);
    h_mcstat_band->SetFillStyle(3004);
    h_mcstat_band->SetLineColor(kGray+2);
    h_mcstat_band->SetLineWidth(1);
    //h_mcstat_band->Draw("E2 same");

    std::cout << "after mc stat band" << std::endl;

    for (size_t i = 0; i < generators.size(); ++i) {
      double val = model_lines[i]->GetBinContent(1);
      std::cout << generators[i] << " cross section: " << val << " (10^-39 cm^2/nucleus)" << std::endl;
      TH1D* flat_line = new TH1D(("line_" + generators[i]).c_str(), "", nBins, xMin, xMax);
      for (int j = 1; j <= nBins; ++j) flat_line->SetBinContent(j, val);
      flat_line->SetLineColor(model_lines[i]->GetLineColor());
      flat_line->SetLineStyle(styles[i]);
      flat_line->SetLineWidth(2);
      flat_line->Draw("same HIST");
    }
    h_line->Draw("same HIST");

    TLegend* legend = new TLegend(0.58, 0.65, 0.85, 0.85);
    legend->AddEntry(dataStats_xsec_units_hist, "Unfolded Nuwro Fake Data", "lpe");
    legend->AddEntry(h_line, "MicroBooNE Prediction", "l");
    //legend->AddEntry(h_mcstat_band, "MC Stat. Uncertainty", "f");
    for (size_t i = 0; i < model_lines.size(); ++i) {
      legend->AddEntry(model_lines[i], titles[i].c_str(), "l");
    }
    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->Draw();

    c->SaveAs("single_bin_xsec_nuwro_FD_onlyxsec_and_Data_unc_applied.pdf");
}

int main() {
  extract_and_style_histogram();
  return 0;
}

