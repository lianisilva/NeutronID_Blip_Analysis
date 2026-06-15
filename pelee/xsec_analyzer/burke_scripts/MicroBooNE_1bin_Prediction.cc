#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TParameter.h"
#include "TLegend.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

const double NUWRO_TOTAL_TRIGGERS = 278921486;
const double NUWRO_TOTAL_POT = 1.154657e21;
const double MicroBooNE_Closure_Test_total_pot = 6.931399e21;
const double EXT_TRIGGERS = 772546236;

std::vector<double> nuwro_triggers = {
    74298913,
    33182039,
    76943858,
    83435996,
    11060680
};

std::vector<double> nuwro_pot = {
    3.08558e20,
    1.37908e20,
    3.153737e20,
    3.47069e20,
    4.5749e19
};

std::vector<double> ext_triggers = {
    65744587,
    153905891,
    207000000,
    238429356,
    107466402,
};

std::vector<double> MicroBooNE_Closure_Test_triggers = {
    314833304,
    244890419,
    325206401,
    553974587,
    234230189
};

std::vector<double> MicroBooNE_Closure_Test_pot = {
    1.30748e21,
    1.01779e21,
    1.33294e21,
    2.30437e21,
    9.68819e20
};

// Overlay files
std::vector<std::string> input_files = {
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_overlay_POST_multisim.root"
};

// Dirt files
std::vector<std::string> dirt_files = {
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_dirt_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_dirt_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_dirt_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_dirt_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_dirt_POST_multisim.root",
};

// NuWro fake data files
std::vector<std::string> nuwro_files = {
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2b_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root",
  /*"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_overlay_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_overlay_POST_multisim.root"*/
};

// EXT files
std::vector<std::string> ext_files = {
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_EXT_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_EXT_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_EXT_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_EXT_POST_multisim.root",
  "/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_EXT_POST_multisim.root",
};

double SumPOTFromFiles(const std::vector<std::string>& file_paths) {
  double total_pot = 0.0;
  for (const auto& file_path : file_paths) {
    TFile* file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "[ERROR] Could not open file: " << file_path << std::endl;
      continue;
    }

    TParameter<float>* pot_param =
      dynamic_cast<TParameter<float>*>(file->Get("summed_pot"));
    if (pot_param) {
      double pot_val = pot_param->GetVal();
      std::cout << "[INFO] File: " << file_path << " | POT: " << pot_val << std::endl;
      total_pot += pot_val;
    } else {
      std::cerr << "[WARN] 'summed_pot' not found in " << file_path << std::endl;
    }

    file->Close();
    delete file;
  }

  std::cout << "[INFO] Total POT summed across all files: " << total_pot << std::endl;
  return total_pot;
}

void plot_cross_section_summary(double selected_counter,
                                double selected_signal_counter,
                                double background_counter,
                                double truth_counter,
                                double data_counter,
                                double mc_pot,
                                double data_pot)
{
  double N_Ar_FV = 8.05e+29;
  double numu_per_cm2_per_POT_in_AV = 7.3762291e-10;

  double pot_scale = data_pot / mc_pot;
  double norm = 1 / (data_pot * numu_per_cm2_per_POT_in_AV * N_Ar_FV);
  double selected_signal = selected_counter - background_counter;
  double eff = static_cast<double>(selected_signal_counter) / truth_counter;

  double xsec_val = selected_signal * norm / eff;
  double xsec_err = std::sqrt(selected_counter * pot_scale) * norm / eff;

  double xsec_truth = truth_counter * norm;
  double xsec_truth_err = xsec_truth * std::sqrt(truth_counter) / truth_counter;

  std::cout << "Selected (unweighted): " << selected_counter << "\n";
  std::cout << "Selected Signal (weighted): " << selected_signal_counter << "\n";
  std::cout << "Background subtracted data (MC): " << selected_signal << "\n";
  std::cout << "Background (weighted): " << background_counter << "\n";
  std::cout << "Total true signal (weighted): " << truth_counter << "\n";
  std::cout << "Efficiency: " << eff << "\n";
  std::cout << "Norm Factor: " << norm << "\n";
  std::cout << "Cross Section (prediction): " << xsec_val
            << " ± " << xsec_err << " cm^2\n";
  std::cout << "Truth Cross Section: " << xsec_truth
            << " ± " << xsec_truth_err << " cm^2\n";

  double background_subtracted_data = data_counter - background_counter;
  double data_norm = 1.0 / (data_pot * numu_per_cm2_per_POT_in_AV * N_Ar_FV);
  double data_xsec_val = background_subtracted_data * data_norm / eff;
  double data_xsec_err =
      data_xsec_val * std::sqrt(background_subtracted_data) / background_subtracted_data;

  std::cout << "Data Selected (NuWro): " << data_counter << "\n";
  std::cout << "Background subtracted Data: " << background_subtracted_data << "\n";
  std::cout << "Cross Section (NuWro fake data): " << data_xsec_val
            << " ± " << data_xsec_err << " cm^2\n";

  TCanvas* c = new TCanvas("c", "Cross Section Plot", 800, 600);

  TH1F* h_xsec = new TH1F("h_xsec", "", 1, 0, 1);
  h_xsec->SetBinContent(1, xsec_val);
  h_xsec->SetBinError(1, xsec_err);
  h_xsec->SetLineColor(kRed);
  h_xsec->SetLineWidth(2);
  h_xsec->SetStats(0);
  h_xsec->SetMinimum(.75 * std::min(xsec_val, data_xsec_val));
  h_xsec->SetMaximum(1.25 * std::max(xsec_val, data_xsec_val));
  h_xsec->GetXaxis()->SetTitle("CC1#muNp0#piNn total #sigma");
  h_xsec->GetYaxis()->SetTitle("#sigma [cm^{2}]");
  h_xsec->GetXaxis()->SetTickLength(0);
  h_xsec->GetXaxis()->SetLabelSize(0);

  TH1F* h_data_xsec = new TH1F("data_h_xsec", "", 1, 0, 1);
  h_data_xsec->SetBinContent(1, data_xsec_val);
  h_data_xsec->SetBinError(1, data_xsec_err);
  h_data_xsec->SetLineColor(kBlack);
  h_data_xsec->SetLineWidth(2);
  h_data_xsec->SetStats(0);

  h_xsec->Draw("E1");
  h_data_xsec->Draw("E1 SAME");

  TLegend* leg1 = new TLegend(0.6, 0.7, 0.88, 0.85);
  leg1->AddEntry(h_xsec, "MicroBooNE Prediction", "l");
  leg1->AddEntry(h_data_xsec, "NuWro Fake Data", "l");
  leg1->Draw("same");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextSize(0.04);
  text.DrawTextNDC(0.4, 0.91, "MicroBooNE Simulation");

  c->SaveAs("summary_xsec.png");
  c->SaveAs("summary_xsec.pdf");
  c->SaveAs("summary_xsec.C");

  TFile* outfile = new TFile("microboone_xsec_closure_test.root", "RECREATE");
  if (!outfile->IsOpen()) {
    std::cerr << "Error: could not create output file!" << std::endl;
    return;
  }
  h_xsec->Write();           // 1-bin MicroBooNE prediction
  h_data_xsec->Write();      // also store NuWro fake-data xsec for convenience
  outfile->Close();
  delete outfile;

  delete c;
}

void MicroBooNE_1bin_Prediction()
{
  gStyle->SetOptStat(0);
  double truth_counter           = 0.0;  // weighted true signal
  double selected_counter        = 0.0;  // total selected (unweighted)
  double selected_signal_counter = 0.0;  // weighted selected signal
  double background_counter      = 0.0;  // weighted selected background
  double data_selected_counter   = 0.0;  // selected NuWro fake data (unweighted)
  double ext_selected_counter    = 0.0;
  double dirt_selected_counter    = 0.0;

  double overlay_pot = SumPOTFromFiles(input_files);

  for (size_t i = 0; i < input_files.size(); ++i) {

    const auto& file_path = input_files[i];
    std::cout << "[INFO] Overlay file: " << file_path << std::endl;

    TFile* file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open file " << file_path << std::endl;
      continue;
    }

    auto pot_param = dynamic_cast<TParameter<float>*>(file->Get("summed_pot"));
    if (!pot_param) {
      std::cerr << "Warning: POT not found in " << file_path << std::endl;
      file->Close();
      continue;
    }
    double file_pot = pot_param->GetVal();
    //double pot_scale = MicroBooNE_Closure_Test_pot[i] / file_pot;
    double pot_scale = nuwro_pot[i] / file_pot;

    TTree* tree = nullptr;
    file->GetObject("stv_tree", tree);
    if (!tree) {
      std::cerr << "Error: stv_tree not found in " << file_path << std::endl;
      file->Close();
      continue;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<bool>   isSelected(reader, "CC1muNp0piNn_Selected");
    TTreeReaderValue<bool>   isSignal(reader,   "CC1muNp0piNn_MC_Signal");
    TTreeReaderArray<double> genie_weight(reader,  "weight_TunedCentralValue_UBGenie");
    TTreeReaderArray<double> spline_weight(reader, "weight_splines_general_Spline");
    TTreeReaderArray<double> nAr_w(reader,        "weight_neutron_argon_xsec");

    double temp_truth_counter           = 0.0;
    double temp_selected_counter        = 0.0;
    double temp_selected_signal_counter = 0.0;
    double temp_background_counter      = 0.0;

    while (reader.Next()) {

      double cv = spline_weight[0] * genie_weight[0] * nAr_w[0];
      //if (!(0 <= spline_weight[0] <= 30)) continue;
      //if (!(0 <= genie_weight[0]  <= 30)) continue;
      if (std::isinf(spline_weight[0]) || std::isinf(genie_weight[0])) {
	if (*isSignal) {
	  std::cout << "signal event thaqt has an infite weight" << std::endl;
	  ++truth_counter;
	}
        continue;
      }

      if (*isSignal) {
        temp_truth_counter += cv;
      }

      if (!*isSelected) continue;

      temp_selected_counter += cv;  // unweighted selected count

      if (*isSignal) {
        temp_selected_signal_counter += cv;
      } else {
        temp_background_counter += cv;
      }
    }

    truth_counter           += temp_truth_counter           * pot_scale;
    selected_counter        += temp_selected_counter        * pot_scale;
    selected_signal_counter += temp_selected_signal_counter * pot_scale;
    background_counter      += temp_background_counter      * pot_scale;

    file->Close();
  }

  //Dirt
  for (size_t i = 0; i < dirt_files.size(); ++i) {

    const auto& file_path = dirt_files[i];
    std::cout << "[INFO] Dirt file: " << file_path << std::endl;

    TFile* file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open DIRT file " << file_path << std::endl;
      continue;
    }

    auto pot_param = dynamic_cast<TParameter<float>*>(file->Get("summed_pot"));
    if (!pot_param) {
      std::cerr << "Warning: POT not found in " << file_path << std::endl;
      file->Close();
      continue;
    }

    double file_pot  = pot_param->GetVal();
    double pot_scale = nuwro_pot[i] / file_pot;
    //double pot_scale = MicroBooNE_Closure_Test_pot[i] / file_pot;

    TTree* tree = nullptr;
    file->GetObject("stv_tree", tree);
    if (!tree) {
      std::cerr << "Error: stv_tree not found in DIRT file " << file_path << std::endl;
      file->Close();
      continue;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<bool>   isSelected(reader, "CC1muNp0piNn_Selected");
    TTreeReaderArray<double> genie_weight(reader,  "weight_TunedCentralValue_UBGenie");
    TTreeReaderArray<double> spline_weight(reader, "weight_splines_general_Spline");
    TTreeReaderArray<double> nAr_w(reader,        "weight_neutron_argon_xsec");

    double temp_selected_counter   = 0.0;
    double temp_background_counter = 0.0;
    double temp_dirt_counter = 0.0;

    while (reader.Next()) {
      double cv = spline_weight[0] * genie_weight[0] * nAr_w[0];

      if (!*isSelected) continue;

      temp_selected_counter   += cv;
      temp_background_counter += cv;   // all selected dirt is background
      temp_dirt_counter += cv;
    }

    selected_counter   += temp_selected_counter   * pot_scale;
    background_counter += temp_background_counter * pot_scale;
    dirt_selected_counter       += temp_dirt_counter * pot_scale;

    file->Close();
  }

  //Nuwro FD
  for (const auto& file_path : nuwro_files) {

    std::cout << "[INFO] NuWro fake data file: " << file_path << std::endl;

    TFile* file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open file " << file_path << std::endl;
      continue;
    }

    TTree* tree = nullptr;
    file->GetObject("stv_tree", tree);
    if (!tree) {
      std::cerr << "Error: stv_tree not found in " << file_path << std::endl;
      file->Close();
      continue;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<bool> isSelected(reader, "CC1muNp0piNn_Selected");

    while (reader.Next()) {
      if (!*isSelected) continue;
      data_selected_counter += 1.0;  // unweighted NuWro selected
    }

    file->Close();
  }

  //EXT
  for (size_t i = 0; i < ext_files.size(); ++i) {

    const auto& file_path = ext_files[i];
    std::cout << "[INFO] EXT file: " << file_path << std::endl;

    //double trigger_scale = MicroBooNE_Closure_Test_triggers[i] / ext_triggers[i];
    double trigger_scale = nuwro_triggers[i] / ext_triggers[i];

    TFile* file = TFile::Open(file_path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open EXT file " << file_path << std::endl;
      continue;
    }

    TTree* tree = nullptr;
    file->GetObject("stv_tree", tree);
    if (!tree) {
      std::cerr << "Error: stv_tree not found in EXT file " << file_path << std::endl;
      file->Close();
      continue;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<bool> isSelected(reader, "CC1muNp0piNn_Selected");

    double temp_selected_counter   = 0.0;
    double temp_background_counter = 0.0;
    double temp_ext_selected_counter = 0.0;

    while (reader.Next()) {
      if (!*isSelected) continue;
      temp_selected_counter   += 1.0;
      temp_background_counter += 1.0;  // EXT background treated with unit weight
      temp_ext_selected_counter += 1.0;
    }

    selected_counter   += temp_selected_counter   * trigger_scale;
    background_counter += temp_background_counter * trigger_scale;
    ext_selected_counter += temp_ext_selected_counter * trigger_scale;

    file->Close();
  }

  std::cout << "\n[INFO] Final counters (after scaling):\n"
            << "  truth_counter           = " << truth_counter << "\n"
            << "  selected_counter        = " << selected_counter << "\n"
            << "  selected_signal_counter = " << selected_signal_counter << "\n"
            << "  background_counter      = " << background_counter << "\n"
	    << "  MC background_counter      = " << background_counter - ext_selected_counter << "\n"
            << "  data_selected_counter   = " << data_selected_counter + ext_selected_counter << "\n"
            << "  ext_selected_counter    = " << ext_selected_counter << "\n"
            << "  dirt_counter		  = " << dirt_selected_counter << "\n";

  plot_cross_section_summary(selected_counter,
                             selected_signal_counter,
                             background_counter,
                             truth_counter,
                             data_selected_counter,
                             overlay_pot,
                             NUWRO_TOTAL_POT);
}

int main()
{
  MicroBooNE_1bin_Prediction();
  return 0;
}
