#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <map>
#include <algorithm>   // for std::find
#include <functional>  // for std::function
#include <cassert>
#include <numeric>

#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TKey.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "Rtypes.h"
#include "TLatex.h"
#include "TParameter.h"
#include "TLine.h"

enum class Sample { kOverlay, kDirt, kEXT, kNuWro };

namespace {
constexpr const char kLabel_Cosmic[] = "Cosmic";
constexpr const char kLabel_Dirt[]   = "Dirt";
constexpr const char kLabel_EXT[]    = "EXT";
}

const std::map<std::string, std::vector<int>> kGroups = {
  {"#nu_{#mu}CCNnNp0#pi",       {5, 6, 8, 9}},
  {"#nu_{#mu}CC0#pi nonsignal", {1, 2, 3, 4, 7}},
  {"#nu_{#mu}CCN#pi",           {10}},
  {"#nu_{#mu}CCOther",          {11}},
  {"Other",                     {12, 13, 14, 15}}
};

struct RunGroup {
  std::string file;
  int run_index; // 0..4 to align with scaling vectors
};


static const char* kTreeName = "stv_tree";//"nuselection/NeutrinoSelectionFilter";//

// ---------------------------- Scaling reference ----------------------------
//static const std::vector<double> NUWRO_TRIGGERS = {
//    74298913, 33182039, 76943858, 83435996, 11060680
//};
//static const std::vector<double> NUWRO_POT = {
//    3.08558e20, 1.37908e20, 3.153737e20, 3.47069e20, 4.5749e19
//};
static const std::vector<double> BNB_Data_TRIGGERS = {
    33105386, 0, 0, 0, 0  // 37273255, 61882791, 6.14e07, 73669038, 35265730
};
static const std::vector<double> BNB_Data_POT = {
    1.399e20, 0, 0, 0, 0  // 1.67e20, 2.61e20, 2.57e20, 3.20e20, 1.48e20
};
static const std::vector<double> EXT_TRIGGERS = {
    33627260, 0, 0, 0, 0  // 65744587, 153905891, 207000000, 238429356, 107466402
};

// ---------------------------- Volumes --------------------------------------
constexpr double FV_X_MIN =   21.5;
constexpr double FV_X_MAX =  234.85;
constexpr double FV_Y_MIN =  -95.0;
constexpr double FV_Y_MAX =   95.0;
constexpr double FV_Z_MIN =   21.5;
constexpr double FV_Z_MAX =  966.8;

constexpr double PCV_X_MIN =   10.0;
constexpr double PCV_X_MAX =  246.35;
constexpr double PCV_Y_MIN = -106.5;
constexpr double PCV_Y_MAX =  106.5;
constexpr double PCV_Z_MIN =   10.0;
constexpr double PCV_Z_MAX = 1026.8;

static void DrawBox(double xmin, double xmax, double ymin, double ymax,
                    int color, int style, int width = 3) {
  auto* l1 = new TLine(xmin, ymin, xmax, ymin);
  auto* l2 = new TLine(xmax, ymin, xmax, ymax);
  auto* l3 = new TLine(xmax, ymax, xmin, ymax);
  auto* l4 = new TLine(xmin, ymax, xmin, ymin);
  for (auto* L : {l1,l2,l3,l4}) { L->SetLineColor(color); L->SetLineStyle(style); L->SetLineWidth(width); L->Draw(); }
}

static void DrawBoxesXZ() {
  DrawBox(FV_X_MIN,  FV_X_MAX,  FV_Z_MIN,  FV_Z_MAX,  kRed, 1, 4);
  DrawBox(PCV_X_MIN, PCV_X_MAX, PCV_Z_MIN, PCV_Z_MAX, kRed, 2, 3);
}

static void DrawBoxesXY() {
  DrawBox(FV_X_MIN,  FV_X_MAX,  FV_Y_MIN,  FV_Y_MAX,  kRed, 1, 4);
  DrawBox(PCV_X_MIN, PCV_X_MAX, PCV_Y_MIN, PCV_Y_MAX, kRed, 2, 3);
}

static inline double first_or_one(const std::vector<double>* v) {
  if (!v || v->empty()) return 1.0;
  const double x = (*v)[0];
  if (!std::isfinite(x)) return 1.0; 
  if (std::abs(x) > 1e3) return 1.0;   // discard absurd values
  return x;
}

static inline std::string ResolveLegendBucket(Sample st, int cat) {
  if (st == Sample::kEXT) return kLabel_EXT;                         // << here
  if (cat == 0 && (st == Sample::kOverlay || st == Sample::kDirt))
    return kLabel_Cosmic;
  if (st == Sample::kDirt && cat != 0) return kLabel_Dirt;

  if (st == Sample::kOverlay && cat != 0) {
    for (const auto& kv : kGroups) {
      if (std::find(kv.second.begin(), kv.second.end(), cat) != kv.second.end())
        return kv.first;
    }
    return "Other";
  }
  return "Other";
}

static double PerFileScale(Sample s, int run_idx, TFile* f)
{
  if (run_idx < 0 || run_idx >= 5) {
    std::cerr << "[warn] Bad run_idx=" << run_idx << " for scaling. Using scale=1.\n";
    return 1.0;
  }

  if (s == Sample::kEXT) {
    const double num = BNB_Data_TRIGGERS[run_idx];
    const double den = EXT_TRIGGERS[run_idx];
    if (den <= 0) { std::cerr << "[warn] EXT trigger denominator <=0; scale=1.\n"; return 1.0; }
    return num / den;
  }

  if (s == Sample::kOverlay || s == Sample::kDirt) {
    auto* pot_param = dynamic_cast<TParameter<float>*>(f->Get("summed_pot"));
    if (!pot_param) {
      std::cerr << "[warn] 'summed_pot' missing in file '" << f->GetName()
                << "'. Using scale=1.\n";
      return 1.0;
    }
    const double pot_file = pot_param->GetVal();
    if (pot_file <= 0.) {
      std::cerr << "[warn] pot_file<=0 in file '" << f->GetName() << "'. Using scale=1.\n";
      return 1.0;
    }
    return BNB_Data_POT[run_idx] / pot_file;
  }
  return 1.0;
}

using BucketMap = std::map<std::string, std::unique_ptr<TH1F>>;
using GetBucketHistFn = std::function<TH1F*(const std::string&)>;

// Map PDG bucket for Overlay
static inline std::string BucketFrom_true_ncat(int category) {
  const int a = category;
  if (a == 1) return "Primary (n,1p)";
  if (a == 2) return "Primary (n,Np)";
  if (a == 3 || a == 4) return "Secondary (n,Xp)";
  if (a == 5) return "Primary (n,gamma)";
  if (a == 6) return "Secondary (n,gamma)";
  if (a == 7) return "ncapture gamma";
  if (a == -9) return "Cosmics and radiogenics";
  if (a == 0 || a == 8) return "Other (truth-matched)";
  return "Other";
}

static const std::map<std::string, int> kBlipColors = {
  {"Other (truth-matched)",   kCyan+3},
  {"Primary (n,1p)", kGreen-2},
  {"Primary (n,Np)",   kSpring+2},
  {"Dirt",   kAzure-5},
  {"EXT",    kBlue-5},
  {"Secondary (n,Xp)", kPink-5},
  {"Primary (n,gamma)",  kYellow-5},
  {"Secondary (n,gamma)",           kRed-5},
  {"ncapture gamma",           kMagenta-5},
  {"Cosmics and radiogenics",           kViolet-5}
};

static void WriteIntegratedEfficiencyCSV_OverlayOnly(
    const std::vector<RunGroup>& overlay_files,
    const std::vector<std::string>& stage_cuts,  // cumulative order
    const std::string& out_csv,
    Long64_t max_events = -1)
{
  TChain ch(kTreeName);
  for (const auto& rg : overlay_files) ch.Add(rg.file.c_str());

  ch.SetBranchStatus("*", 0);

  Bool_t mc_signal = kFALSE;
  ch.SetBranchStatus("CC1muNp0piNn_MC_Signal", 1);
  ch.SetBranchAddress("CC1muNp0piNn_MC_Signal", &mc_signal);

  const size_t S = stage_cuts.size();
  std::vector<std::unique_ptr<Bool_t>> stage_vals(S);
  for (size_t i = 0; i < S; ++i) {
    stage_vals[i] = std::make_unique<Bool_t>(kFALSE);
    ch.SetBranchStatus(stage_cuts[i].c_str(), 1);
    ch.SetBranchAddress(stage_cuts[i].c_str(), stage_vals[i].get());
  }

  std::vector<double>* w_spline = nullptr;
  std::vector<double>* w_genie  = nullptr;
  std::vector<double>* w_nxsec  = nullptr;
  //ch.SetBranchStatus("weight_splines_general_Spline", 1);
  //ch.SetBranchStatus("weight_TunedCentralValue_UBGenie", 1);
  //ch.SetBranchStatus("weight_neutron_argon_xsec", 1);
  //ch.SetBranchAddress("weight_splines_general_Spline", &w_spline);
  //ch.SetBranchAddress("weight_TunedCentralValue_UBGenie", &w_genie);
  //ch.SetBranchAddress("weight_neutron_argon_xsec", &w_nxsec);

  auto w0 = [](const std::vector<double>* v)->double {
    return (v && !v->empty() && std::isfinite((*v)[0])) ? (*v)[0] : 1.0;
  };

  std::vector<double> num_sel_sig(S, 0.0);
  double denom_true_sig = 0.0;

  const Long64_t nent = ch.GetEntries();
  const Long64_t lim  = (max_events>0 && max_events<nent) ? max_events : nent;

  for (Long64_t i = 0; i < lim; ++i) {
    ch.GetEntry(i);

    const double w = w0(w_spline) * w0(w_genie) * w0(w_nxsec);

    if (mc_signal) denom_true_sig += w;

    bool pass = true;
    for (size_t s = 0; s < S; ++s) {
      pass = pass && (*stage_vals[s]);
      if (!pass) break;
      if (mc_signal) num_sel_sig[s] += w;
    }
  }

  std::ofstream out(out_csv);
  out << "stage_index,stage_label,selected_signal_weighted,true_signal_total_weighted,integrated_efficiency\n";
  for (size_t s = 0; s < S; ++s) {
    const double eff = (denom_true_sig > 0.0) ? (num_sel_sig[s] / denom_true_sig) : 0.0;
    out << s << ",\"" << stage_cuts[s] << "\"," << num_sel_sig[s] << "," << denom_true_sig << "," << eff << "\n";
  }
  out.close();

  std::cout << "[INFO] Wrote integrated efficiency CSV: " << out_csv << std::endl;

  ch.SetBranchStatus("*", 1);
  ch.ResetBranchAddresses();
}

using BlipBucketMap   = std::map<std::string, std::unique_ptr<TH1F>>;
using GetBlipBucketFn = std::function<TH1F*(const std::string&)>;

// Per file: apply event cuts, then read a SINGLE BLIP (by index_branch) and fill var from var_branch[idx]
static void FillBlip_PerFileScaled1D(const RunGroup& rg,
                                    Sample s,
                                    const GetBlipBucketFn& get_bucket_hist,
                                    const char* index_branch,            // single int index (fallback)
                                    const char* var_branch_vec,          // e.g. "blip_energy" or blip-domain var
                                    const char* blip_true_ncategory,     // true particle that produced the blip (overlay only, blip-domain)
                                    const std::vector<std::string>& cuts,
                                    int nbins, double xmin, double xmax,
                                    Long64_t max_events = -1,
                                    bool all_blips = false,
                                    const char* indices_branch_vec = nullptr)
{
  std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(), "READ"));
  if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; return; }
  TTree* tr = (TTree*)f->Get(kTreeName);
  if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; return; }

  tr->SetBranchStatus("*", 0);

  // Event-level cuts
  const size_t ncuts = cuts.size();
  std::unique_ptr<Bool_t[]> cut_vals;
  if (ncuts > 0) {
    cut_vals.reset(new Bool_t[ncuts]);
    for (size_t i = 0; i < ncuts; ++i) {
      cut_vals[i] = kFALSE;
      tr->SetBranchStatus(cuts[i].c_str(), 1);
      tr->SetBranchAddress(cuts[i].c_str(), &cut_vals[i]);
    }
  }

  // Index sources
  int blip_idx = -1;                            // single index branch (fallback)
  if (!all_blips && !indices_branch_vec) {
    tr->SetBranchStatus(index_branch, 1);
    tr->SetBranchAddress(index_branch, &blip_idx);
  }
  std::vector<int>* indices_v = nullptr;       // per-event vector<int> of indices
  if (indices_branch_vec && tr->GetBranch(indices_branch_vec)) {
    tr->SetBranchStatus(indices_branch_vec, 1);
    tr->SetBranchAddress(indices_branch_vec, &indices_v);
  }

  // Value vector to plot (domain must match how you intend to index)
  std::vector<float>* var_v = nullptr;
  tr->SetBranchStatus(var_branch_vec, 1);
  tr->SetBranchAddress(var_branch_vec, &var_v);

  // PDG (overlay only; Blip domain). Dirt/EXT do not use PDG.
  std::vector<int>* category_v = nullptr;
  const bool is_overlay = (s == Sample::kOverlay);
  if (is_overlay && tr->GetBranch(blip_true_ncategory)) {
    tr->SetBranchStatus(blip_true_ncategory, 1);
    tr->SetBranchAddress(blip_true_ncategory, &category_v);
  }

  // MC weights (overlay/dirt)
  std::vector<double>* w_spline = nullptr;
  std::vector<double>* w_genie  = nullptr;
  std::vector<double>* w_nxsec  = nullptr;
  const bool is_mc = (s == Sample::kOverlay || s == Sample::kDirt);
  //if (is_mc) {
  //  tr->SetBranchStatus("weight_splines_general_Spline", 1);
  //  tr->SetBranchStatus("weight_TunedCentralValue_UBGenie", 1);
  //  tr->SetBranchStatus("weight_neutron_argon_xsec", 1);
  //  tr->SetBranchAddress("weight_splines_general_Spline", &w_spline);
  //  tr->SetBranchAddress("weight_TunedCentralValue_UBGenie", &w_genie);
  //  tr->SetBranchAddress("weight_neutron_argon_xsec", &w_nxsec);
  //}

  // Per-file temporaries by legend bucket
  std::map<std::string, std::unique_ptr<TH1F>> tmp_by_bucket;
  auto get_tmp = [&](const std::string& label)->TH1F*{
    auto it = tmp_by_bucket.find(label);
    if (it != tmp_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(("tmp_blip_"+label).c_str(), "", nbins, xmin, xmax);
    h->Sumw2(true); h->SetDirectory(nullptr);
    TH1F* raw = h.get();
    tmp_by_bucket[label] = std::move(h);
    return raw;
  };

  const Long64_t nent = tr->GetEntries();
  const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

  for (Long64_t i=0; i<lim; ++i) {
    tr->GetEntry(i);

    // Event-level cuts
    bool pass = true;
    for (size_t j=0; j<ncuts; ++j) { if (!cut_vals[j]) { pass=false; break; } }
    if (!pass) continue;

    if (!var_v) continue;
    const int n_vals = (int)var_v->size();

    // Build per-event index list (precedence: all_blips > indices_branch_vec > index_branch)
    std::vector<int> idx_list;
    if (all_blips) {
      idx_list.resize(n_vals);
      std::iota(idx_list.begin(), idx_list.end(), 0);
    } else if (indices_v && !indices_v->empty()) {
      idx_list = *indices_v;
    } else {
      if (blip_idx >= 0) idx_list.push_back(blip_idx);
    }

    if (idx_list.empty()) continue;

    std::sort(idx_list.begin(), idx_list.end());
    idx_list.erase(std::unique(idx_list.begin(), idx_list.end()), idx_list.end());

    // Event weight (same for all indices in this event)
    double w = 1.0;
    // weights disabled for testing — uncomment and restore guard when ready
    //if (is_mc) {
    //  if (!w_spline || !w_genie || !w_nxsec) continue;
    //  w *= (*w_spline)[0] * (*w_genie)[0];
    //  w *= (*w_nxsec)[0];
    //  if (!std::isfinite(w) || std::abs(w) > 1e6) continue;
    //}

    for (int idx : idx_list) {
      if (idx < 0 || idx >= n_vals) continue;

      // Overlay-only true_blip_ncategory. If missing/OOB, default to "Other".
      int category = 0;
      if (is_overlay) {
        if (category_v && idx < (int)category_v->size()) category = (*category_v)[idx];
        else category = 0;
      }

      const float val = (*var_v)[idx];

      std::string label;
      if (s == Sample::kOverlay)   label = BucketFrom_true_ncat(category);
      else if (s == Sample::kDirt) label = "Dirt";
      else                         label = "EXT"; 

      get_tmp(label)->Fill(val, w);
    }
  }

  // Per-file normalization to NuWro: triggers for EXT, POT for Overlay/Dirt
  const double scale = PerFileScale(s, rg.run_index, f.get());
  for (auto& kv : tmp_by_bucket) {
    kv.second->Scale(scale);
    get_bucket_hist(kv.first)->Add(kv.second.get());
  }

  tr->SetBranchStatus("*", 1);
  tr->ResetBranchAddresses();
}

static void FillBlip_Data_FromFiles1D(const std::vector<RunGroup>& nuwro_files,
                                      const char* index_branch,
                                      const char* var_branch_vec,
                                      const std::vector<std::string>& cuts,
                                      TH1F& h_out,
                                      Long64_t max_events = -1,
                                      bool all_blips = false,
                                      const char* indices_branch_vec = nullptr)
{
  h_out.Sumw2(true);

  for (const auto& rg : nuwro_files) {
    std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(), "READ"));
    if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; continue; }
    TTree* tr = (TTree*)f->Get(kTreeName);
    if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; continue; }

    tr->SetBranchStatus("*", 0);

    const size_t ncuts = cuts.size();
    std::unique_ptr<Bool_t[]> cut_vals;
    if (ncuts > 0) {
      cut_vals.reset(new Bool_t[ncuts]);
      for (size_t i = 0; i < ncuts; ++i) {
        cut_vals[i] = kFALSE;
        tr->SetBranchStatus(cuts[i].c_str(), 1);
        tr->SetBranchAddress(cuts[i].c_str(), &cut_vals[i]);
      }
    }

    int blip_idx = -1;
    if (!all_blips && !indices_branch_vec) {
      tr->SetBranchStatus(index_branch, 1);
      tr->SetBranchAddress(index_branch, &blip_idx);
    }
    std::vector<int>* indices_v = nullptr;
    if (indices_branch_vec && tr->GetBranch(indices_branch_vec)) {
      tr->SetBranchStatus(indices_branch_vec, 1);
      tr->SetBranchAddress(indices_branch_vec, &indices_v);
    }

    std::vector<float>* var_v = nullptr;
    tr->SetBranchStatus(var_branch_vec, 1);
    tr->SetBranchAddress(var_branch_vec, &var_v);

    const Long64_t nent = tr->GetEntries();
    const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

    for (Long64_t i=0; i<lim; ++i) {
      tr->GetEntry(i);

      bool pass = true;
      for (size_t j=0; j<ncuts; ++j) { if (!cut_vals[j]) { pass=false; break; } }
      if (!pass) continue;

      if (!var_v) continue;
      const int n_vals = (int)var_v->size();

      std::vector<int> idx_list;
      if (all_blips) {
        idx_list.resize(n_vals);
        std::iota(idx_list.begin(), idx_list.end(), 0);
      } else if (indices_v && !indices_v->empty()) {
        idx_list = *indices_v;
      } else {
        if (blip_idx >= 0) idx_list.push_back(blip_idx);
      }
      if (idx_list.empty()) continue;

      std::sort(idx_list.begin(), idx_list.end());
      idx_list.erase(std::unique(idx_list.begin(), idx_list.end()), idx_list.end());

      for (int idx : idx_list) {
        if (idx < 0 || idx >= n_vals) continue;
        h_out.Fill((*var_v)[idx], 1.0);
      }
    }

    tr->SetBranchStatus("*", 1);
    tr->ResetBranchAddresses();
  }
}


static void FillBlip_FromFiles1D(const std::vector<RunGroup>& files,
                                Sample s,
                                const GetBlipBucketFn& get_bucket_hist,
                                const char* index_branch,
                                const char* var_branch_vec,
                                const char* blip_true_ncategory,
                                int nbins, double xmin, double xmax,
                                const std::vector<std::string>& cuts,
                                Long64_t max_events = -1,
				bool all_blips = false,
				const char* indices_branch_vec = nullptr)
{
  for (const auto& rg : files)
    FillBlip_PerFileScaled1D(rg, s, get_bucket_hist,
                            index_branch, var_branch_vec, blip_true_ncategory,
                            cuts, nbins, xmin, xmax, max_events, all_blips, indices_branch_vec);
}

static void MakeBlipStack1D(const std::vector<RunGroup>& overlay_files,
                           const std::vector<RunGroup>& dirt_files,
                           const std::vector<RunGroup>& ext_files,
                           const std::vector<RunGroup>& data_files,  // NEW
                           const char* index_branch,                   // e.g. "blip_trkid"
                           const char* var_branch_vec,                 // e.g. "blip_energy"
                           const char* blip_true_ncategory,            // e.g. "blip_true_ncategory"
                           const char* axis_title,
                           int nbins, double xmin, double xmax,
                           const std::vector<std::string>& cuts,
                           Long64_t max_events = -1,
			   bool all_blips = false,
			   const char* indices_branch_vec = nullptr)
{
  std::string cutlabel;
  for (size_t i=0;i<cuts.size();++i) {
    if(i) { cutlabel+=""; }
    cutlabel+=cuts[i];
  }
  const std::string base = std::string(var_branch_vec) + "_idx_" + index_branch + "_" + (cutlabel.empty()?"NoCuts":cutlabel) + "_" + (indices_branch_vec?indices_branch_vec:"");

  BlipBucketMap h_by_bucket;
  auto get_bucket_hist = [&](const std::string& label)->TH1F*{
    auto it = h_by_bucket.find(label);
    if (it != h_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(
      (std::string("h_blip_")+label+"_"+base).c_str(),
      (std::string(";")+axis_title+";Events (scaled to NuWro)").c_str(),
      nbins, xmin, xmax);
    h->Sumw2();
    h->SetDirectory(nullptr);
    const int col = (kBlipColors.count(label) ? kBlipColors.at(label) : kGray+2);
    h->SetFillColor(col);
    h->SetLineColor(col);
    if (label == "EXT") {
      h->SetFillStyle(3005);
      h->SetLineWidth(2);
    }
    TH1F* raw = h.get();
    h_by_bucket[label] = std::move(h);
    return raw;
  };

  FillBlip_FromFiles1D(overlay_files, Sample::kOverlay, get_bucket_hist,
                      index_branch, var_branch_vec, blip_true_ncategory,
                      nbins, xmin, xmax, cuts, max_events, all_blips, indices_branch_vec);
  FillBlip_FromFiles1D(dirt_files,    Sample::kDirt,    get_bucket_hist,
                      index_branch, var_branch_vec, blip_true_ncategory,
                      nbins, xmin, xmax, cuts, max_events, all_blips, indices_branch_vec);
  FillBlip_FromFiles1D(ext_files,     Sample::kEXT,     get_bucket_hist,
                      index_branch, var_branch_vec, blip_true_ncategory,
                      nbins, xmin, xmax, cuts, max_events, all_blips, indices_branch_vec);

  TH1F h_BNB_data(("h_blip_nuwro_"+base).c_str(),
               (std::string(";")+axis_title+";Events (NuWro fake data)").c_str(),
               nbins, xmin, xmax);
  h_BNB_data.SetLineColor(kBlack);
  h_BNB_data.SetLineWidth(3);
  h_BNB_data.SetMarkerStyle(20);
  h_BNB_data.SetMarkerColor(kBlack);

  FillBlip_Data_FromFiles1D(data_files, index_branch, var_branch_vec,
                             cuts, h_BNB_data, max_events, all_blips, indices_branch_vec);

  /*if (auto it = h_by_bucket.find("EXT"); it != h_by_bucket.end()) {
    std::cout << "Adding EXT distribution to nuwro histogram" << std::endl;
    h_BNB_data.Add(it->second.get());
  }

  if (auto it = h_by_bucket.find("Dirt"); it != h_by_bucket.end()) {
    h_BNB_data.Add(it->second.get());
  }*/

  THStack st(("st_blip_" + base).c_str(), (std::string(";") + axis_title + ";Events").c_str());
  
  //std::vector<std::string> order = {"Muon","Proton","Pion","Other","Dirt","EXT"};
  std::vector<std::string> order =  {"Primary (n,1p)", "Primary (n,Np)", "Dirt", "EXT", "Secondary (n,Xp)", "Primary (n,gamma)", "Secondary (n,gamma)", "ncapture gamma", "Cosmics and radiogenics", "Other (truth-matched)"};

  for (const auto& name : order) {
    auto it = h_by_bucket.find(name);
    if (it != h_by_bucket.end() && it->second->Integral() > 0)
      st.Add(it->second.get());
  }


  TCanvas c(("c_blip_" + base).c_str(), ("Blip stack: " + base).c_str(), 1100, 750);
  //c.SetLogy();
  gStyle->SetOptStat(0); gPad->SetLeftMargin(0.15);
  double max_stack = st.GetMaximum();
  double max_data  = h_BNB_data.GetMaximum();
  double ymax = 1.25 * std::max(max_stack, max_data);
  st.SetMaximum(ymax);
  //st.SetMaximum(10e7);
  //st.SetMinimum(.01);
  st.Draw("hist");
  st.GetXaxis()->SetTitle(axis_title);
  st.GetYaxis()->SetTitle("Candidates");
  //st.GetYaxis()->SetTitleOffset(1.6);
  h_BNB_data.Draw("E1 SAME");

  double total_mc = 0.0;
  for (auto& kv : h_by_bucket) total_mc += kv.second->Integral();

  TLegend leg(0.30, 0.55, 0.90, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.022);
  leg.SetNColumns(2); // helps when all 6 are present

  for (const auto& name : order) {
    auto it = h_by_bucket.find(name);
    if (it == h_by_bucket.end()) continue;
    int Nbins = it->second->GetNbinsX();
    const double val = it->second->Integral(0, Nbins+1);
    if (val <= 0) continue;
    std::string lab = name;
    if (total_mc > 0) lab += Form(" (%.1f%%)", 100.0*val/total_mc);
    leg.AddEntry(it->second.get(), lab.c_str(), "f");
  }
  leg.AddEntry(&h_BNB_data, "Open Data (1.4e20 POT)", "lep");  //"Nuwro Fake Data (1.15e20 POT)"
  leg.Draw();

  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.035);          // <-- explicit size!
  lat.SetTextAlign(13);           // left/top-ish
  lat.DrawLatexNDC(0.1, 0.98, "MicroBooNE Simulation");
  lat.SetTextSize(0.035);
  //lat.DrawLatexNDC(0.4, 0.98, "Simulation events: 5539664");
  //lat.DrawLatexNDC(0.7, 0.98, "Data Events: 4341548");

  const std::string png = base + "_blip.png";
  const std::string pdf = base + "_blip.pdf";
  c.SaveAs(png.c_str());
  c.SaveAs(pdf.c_str());
  std::cout << "[saved] " << png << " and " << pdf << "\n";
}

static void FillBuckets_PerFileScaled1D(const RunGroup& rg,
                                        Sample s,
                                        const GetBucketHistFn& get_bucket_hist,
                                        const char* var_branch,
                                        const std::vector<std::string>& cuts,
                                        int nbins, double xmin, double xmax,
                                        Long64_t max_events,
					Long64_t& n_events_out)
{
  std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(), "READ"));
  if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; return; }
  TTree* tr = (TTree*)f->Get(kTreeName);
  if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; return; }

  tr->SetBranchStatus("*", 0);

  float var_val = 0.f;
  tr->SetBranchStatus(var_branch, 1);
  tr->SetBranchAddress(var_branch, &var_val);

  int cat = 0;
  /*double ev_cat_12 = 0;
  double ev_cat_13 = 0;
  double ev_cat_14 = 0;
  double ev_cat_15 = 0;*/

  tr->SetBranchStatus("CC1muNp0piNn_EventCategory", 1);
  tr->SetBranchAddress("CC1muNp0piNn_EventCategory", &cat);

  const size_t ncuts = cuts.size();
  std::unique_ptr<Bool_t[]> cut_vals;
  if (ncuts > 0) {
    cut_vals.reset(new Bool_t[ncuts]);
    for (size_t i = 0; i < ncuts; ++i) {
      cut_vals[i] = kFALSE;
      tr->SetBranchStatus(cuts[i].c_str(), 1);
      tr->SetBranchAddress(cuts[i].c_str(), &cut_vals[i]);
    }
  }

  std::vector<double>* w_spline = nullptr;
  std::vector<double>* w_genie  = nullptr;
  std::vector<double>* w_nxsec  = nullptr;
  //if (s == Sample::kOverlay || s == Sample::kDirt) {
  //  tr->SetBranchStatus("weight_splines_general_Spline", 1);
  //  tr->SetBranchStatus("weight_TunedCentralValue_UBGenie", 1);
  //  tr->SetBranchStatus("weight_neutron_argon_xsec", 1);
  //  tr->SetBranchAddress("weight_splines_general_Spline", &w_spline);
  //  tr->SetBranchAddress("weight_TunedCentralValue_UBGenie", &w_genie);
  //  tr->SetBranchAddress("weight_neutron_argon_xsec", &w_nxsec);
  //}

  // Per-file temporaries per legend bucket
  std::map<std::string, std::unique_ptr<TH1F>> tmp_by_bucket;
  auto get_tmp = [&](const std::string& label)->TH1F*{
    auto it = tmp_by_bucket.find(label);
    if (it != tmp_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(("tmp_"+label).c_str(), "", nbins, xmin, xmax);
    h->Sumw2(true); h->SetDirectory(nullptr);
    TH1F* raw = h.get();
    tmp_by_bucket[label] = std::move(h);
    return raw;
  };

  const Long64_t nent = tr->GetEntries();
  n_events_out += nent;
  const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

  for (Long64_t i=0;i<lim;++i) {
    tr->GetEntry(i);
    bool pass = true;
    for (size_t j=0; j<ncuts; ++j) { if (!cut_vals[j]) { pass=false; break; } }
    if (!pass) continue;

    double w = 1.0;
    if (s == Sample::kOverlay || s == Sample::kDirt) {
      // weights disabled for testing — just skip the dereference entirely
      //w *= (*w_spline)[0];
      //w *= (*w_genie)[0];
      //w *= (*w_nxsec)[0];
    }

    if (!std::isfinite(w) || std::abs(w) > 1e6) {
      continue;
    }

    /*if (cat == 12) {ev_cat_12 += w;}
    else if (cat == 13) {ev_cat_13 += w;}
    else if (cat == 14) {ev_cat_14 += w;}
    else if (cat == 15) {ev_cat_15 += w;}*/

    const std::string label = ResolveLegendBucket(s, cat);
    get_tmp(label)->Fill(var_val, w);
  }

  const double scale = PerFileScale(s, rg.run_index, f.get());
  for (auto& kv : tmp_by_bucket) {
    kv.second->Scale(scale);
    get_bucket_hist(kv.first)->Add(kv.second.get());
  }

  //std::cout << " Event Category counts. 12: " << ev_cat_12 << "; 13: " << ev_cat_13 << "; 14: " << ev_cat_14 << "; 15: " << ev_cat_15 << std::endl;

  tr->SetBranchStatus("*", 1);
  tr->ResetBranchAddresses();
}

static void FillBuckets_FromFiles1D(const std::vector<RunGroup>& files,
                                    Sample s,
                                    const GetBucketHistFn& get_bucket_hist,
                                    const char* var_branch,
                                    const std::vector<std::string>& cuts,
                                    int nbins, double xmin, double xmax,
                                    Long64_t max_events,
				    Long64_t& n_events_out)
{
  for (const auto& rg : files)
    FillBuckets_PerFileScaled1D(rg, s, get_bucket_hist,
                                var_branch, cuts, nbins, xmin, xmax, max_events, n_events_out);
}

static void FillData_FromFiles1D(const std::vector<RunGroup>& files,
                                  TH1F& h_out,
                                  const char* var_branch,
                                  const std::vector<std::string>& cuts,
                                  int nbins, double xmin, double xmax,
                                  Long64_t max_events,
				  Long64_t& n_selected_out)
{
  h_out.Sumw2(true);
  n_selected_out = 0;

  for (const auto& rg : files) {
    std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(), "READ"));
    if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; continue; }
    TTree* tr = (TTree*)f->Get(kTreeName);
    if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; continue; }

    tr->SetBranchStatus("*", 0);

    float var_val = 0.f;
    tr->SetBranchStatus(var_branch, 1);
    tr->SetBranchAddress(var_branch, &var_val);

    const size_t ncuts = cuts.size();
    std::unique_ptr<Bool_t[]> cut_vals(new Bool_t[ncuts]);
    for (size_t i = 0; i < ncuts; ++i) {
      cut_vals[i] = kFALSE;
      tr->SetBranchStatus(cuts[i].c_str(), 1);
      tr->SetBranchAddress(cuts[i].c_str(), &cut_vals[i]);
    }

    const Long64_t nent = tr->GetEntries();
    n_selected_out += nent;
    const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;
    for (Long64_t i=0;i<lim;++i) {
      tr->GetEntry(i);
      bool pass = true;
      for (size_t j=0; j<ncuts; ++j) { if (!cut_vals[j]) { pass=false; break; } }
      if (!pass) continue;
      h_out.Fill(var_val, 1.0); // NuWro: no scaling, no weights
    }

    tr->SetBranchStatus("*", 1);
    tr->ResetBranchAddresses();
  }
}

static void MakeStack1DForFlag(const std::vector<RunGroup>& overlay_files,
                               const std::vector<RunGroup>& dirt_files,
                               const std::vector<RunGroup>& ext_files,
                               const std::vector<RunGroup>& data_files,
                               const char* var_branch,
                               const char* axis_title,
                               int nbins, double xmin, double xmax,
                               const std::vector<std::string>& cuts,
                               Long64_t max_events)
{
  // 0) Label for this plot
  std::string cutlabel;
  for (size_t i=0;i<cuts.size();++i){ if(i) cutlabel+="_AND_"; cutlabel+=cuts[i]; }
  const std::string base = std::string(var_branch) + "_" + (cutlabel.empty()?"NoCuts":cutlabel);

  static const std::map<std::string, int> kGroupColors = {
    {"#nu_{#mu}CCNnNp0#pi",          633},
    {"#nu_{#mu}CC0#pi nonsignal",    417},
    {"#nu_{#mu}CCOther",             920},
    {"#nu_{#mu}CCN#pi",               394},
    {"Other",                        616},
    {kLabel_Dirt,                    28},
    {kLabel_Cosmic,                  9},
    {kLabel_EXT,                     28} // hatched below
  };

  // 2) Bucket map + factory
  BucketMap h_by_bucket;
  auto get_bucket_hist = [&](const std::string& label) -> TH1F* {
    auto it = h_by_bucket.find(label);
    if (it != h_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(
      (std::string("h_")+label+"_"+base).c_str(),
      (std::string(";")+axis_title+";Events (scaled to NuWro)").c_str(),
      nbins, xmin, xmax);
    h->Sumw2(); h->SetDirectory(nullptr);
    const int col = (kGroupColors.count(label) ? kGroupColors.at(label) : kGray+2);
    h->SetFillColor(col); h->SetLineColor(col);
    if (label == kLabel_EXT) { h->SetFillStyle(3005); h->SetLineWidth(2); }
    TH1F* raw = h.get();
    h_by_bucket[label] = std::move(h);
    return raw;
  };

  // 3) Data points
  TH1F h_BNB_data(("h_BNB_data_" + base).c_str(),
               (std::string(";") + axis_title + ";Events (NuWro fake data)").c_str(),
               nbins, xmin, xmax);
  h_BNB_data.SetLineColor(kBlack);
  h_BNB_data.SetLineWidth(3);
  h_BNB_data.SetMarkerStyle(20);
  h_BNB_data.SetMarkerColor(kBlack);

  Long64_t n_overlay_events = 0;
  Long64_t n_dirt_events    = 0;
  Long64_t n_ext_events     = 0;

  // 4) Fill MC buckets (properly scaled per file)
  FillBuckets_FromFiles1D(overlay_files, Sample::kOverlay,
                          get_bucket_hist, var_branch, cuts, nbins, xmin, xmax, max_events, n_overlay_events);
  FillBuckets_FromFiles1D(dirt_files,    Sample::kDirt,
                          get_bucket_hist, var_branch, cuts, nbins, xmin, xmax, max_events, n_dirt_events);
  FillBuckets_FromFiles1D(ext_files,     Sample::kEXT,
                          get_bucket_hist, var_branch, cuts, nbins, xmin, xmax, max_events, n_ext_events);

  // 5) Fill NuWro points (no scaling)
  Long64_t n_data_events = 0;
  FillData_FromFiles1D(data_files, h_BNB_data, var_branch, cuts, nbins, xmin, xmax, max_events, n_data_events);

  /*if (auto it = h_by_bucket.find(kLabel_EXT); it != h_by_bucket.end()) {
    h_BNB_data.Add(it->second.get());  // NuWro fake data := NuWro + EXT
  }

  if (auto it = h_by_bucket.find(kLabel_Dirt); it != h_by_bucket.end()) {
    h_BNB_data.Add(it->second.get());  // NuWro fake data := NuWro + EXT
  }*/

  std::cout << "\n=== Bucket integrals (debug) ===\n";
  double total_mc_debug = 0.0;
  for (auto& kv : h_by_bucket) {
    double integral = kv.second->Integral();
    std::cout << kv.first << " : " << integral << std::endl;
    total_mc_debug += integral;
  }
  std::cout << "Total MC = " << total_mc_debug << "\n";

  // 6) Make stack
  THStack st(("st_" + base).c_str(), (std::string(";") + axis_title + ";Events").c_str());
  auto add_if = [&](const std::string& name){
    auto it = h_by_bucket.find(name);
    if (it != h_by_bucket.end() && it->second->Integral()>0) st.Add(it->second.get());
  };
  add_if("#nu_{#mu}CCNnNp0#pi");
  add_if("#nu_{#mu}CC0#pi nonsignal");
  add_if("#nu_{#mu}CCOther");
  add_if("#nu_{#mu}CCN#pi");
  add_if("Other");
  add_if(kLabel_Dirt);
  add_if(kLabel_Cosmic);
  add_if(kLabel_EXT);

  // 7) Draw
  TCanvas c(("c_" + base).c_str(), ("Stack (scaled to NuWro): " + base).c_str(), 1100, 750);
  gStyle->SetOptStat(0);
  st.Draw("hist");
  st.GetXaxis()->SetTitle(axis_title);
  st.GetYaxis()->SetTitleOffset(1.5);
  st.GetYaxis()->SetTitle("Events");
  h_BNB_data.Draw("E1 SAME");

  // 8) Legend with fractions (of total MC)
  double total_mc = 0.0;
  for (auto& kv : h_by_bucket) {
    const double val = kv.second->Integral();
    if (val > 0) total_mc += val;
  }

  TLegend leg(0.10, 0.55, 0.85, 0.85);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42); leg.SetTextSize(0.030);

  auto add_leg = [&](const std::string& name, const char* opt="f"){
    auto it = h_by_bucket.find(name);
    if (it==h_by_bucket.end()) return;
    const double val = it->second->Integral();
    if (val<=0 || total_mc <= 0) return;
    const double frac = 100.0*val/total_mc;
    std::string lab = name + Form(" (%.1f%%)", frac);
    leg.AddEntry(it->second.get(), lab.c_str(), opt);
  };
  add_leg("#nu_{#mu}CCNnNp0#pi");
  add_leg("#nu_{#mu}CC0#pi nonsignal");
  add_leg("#nu_{#mu}CCOther");
  add_leg("#nu_{#mu}CCN#pi");
  add_leg("Other");
  add_leg(kLabel_Dirt);
  add_leg(kLabel_Cosmic);
  add_leg(kLabel_EXT);
  leg.AddEntry(&h_BNB_data, "BNB data (1.15e20 POT)", "lep");
  leg.Draw();

  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextSize(0.04);          // <-- explicit size!
  lat.SetTextAlign(13);           // left/top-ish
  lat.DrawLatexNDC(0.15, 0.88, "MicroBooNE Simulation");
  lat.SetTextSize(0.035);
  lat.DrawLatexNDC(0.15, 0.83, Form("Simulation events: %lld", n_overlay_events));
  lat.DrawLatexNDC(0.15, 0.78, Form("Data events: %lld", n_data_events));

  const std::string png = base + ".png";
  const std::string pdf = base + ".pdf";
  c.SaveAs(png.c_str());
  c.SaveAs(pdf.c_str());
  std::cout << "[saved] " << png << " and " << pdf << "\n";
}

static void Make2DOverlayVertexPlot(const std::vector<RunGroup>& overlay_files,
                                    const char* x_branch, const char* y_branch,
                                    const char* x_title,  const char* y_title,
                                    int nx, double xmin, double xmax,
                                    int ny, double ymin, double ymax,
                                    const std::vector<std::string>& cuts,
                                    bool draw_boxes_xz = false,
                                    bool draw_boxes_xy = false,
                                    Long64_t max_events = -1)
{
  TChain ch(kTreeName);
  for (const auto& rg : overlay_files) ch.Add(rg.file.c_str());

  ch.SetBranchStatus("*", 0);

  float xv = 0.f, yv = 0.f;
  ch.SetBranchStatus(x_branch, 1);
  ch.SetBranchStatus(y_branch, 1);
  ch.SetBranchAddress(x_branch, &xv);
  ch.SetBranchAddress(y_branch, &yv);
  const size_t ncuts = cuts.size();
  std::unique_ptr<Bool_t[]> cut_vals;
  if (ncuts > 0) {
    cut_vals.reset(new Bool_t[ncuts]);
    for (size_t i = 0; i < ncuts; ++i) {
      cut_vals[i] = kFALSE;
      ch.SetBranchStatus(cuts[i].c_str(), 1);
      ch.SetBranchAddress(cuts[i].c_str(), &cut_vals[i]);
    }
  }
  TH2F h2("h2_overlay_vertices",
          (std::string(";") + x_title + ";" + y_title).c_str(),
          nx, xmin, xmax, ny, ymin, ymax);
  h2.Sumw2(false);
  h2.SetStats(false);

  const Long64_t nent = ch.GetEntries();
  const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

  for (Long64_t i = 0; i < lim; ++i) {
    ch.GetEntry(i);

    bool pass = true;
    for (size_t j = 0; j < ncuts; ++j) { if (!cut_vals[j]) { pass = false; break; } }
    if (!pass) continue;

    h2.Fill(xv, yv, 1.0);
  }

  std::string cutlabel;
  for (size_t i=0;i<cuts.size();++i){ if (i) cutlabel+="_AND_"; cutlabel+=cuts[i]; }
  const std::string base = std::string(x_branch) + "_vs_" + y_branch + "_" +
                           (cutlabel.empty() ? "NoCuts" : cutlabel);

  TCanvas c(("c2_overlay_"+base).c_str(), ("Overlay Vertices: "+base).c_str(), 900, 750);
  gStyle->SetOptStat(0);
  gPad->SetRightMargin(0.14);
  gPad->SetLeftMargin(0.14);
  h2.Draw("COLZ");

  if (draw_boxes_xz) DrawBoxesXZ();
  if (draw_boxes_xy) DrawBoxesXY();

  const std::string png = base + "_overlay2D.png";
  const std::string pdf = base + "_overlay2D.pdf";

  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.04);          // <-- explicit size!
  lat.SetTextAlign(13);           // left/top-ish
  lat.DrawLatexNDC(0.35, 0.5, "MicroBooNE Simulation");

  c.SaveAs(png.c_str());
  c.SaveAs(pdf.c_str());
  std::cout << "[saved] " << png << " and " << pdf << "\n";
  ch.SetBranchStatus("*", 1);
  ch.ResetBranchAddresses();
}

static void WritePurityByStageCSV_OverlayDirtEXT(
    const std::vector<RunGroup>& overlay_files,
    const std::vector<RunGroup>& dirt_files,
    const std::vector<RunGroup>& ext_files,
    const std::vector<std::string>& cuts,          // ordered, cumulative stages
    const std::string& out_csv,
    Long64_t max_events = -1)
{
  auto accumulate_file = [&](const RunGroup& rg, Sample s,
                             std::vector<double>& num_sig,
                             std::vector<double>& den_all)
  {
    std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(), "READ"));
    if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; return; }
    TTree* tr = (TTree*)f->Get(kTreeName);
    if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; return; }

    tr->SetBranchStatus("*", 0);

    int cat = 0;
    tr->SetBranchStatus("CC1muNp0piNn_EventCategory", 1);
    tr->SetBranchAddress("CC1muNp0piNn_EventCategory", &cat);

    const size_t S = cuts.size();
    std::unique_ptr<Bool_t[]> cut_vals;
    if (S > 0) {
      cut_vals.reset(new Bool_t[S]);
      for (size_t i = 0; i < S; ++i) {
        cut_vals[i] = kFALSE;
        tr->SetBranchStatus(cuts[i].c_str(), 1);
        tr->SetBranchAddress(cuts[i].c_str(), &cut_vals[i]);
      }
    }

    std::vector<double>* w_spline = nullptr;
    std::vector<double>* w_genie  = nullptr;
    std::vector<double>* w_nxsec  = nullptr;
    //if (s == Sample::kOverlay || s == Sample::kDirt) {
    //  tr->SetBranchStatus("weight_splines_general_Spline", 1);
    //  tr->SetBranchStatus("weight_TunedCentralValue_UBGenie", 1);
    //  tr->SetBranchStatus("weight_neutron_argon_xsec", 1);
    //  tr->SetBranchAddress("weight_splines_general_Spline", &w_spline);
    //  tr->SetBranchAddress("weight_TunedCentralValue_UBGenie", &w_genie);
    //  tr->SetBranchAddress("weight_neutron_argon_xsec", &w_nxsec);
    //}

    std::vector<double> file_num_sig(S, 0.0), file_den_all(S, 0.0);

    const Long64_t nent = tr->GetEntries();
    const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

    for (Long64_t i = 0; i < lim; ++i) {
      tr->GetEntry(i);

      double w = 1.0;
      if (s == Sample::kOverlay || s == Sample::kDirt) {
        if (!w_spline || w_spline->empty() || !w_genie || w_genie->empty() || !w_nxsec || w_nxsec->empty()) continue;
        w *= (*w_spline)[0];
        w *= (*w_genie)[0];
        w *= (*w_nxsec)[0];
        if (!std::isfinite(w) || std::abs(w) > 1e6) continue;
      }

      bool cumulative_pass = true;
      for (size_t st = 0; st < S; ++st) {
        cumulative_pass = cumulative_pass && cut_vals[st];
        if (!cumulative_pass) break;
        file_den_all[st] += w;
        if (s == Sample::kOverlay) {
          const std::string label = ResolveLegendBucket(s, cat);
          if (label == "signal") file_num_sig[st] += w;
        }
      }
    }

    const double scale = PerFileScale(s, rg.run_index, f.get());
    for (size_t st = 0; st < cuts.size(); ++st) {
      num_sig[st] += file_num_sig[st] * scale;
      den_all[st] += file_den_all[st] * scale;
    }

    tr->SetBranchStatus("*", 1);
    tr->ResetBranchAddresses();
  };

  const size_t S = cuts.size();
  std::vector<double> num_sig(S, 0.0), den_all(S, 0.0);

  for (const auto& rg : overlay_files) accumulate_file(rg, Sample::kOverlay, num_sig, den_all);
  for (const auto& rg : dirt_files)    accumulate_file(rg, Sample::kDirt,    num_sig, den_all);
  for (const auto& rg : ext_files)     accumulate_file(rg, Sample::kEXT,     num_sig, den_all);

  std::ofstream out(out_csv);
  out << "stage_index,stage_label,selected_signal_weighted,selected_total_weighted,purity\n";
  for (size_t st = 0; st < S; ++st) {
    const double pur = (den_all[st] > 0.0) ? (num_sig[st] / den_all[st]) : 0.0;
    out << st << ",\"" << cuts[st] << "\"," << num_sig[st] << "," << den_all[st] << "," << pur << "\n";
  }
  out.close();
  std::cout << "[INFO] Wrote purity-by-stage CSV: " << out_csv << std::endl;
}


int main() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  Long64_t max_events = -1;

/*  std::vector<RunGroup> overlay_files = {
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_overlay_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_overlay_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_overlay_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_overlay_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_overlay_POST_multisim.root", 4},
  };
*/

  std::vector<RunGroup> overlay_files = {
    //{"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root", 0},
    {"../overlay_out.root", 0},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run2_CC1muNp0piNn_overlay_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_overlay_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_overlay_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_overlay_POST_multisim.root", 4},
  };

/*  std::vector<RunGroup> dirt_files = {
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_dirt_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_dirt_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_dirt_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_dirt_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_dirt_POST_multisim.root", 4},
  };
*/

  std::vector<RunGroup> dirt_files = {
    //{"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root", 0},//MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260224.root
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run2_CC1muNp0piNn_dirt_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_dirt_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_dirt_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_dirt_POST_multisim.root", 4},
  };


/*  std::vector<RunGroup> ext_files = {
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_EXT_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_EXT_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_EXT_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_EXT_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_EXT_POST_multisim.root", 4},
  };
*/

  std::vector<RunGroup> ext_files = {
    //{"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_beamoffEXT_surprise_blipPlus_20260522.root", 0},
    //{"../EXT_out.root", 0},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run2_CC1muNp0piNn_EXT_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_EXT_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_EXT_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_EXT_POST_multisim.root", 4},
  };

/*  std::vector<RunGroup> nuwro_files = {
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2b_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 4},
  };
*/

  std::vector<RunGroup> nuwro_files = {
    //{"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root", 0},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run2b_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 4},
  };

  std::vector<RunGroup> data_files = {
    //{"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_11_BNB_opendata_surprise_blipPlus_20260522.root", 0},
    {"../openData_out.root", 0},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run2_Beam_On_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run3_Beam_On_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run4_Beam_On_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run5_Beam_On_POST_multisim.root", 4},
  };

  std::cout << "[info] Using per-file normalization to NuWro (POT for overlay/dirt; triggers for EXT)\n";

  const char* cuts[] = {
    "MC_Signal",
    "nu_mu_cc",
    "no_reco_showers",
    "muon_passed_mom_cuts",
    "muon_contained",
    "muon_quality_ok",
    "has_p_candidate",
    "passed_proton_pid_cut",
    "protons_contained",
    "lead_p_passed_mom_cuts",
    "has_secondary_proton_cand",
    "Selected"
  };

  //std::vector<std::string> cuts_vec = { "CC1muNp0piNn_nu_mu_cc", "CC1muNp0piNn_passed_proton_pid_cut", "CC1muNp0piNn_lead_p_passed_mom_cuts", "CC1muNp0piNn_protons_contained" };
  //Bool_t b = blip_touchtrk[0];
  std::vector<std::string> cuts_vec = {"Blip_n_Np_nu_mu_cc"};//{ "b", "b", "b", "b" };

  //std::vector<std::string> eff_stages = { "CC1muNp0piNn_nu_mu_cc", "CC1muNp0piNn_no_reco_showers", "CC1muNp0piNn_has_p_candidate", "CC1muNp0piNn_lead_p_passed_mom_cuts", "CC1muNp0piNn_Selected" };

  //WriteIntegratedEfficiencyCSV_OverlayOnly(overlay_files, eff_stages,
  //                                       "integrated_efficiency_by_stage.csv",
  //                                       /*max_events=*/-1);

  /*WritePurityByStageCSV_OverlayDirtEXT(
    overlay_files,
    dirt_files,
    ext_files,
    eff_stages,                                  // same ordered cumulative list
    "purity_by_stage.csv",
    -1
  );*/

  //MakeStack1DForFlag(overlay_files, dirt_files, ext_files, data_files, "reco_nu_vtx_sce_x", "reco #nu vtx SCE-corrected X [cm]", 52, 0.0, 260.0, cuts_vec, max_events);

  //MakeStack1DForFlag(overlay_files, dirt_files, ext_files, data_files, "reco_nu_vtx_sce_y", "reco #nu vtx SCE-corrected Y [cm]", 58, -116.0, 116.0, cuts_vec, max_events);

  //MakeStack1DForFlag(overlay_files, dirt_files, ext_files, data_files, "reco_nu_vtx_sce_z", "reco #nu vtx SCE-corrected Z [cm]", 52, 0.0, 1040.0, cuts_vec, max_events);

  //MakeStack1DForFlag(overlay_files, dirt_files, ext_files, data_files, "topological_score", "Neutrino Slice Topo Score", 40, 0.0, 1.0, cuts_vec, max_events);

/*
  Make2DOverlayVertexPlot(
    overlay_files,
    "reco_nu_vtx_sce_x", "reco_nu_vtx_sce_z",
    "reco #nu vtx SCE-corrected X [cm]", "reco #nu vtx SCE-corrected Z [cm]",
    52, 0.0, 260.0,   52, 0.0, 1040.0,
    {}, true, false,
    -1);
  + 1 for "reco_nu_vtx_sce_x", "reco_nu_vtx_sce_y",
*/

/*
  MakeBlipStack1D(
    overlay_files, dirt_files, ext_files, data_files,
    "CC1muNp0piNn_muon_candidate_idx",
    "trk_score_v",
    "backtracked_pdg",
    "Track/Shower Score",
    50, 0.0, 1.0,
    cuts_vec,
    -1,
    true
  );
  + 1 for trk_llr_pid_score_v & 1 for trk_len_v & 1 for trk_distance_v
*/

  ///// data files ///
/*
  MakeBlipStack1D(
    overlay_files, dirt_files, ext_files, data_files,
    "CC1muNp0piNn_lead_p_candidate_idx",   // requires this branch to be and int not a vector<int> via FillBlip_PerFileScaled1D
    "trk_energy_proton_v",
    "backtracked_pdg",
    "Proton Energy from Track (GeV)",
    100, 0.0, .433,
    cuts_vec,
    -1
  );
  + 1 for trk_llr_pid_score
*/

  MakeBlipStack1D(
    overlay_files, dirt_files, ext_files, data_files,   //nuwro_files,
    "",   // this only works because we input the last parameter, index vector, which trumps this param, index int
    "blip_energy",//"trk_len_v",
    "blip_true_ncategory",//"backtracked_pdg",
    "Reconstructed Blip Energy [MeVee]",//"Secondary Proton Candidate Track Length (cm)",
    20, 0.0, 10.,
    cuts_vec,
    -1,
    false,     // true =  all_blips — iterate over every blip
    "blip_trkid"   //nullptr = no index vec needed
  );

  MakeBlipStack1D(
    overlay_files, dirt_files, ext_files, data_files, //nuwro_files,
    "",
    "blip_x",//"trk_llr_pid_score_v",
    "blip_true_ncategory",//"backtracked_pdg",
    "Blip X Position",//"Secondary Proton Candidate LLRPID",
    52, 0.0, 260.0,
    cuts_vec,
    -1,
    false,
    "blip_trkid"//"sp_trk_ind"
  );

  std::cout << "[done] Saved PNG/PDF figures normalized run-by-run to NuWro.\n";
  return 0;
}
