#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <map>
#include <algorithm>
#include <limits>
#include <unordered_map>
#include <functional>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TParameter.h"

enum class Sample { kOverlay, kDirt, kEXT, kNuWro };

struct RunGroup {
  std::string file;
  int run_index;
};

static const char* kTreeName = "nuselection/NeutrinoSelectionFilteri";//"stv_tree";

static const std::vector<double> NUWRO_TRIGGERS    = {74298913,  33182039,  76943858,   83435996,  11060680};
static const std::vector<double> NUWRO_POT         = {3.08558e20,1.37908e20,3.153737e20,3.47069e20,4.5749e19};
static const std::vector<double> EXT_TRIGGERS      = {65744587,  153905891, 207000000,  238429356, 107466402};
static const std::vector<double> BNB_Data_TRIGGERS = {37273255,  61882791,  61375491,   73669038,  35265730};
static const std::vector<double> BNB_Data_POT      = {1.67e20,   2.61e20,   2.57e20,    3.20e20,   1.48e20};

enum class PFPBucket {
  kSignal,
  kNIProton,
  kPrimary,
  kPIProton,
  kOther,
  kDirt,
  kCosmicEM,
  kEXT
};

using BucketMap    = std::map<std::string, std::unique_ptr<TH1F>>;
using GetBucketFn  = std::function<TH1F*(const std::string&)>;

static inline int CosmicEM() { return std::numeric_limits<int>::lowest(); }

static inline const char* BucketName(PFPBucket b){
  switch (b){
    case PFPBucket::kSignal:    return "Signal protons";
    case PFPBucket::kNIProton:  return "NI protons (non-sig.)";
    case PFPBucket::kPrimary:   return "Primary particles";
    case PFPBucket::kPIProton:  return "PI protons";
    case PFPBucket::kOther:     return "Other particles";
    case PFPBucket::kDirt:      return "Dirt";
    case PFPBucket::kCosmicEM:  return "cosmic/EM";
    case PFPBucket::kEXT:       return "EXT";
  }
  return "Other particles";
}

static inline std::string tolower_copy(std::string s){
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return s;
}

static inline bool contains_nocase(const std::string& hay, const std::string& needle){
  auto h = tolower_copy(hay);
  auto n = tolower_copy(needle);
  return h.find(n) != std::string::npos;
}

static inline std::unordered_map<int,int> BuildMCIndexMap(const std::vector<int>& trkid)
{
  std::unordered_map<int,int> m;
  m.reserve(trkid.size()*2);
  for (int i=0;i<(int)trkid.size();++i) m[trkid[i]] = i;
  return m;
}

static std::pair<bool,int> IsFromFSNeutron(const std::vector<int>& trkid,
                const std::vector<int>& pdg,
                const std::vector<int>& mother,
                int mother_trkid)
{
  int cur = mother_trkid;
  int guard=0;
  while (cur!=0 && guard++<10000){
    auto it = std::find(trkid.begin(), trkid.end(), cur);
    if (it==trkid.end()) return {false,0};
    int idx = (int)std::distance(trkid.begin(), it);
    int this_pdg = pdg[idx];
    int mom      = mother[idx];

    if (this_pdg!=2112) return {false,0};
    if (this_pdg==2112 && mom==0) return {true, cur};
    cur = mom;
  }
  return {false,0};
}

static PFPBucket ClassifyTruthForProtonPFP(
    int backtracked_tid,
    const std::vector<int>&   mc_trkid,
    const std::vector<int>&   mc_pdg,
    const std::vector<int>&   mc_mother,
    const std::vector<std::string>& mc_process)
{
  // find the matched truth index
  int idx = -1;
  // Optional: speed up with a map if you call this many times per event
  // auto idx_map = BuildMCIndexMap(mc_trkid); auto it=idx_map.find(backtracked_tid);
  // if (it==idx_map.end()) return PFPBucket::kOther; else idx = it->second;
  for (int i=0; i<(int)mc_trkid.size(); ++i) if (mc_trkid[i]==backtracked_tid){ idx=i; break; }
  if (idx<0) return PFPBucket::kCosmicEM; //If there is a backtracked track ID that does not find a match in the trk id list, it is from a EM pfp

  const int abs_pdg   = std::abs(mc_pdg[idx]);
  const int mother_id = mc_mother[idx];
  const std::string& proc = (idx < (int)mc_process.size()) ? mc_process[idx] : std::string();
  const std::string  p    = tolower_copy(proc);

  // We're classifying *proton-like* PFPs by truth, but mis-ID can happen
  // Handle protons with process tags first
  if (abs_pdg == 2212){
    if (contains_nocase(p, "neutroninelastic")){
      auto born = IsFromFSNeutron(mc_trkid, mc_pdg, mc_mother, mother_id);
      return born.first ? PFPBucket::kSignal : PFPBucket::kNIProton;
    }
    if (contains_nocase(p, "protoninelastic")){
      return PFPBucket::kPIProton;
    }
    if (mother_id == 0){
      return PFPBucket::kPrimary;
    }
    return PFPBucket::kOther;
  }

  // Non-proton truth: if it's produced at the neutrino vertex (mother==0) mark as Primary,
  // else dump to Other.
  if (mother_id == 0) return PFPBucket::kPrimary;
  return PFPBucket::kOther;
}

static PFPBucket ClassifyPFPBucket(
    Sample sample,
    int backtracked_tid,
    const std::vector<int>&   mc_trkid,
    const std::vector<int>&   mc_pdg,
    const std::vector<int>&   mc_mother,
    const std::vector<std::string>& mc_process)
{
  if (sample == Sample::kDirt) return PFPBucket::kDirt;
  if (sample == Sample::kEXT)  return PFPBucket::kEXT;

  if (sample == Sample::kOverlay) {
    if (backtracked_tid == CosmicEM()) {
      return PFPBucket::kCosmicEM;
    }
    if (backtracked_tid > 0)
      return ClassifyTruthForProtonPFP(backtracked_tid, mc_trkid, mc_pdg, mc_mother, mc_process);
    return PFPBucket::kCosmicEM;
  }

  return PFPBucket::kOther;
}

static inline bool InPCV(double x,double y,double z){
  return (10.0 < x && x < 246.35) && (-106.5 < y && y < 106.5) && (10.0 < z && z < 1026.80);
}
static inline bool FullyContained(double sx,double sy,double sz,double ex,double ey,double ez){
  bool s_ok = InPCV(sx,sy,sz) && (sz<700.0 || sz>750.0);
  bool e_ok = InPCV(ex,ey,ez) && (ez<700.0 || ez>750.0);
  return s_ok && e_ok;
}
static inline double dist2(double x1,double y1,double z1,double x2,double y2,double z2){
  const double dx=x1-x2, dy=y1-y2, dz=z1-z2; return dx*dx+dy*dy+dz*dz;
}

static double ComputePFPProximity(const std::vector<float>* sx_v,
                                  const std::vector<float>* sy_v,
                                  const std::vector<float>* sz_v,
                                  const std::vector<float>* ex_v,
                                  const std::vector<float>* ey_v,
                                  const std::vector<float>* ez_v,
                                  float sec_sx,
                                  float sec_sy,
                                  float sec_sz,
                                  float sec_ex,
                                  float sec_ey,
                                  float sec_ez)
{
  if (!sx_v || !sy_v || !sz_v || !ex_v || !ey_v || !ez_v) return std::numeric_limits<double>::infinity();
  const int n = (int)sx_v->size();

  constexpr double eps  = 1e-4;          // ~0.1 mm tolerance
  constexpr double eps2 = eps*eps;

  double prox2 = std::numeric_limits<double>::infinity();

  for (int i = 0; i < n; ++i) {
    const double sx = sx_v->at(i), sy = sy_v->at(i), sz = sz_v->at(i);
    const double ex = ex_v->at(i), ey = ey_v->at(i), ez = ez_v->at(i);
    // skip invalid PFPs (sentinel -500) or missing coords
    if (sx <= -500 || sy <= -500 || sz <= -500) continue;

    // exclude the candidate itself by matching both endpoints within eps
    const bool same_start = (std::abs(sx - sec_sx) < eps) && (std::abs(sy - sec_sy) < eps) && (std::abs(sz - sec_sz) < eps);
    const bool same_end   = (std::abs(ex - sec_ex) < eps) && (std::abs(ey - sec_ey) < eps) && (std::abs(ez - sec_ez) < eps);
    if (same_start && same_end) continue;

    const double dss = dist2(sx,sy,sz, sec_sx,sec_sy,sec_sz);
    const double des = dist2(ex,ey,ez, sec_sx,sec_sy,sec_sz);
    const double dse = dist2(sx,sy,sz, sec_ex,sec_ey,sec_ez);
    const double dee = dist2(ex,ey,ez, sec_ex,sec_ey,sec_ez);

    double dmin2 = std::min(std::min(dss,des), std::min(dse,dee));
    if (dmin2 <= eps2) continue;                 // ignore numerically-identical overlaps

    if (dmin2 < prox2) prox2 = dmin2;
  }

  return (prox2 == std::numeric_limits<double>::infinity()) ? prox2 : std::sqrt(prox2);
}

#include <cmath>

static inline float ComputeRelativeDirection(float vtx_x, float vtx_y, float vtx_z,
                                             float sx,    float sy,    float sz,
                                             float ex,    float ey,    float ez)
{
  const float ax = sx - vtx_x;
  const float ay = sy - vtx_y;
  const float az = sz - vtx_z;

  const float bx = ex - sx;
  const float by = ey - sy;
  const float bz = ez - sz;

  const float a2 = ax*ax + ay*ay + az*az;
  const float b2 = bx*bx + by*by + bz*bz;
  if (a2 <= 0.0f || b2 <= 0.0f) return 0.0f;

  float cosv = (ax*bx + ay*by + az*bz) / std::sqrt(a2 * b2);

  if (cosv >  1.0f) cosv = 1.0f;
  if (cosv < -1.0f) cosv = -1.0f;
  return cosv;
}

static double PerFileScale(Sample s, int run_idx, TFile* f)
{
  if (run_idx < 0 || run_idx >= 5) return 1.0;

  if (s == Sample::kEXT) {
    const double num = NUWRO_TRIGGERS[run_idx];
    const double den = EXT_TRIGGERS[run_idx];
    if (den <= 0) return 1.0;
    return num / den;
  }

  if (s == Sample::kOverlay || s == Sample::kDirt) {
    auto* pot_param = dynamic_cast<TParameter<float>*>(f->Get("summed_pot"));
    if (!pot_param) return 1.0;
    const double pot_file = pot_param->GetVal();
    if (pot_file <= 0) return 1.0;
    return NUWRO_POT[run_idx] / pot_file;
  }

  return 1.0;
}

struct SelCfg {
  double min_track_len = 0.0;
  double max_pid_score = -0.15;
  bool   require_same_slice = true;
  double min_distance = 10.0;
  double max_distance = 200.0;
  double min_trk_score = 0.5;
  double min_prox = 21.0;
  double min_direction = 0.2;
};

static void FillCandidates_PerFile(const RunGroup& rg, Sample s,
				   const GetBucketFn& get_bucket_hist_slc,
                                   const GetBucketFn& get_bucket_hist_vtx,
                                   const GetBucketFn& get_bucket_hist_len,
                                   const GetBucketFn& get_bucket_hist_pid,
                                   const GetBucketFn& get_bucket_hist_trkswr,
                                   const GetBucketFn& get_bucket_hist_prox,
                                   const GetBucketFn& get_bucket_hist_dir,
                                   const SelCfg& cfg,
                                   Long64_t max_events = -1)
{
  std::string sample_name;
  switch (s) {
    case Sample::kOverlay: sample_name = "Overlay"; break;
    case Sample::kDirt:    sample_name = "Dirt";    break;
    case Sample::kEXT:     sample_name = "EXT";     break;
    case Sample::kNuWro:   sample_name = "NuWro";   break;
    default:               sample_name = "Unknown"; break;
  }
  std::cout << "[Processing] " << sample_name << " file, run index "
            << rg.run_index << " -> " << rg.file << std::endl;

  std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(),"READ"));
  if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; return; }
  TTree* tr = (TTree*)f->Get(kTreeName);
  if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; return; }

  const double file_scale = PerFileScale(s, rg.run_index, f.get());

  tr->SetBranchStatus("*",0);

  int slice_id = 0;
  float vx = 0, vy = 0, vz = 0;
  Bool_t numucc = kFALSE;

  std::vector<float> *sx_v=nullptr, *sy_v=nullptr, *sz_v=nullptr;
  std::vector<float> *ex_v=nullptr, *ey_v=nullptr, *ez_v=nullptr;
  std::vector<float> *len_v=nullptr, *pid_v=nullptr, *dist_v=nullptr, *trkscore_v=nullptr;
  std::vector<int>   *slc_v=nullptr, *gen_v=nullptr, *backtrack_tid_v=nullptr;

  std::vector<int>*         all_mc_trkid   = nullptr;
  std::vector<int>*         all_mc_pdg     = nullptr;
  std::vector<int>*         all_mc_mother  = nullptr;
  std::vector<std::string>* all_mc_process = nullptr;

  tr->SetBranchStatus("slice_id",1);                tr->SetBranchAddress("slice_id",&slice_id);
  tr->SetBranchStatus("reco_nu_vtx_sce_x",1);      tr->SetBranchAddress("reco_nu_vtx_sce_x",&vx);
  tr->SetBranchStatus("reco_nu_vtx_sce_y",1);      tr->SetBranchAddress("reco_nu_vtx_sce_y",&vy);
  tr->SetBranchStatus("reco_nu_vtx_sce_z",1);      tr->SetBranchAddress("reco_nu_vtx_sce_z",&vz);
  tr->SetBranchStatus("CC1muNp0piNn_nu_mu_cc",1);  tr->SetBranchAddress("CC1muNp0piNn_nu_mu_cc",&numucc);

  tr->SetBranchStatus("nonprim_trk_sce_start_x_v",1);   tr->SetBranchAddress("nonprim_trk_sce_start_x_v",&sx_v);
  tr->SetBranchStatus("nonprim_trk_sce_start_y_v",1);   tr->SetBranchAddress("nonprim_trk_sce_start_y_v",&sy_v);
  tr->SetBranchStatus("nonprim_trk_sce_start_z_v",1);   tr->SetBranchAddress("nonprim_trk_sce_start_z_v",&sz_v);
  tr->SetBranchStatus("nonprim_trk_sce_end_x_v",1);     tr->SetBranchAddress("nonprim_trk_sce_end_x_v",&ex_v);
  tr->SetBranchStatus("nonprim_trk_sce_end_y_v",1);     tr->SetBranchAddress("nonprim_trk_sce_end_y_v",&ey_v);
  tr->SetBranchStatus("nonprim_trk_sce_end_z_v",1);     tr->SetBranchAddress("nonprim_trk_sce_end_z_v",&ez_v);
  tr->SetBranchStatus("nonprim_slc_id_v",1);            tr->SetBranchAddress("nonprim_slc_id_v",&slc_v);
  tr->SetBranchStatus("nonprim_trk_len_v",1);           tr->SetBranchAddress("nonprim_trk_len_v",&len_v);
  tr->SetBranchStatus("nonprim_trk_llr_pid_score_v",1); tr->SetBranchAddress("nonprim_trk_llr_pid_score_v",&pid_v);
  tr->SetBranchStatus("nonprim_trk_distance_v",1);      tr->SetBranchAddress("nonprim_trk_distance_v",&dist_v);
  tr->SetBranchStatus("nonprim_trk_score_v",1);         tr->SetBranchAddress("nonprim_trk_score_v",&trkscore_v);
  tr->SetBranchStatus("nonprim_pfp_generation",1);     tr->SetBranchAddress("nonprim_pfp_generation",&gen_v);
  tr->SetBranchStatus("nonprim_backtracked_tid",1);     tr->SetBranchAddress("nonprim_backtracked_tid",&backtrack_tid_v);

  if (s == Sample::kOverlay) {
    if (tr->GetBranch("all_mc_trkid"))   { tr->SetBranchStatus("all_mc_trkid",1);   tr->SetBranchAddress("all_mc_trkid",&all_mc_trkid); }
    if (tr->GetBranch("all_mc_pdg"))     { tr->SetBranchStatus("all_mc_pdg",1);     tr->SetBranchAddress("all_mc_pdg",&all_mc_pdg); }
    if (tr->GetBranch("all_mc_mother"))  { tr->SetBranchStatus("all_mc_mother",1);  tr->SetBranchAddress("all_mc_mother",&all_mc_mother); }
    if (tr->GetBranch("all_mc_process")) { tr->SetBranchStatus("all_mc_process",1); tr->SetBranchAddress("all_mc_process",&all_mc_process); }
  }

  std::vector<double>* w_spline=nullptr; std::vector<double>* w_genie=nullptr; std::vector<double>* w_nxsec=nullptr;
  const bool is_mc = (s==Sample::kOverlay || s==Sample::kDirt);
  if (is_mc) {
    tr->SetBranchStatus("weight_splines_general_Spline",1);
    tr->SetBranchStatus("weight_TunedCentralValue_UBGenie",1);
    tr->SetBranchStatus("weight_neutron_argon_xsec",1);
    tr->SetBranchAddress("weight_splines_general_Spline",&w_spline);
    tr->SetBranchAddress("weight_TunedCentralValue_UBGenie",&w_genie);
    tr->SetBranchAddress("weight_neutron_argon_xsec",&w_nxsec);
  }

  const Long64_t nent = tr->GetEntries();
  const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

  for (Long64_t i=0;i<lim;++i) {
    tr->GetEntry(i);
    if (!numucc) continue;
    if (!sx_v || !len_v || !pid_v || !slc_v ) continue;
    const int np = (int)sx_v->size();

    double w = 1.0;
    if (is_mc) {
      if (!w_spline || !w_genie || !w_nxsec) continue;
      //w *= (*w_spline)[0] * (*w_genie)[0] * (*w_nxsec)[0];
      if (!std::isfinite(w) || std::abs(w)>1e6) continue;
    }
    w *= file_scale;

    for (int ip=0; ip<np; ++ip) {

      const double sx=sx_v->at(ip), sy=sy_v->at(ip), sz=sz_v->at(ip);
      const double ex=ex_v->at(ip), ey=ey_v->at(ip), ez=ez_v->at(ip);
      const bool in_slice = (slc_v && ip < (int)slc_v->size() && slc_v->at(ip) == slice_id);
      const int bt_tid_slc = (backtrack_tid_v && ip < (int)backtrack_tid_v->size())
                       ? backtrack_tid_v->at(ip) : 0;
      const PFPBucket b_slc = ClassifyPFPBucket(
          s, bt_tid_slc, *all_mc_trkid, *all_mc_pdg, *all_mc_mother, *all_mc_process);

      if (!FullyContained(sx,sy,sz,ex,ey,ez)) continue;

      const int gen = gen_v->at(ip);
      if (gen != 2) continue;

      if (get_bucket_hist_slc) {
        get_bucket_hist_slc(BucketName(b_slc))->Fill(in_slice ? 2.0 : 1.0, w);
      }

      if (cfg.require_same_slice && slc_v->at(ip)!=slice_id) continue;

      const double L = len_v->at(ip);
      const double trk_score = trkscore_v->at(ip);
      const double pid = pid_v->at(ip);
      //const double dist = dist_v->at(ip);
      const double dist = std::sqrt(dist2(sx,sy,sz,vx,vy,vz));
      const double prox = ComputePFPProximity(sx_v,sy_v,sz_v,ex_v,ey_v,ez_v, sx, sy, sz, ex, ey, ez);
      const float dir = ComputeRelativeDirection(vx,vy,vz, sx,sy,sz, ex,ey,ez);

      if (L < cfg.min_track_len) continue;
      if (trk_score < cfg.min_trk_score) continue;
      if (pid > cfg.max_pid_score) continue;
      if (dist < cfg.min_distance) continue;
      if (dist > cfg.max_distance) continue;
      if (prox < cfg.min_prox) continue;
      if (dir < cfg.min_direction) continue;

      const int bt_tid = (backtrack_tid_v && ip < (int)backtrack_tid_v->size())
                          ? backtrack_tid_v->at(ip) : 0;

      const PFPBucket b = (s==Sample::kOverlay)
        ? ClassifyPFPBucket(s, bt_tid,
                            (all_mc_trkid?*all_mc_trkid:std::vector<int>{}),
                            (all_mc_pdg?*all_mc_pdg:std::vector<int>{}),
                            (all_mc_mother?*all_mc_mother:std::vector<int>{}),
                            (all_mc_process?*all_mc_process:std::vector<std::string>{}))
        : ClassifyPFPBucket(s, bt_tid,
                            std::vector<int>{}, std::vector<int>{},
                            std::vector<int>{}, std::vector<std::string>{});

      TH1F* hvtx = get_bucket_hist_vtx(BucketName(b));
      TH1F* hlen = get_bucket_hist_len(BucketName(b));
      hvtx->Fill(dist, w);
      hlen->Fill(L,    w);
      get_bucket_hist_pid(BucketName(b))->Fill(pid, w);
      get_bucket_hist_trkswr(BucketName(b))->Fill(trk_score, w);
      get_bucket_hist_prox(BucketName(b))->Fill(prox, w);
      get_bucket_hist_dir(BucketName(b))->Fill(dir, w);
    }
  }

  tr->SetBranchStatus("*",1);
  tr->ResetBranchAddresses();
}


static void MakeSecondaryProtonPlots(const std::vector<RunGroup>& overlay_files,
                                     const std::vector<RunGroup>& dirt_files,
                                     const std::vector<RunGroup>& ext_files,
                                     const std::vector<RunGroup>& nuwro_files,
                                     const SelCfg& cfg,
                                     const std::string& out_prefix,
                                     bool add_ext_dirt_to_nuwro = false,
                                     Long64_t max_events = -1)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  BucketMap slc_by_bucket;
  BucketMap vtx_by_bucket;
  BucketMap len_by_bucket;
  BucketMap pid_by_bucket;
  BucketMap trkswr_by_bucket;
  BucketMap prox_by_bucket;
  BucketMap dir_by_bucket;

  static const std::map<std::string,int> kBucketColors = {
    {"Signal protons",           633},
    {"NI protons (non-sig.)",  417},
    {"Primary particles",        920},
    {"PI protons",               394},
    {"Other particles",          616},
    {"Dirt",                     28},
    {"cosmic/EM",                9},
    {"EXT",                      kBlack}
  };

  auto get_bucket_hist_slc = [&](const std::string& name)->TH1F*{
    auto it = slc_by_bucket.find(name);
    if (it != slc_by_bucket.end()) return it->second.get();

    auto h = std::make_unique<TH1F>(("h_slc_"+name).c_str(),
                                  ";PFP presence relative to #nu slice;Candidates",
                                  2, 0.5, 2.5);
    h->Sumw2(); h->SetDirectory(nullptr);
    int col = kBucketColors.count(name) ? kBucketColors.at(name) : kGray+2;
    h->SetFillColor(col); h->SetLineColor(col);
    if (name=="EXT") { h->SetFillStyle(3005); h->SetLineWidth(2); }
    h->GetXaxis()->SetBinLabel(1, "Slice #nu Not In");
    h->GetXaxis()->SetBinLabel(2, "Slice #nu In");

    TH1F* raw = h.get();
    slc_by_bucket[name] = std::move(h);
    return raw;
  };

  auto get_bucket_hist_vtx = [&](const std::string& name)->TH1F*{
    auto it = vtx_by_bucket.find(name);
    if (it != vtx_by_bucket.end()) return it->second.get();
    const std::string hname = "h_vtx_" + name;
    auto h = std::make_unique<TH1F>(hname.c_str(),
                                  ";Vertex #nu Candidate Distance to Reco Vertex (cm);Candidates",
                                  /*nbins*/ 22, /*xmin*/ 0.0, /*xmax*/ 220.0);
    h->Sumw2(); h->SetDirectory(nullptr);
    int col = kBucketColors.count(name) ? kBucketColors.at(name) : kGray+2;
    h->SetFillColor(col); h->SetLineColor(col);
    if (name=="EXT") { h->SetFillStyle(3005); h->SetLineWidth(2); }
    TH1F* raw = h.get();
    vtx_by_bucket[name] = std::move(h);
    return raw;
  };

  auto get_bucket_hist_len = [&](const std::string& name)->TH1F*{
    auto it = len_by_bucket.find(name);
    if (it != len_by_bucket.end()) return it->second.get();
    const std::string hname = "h_len_" + name;
    auto h = std::make_unique<TH1F>(hname.c_str(),
                                  ";Candidate Track Length (cm);Candidates",
                                  /*nbins*/ 14, /*xmin*/ 0.0, /*xmax*/ 70.0);
    h->Sumw2(); h->SetDirectory(nullptr);
    int col = kBucketColors.count(name) ? kBucketColors.at(name) : kGray+2;
    h->SetFillColor(col); h->SetLineColor(col);
    if (name=="EXT") { h->SetFillStyle(3005); h->SetLineWidth(2); }
    TH1F* raw = h.get();
    len_by_bucket[name] = std::move(h);
    return raw;
  };

  auto get_bucket_hist_pid = [&](const std::string& name)->TH1F*{
    auto it = pid_by_bucket.find(name);
    if (it != pid_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(("h_pid_"+name).c_str(),
                                  ";Candidate LLR PID Score;Candidates",
                                  36,-1.0,1.0);
    h->Sumw2(); h->SetDirectory(nullptr);
    int col = kBucketColors.count(name) ? kBucketColors.at(name) : kGray+2;
    h->SetFillColor(col); h->SetLineColor(col);
    if (name=="EXT") { h->SetFillStyle(3005); h->SetLineWidth(2); }
    TH1F* raw = h.get(); pid_by_bucket[name] = std::move(h); return raw;
  };

  auto get_bucket_hist_trkswr = [&](const std::string& name)->TH1F*{
    auto it = trkswr_by_bucket.find(name);
    if (it != trkswr_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(("h_trkswr_"+name).c_str(),
                                  ";Candidate Track/Shower Score;Candidates",
                                  20,0.0,1.0);
    h->Sumw2(); h->SetDirectory(nullptr);
    int col = kBucketColors.count(name) ? kBucketColors.at(name) : kGray+2;
    h->SetFillColor(col); h->SetLineColor(col);
    if (name=="EXT") { h->SetFillStyle(3005); h->SetLineWidth(2); }
    TH1F* raw = h.get(); trkswr_by_bucket[name] = std::move(h); return raw;
  };

  auto get_bucket_hist_prox = [&](const std::string& name)->TH1F*{
    auto it = prox_by_bucket.find(name);
    if (it != prox_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(("h_prox_"+name).c_str(),
                                  ";Candidate Proximity to Nearest PFP (cm);Candidates",
                                  20,0.0,50.0);
    h->Sumw2(); h->SetDirectory(nullptr);
    int col = kBucketColors.count(name) ? kBucketColors.at(name) : kGray+2;
    h->SetFillColor(col); h->SetLineColor(col);
    if (name=="EXT") { h->SetFillStyle(3005); h->SetLineWidth(2); }
    TH1F* raw = h.get(); prox_by_bucket[name] = std::move(h); return raw;
  };

  auto get_bucket_hist_dir = [&](const std::string& name)->TH1F*{
    auto it = dir_by_bucket.find(name);
    if (it != dir_by_bucket.end()) return it->second.get();
    auto h = std::make_unique<TH1F>(("h_dir_"+name).c_str(),
                                  ";cos(#chi);Candidates",
                                  40,-1.0,1.0);
    h->Sumw2(); h->SetDirectory(nullptr);
    int col = kBucketColors.count(name) ? kBucketColors.at(name) : kGray+2;
    h->SetFillColor(col); h->SetLineColor(col);
    if (name=="EXT") { h->SetFillStyle(3005); h->SetLineWidth(2); }
    TH1F* raw = h.get(); dir_by_bucket[name] = std::move(h); return raw;
  };

  TH1F h_nuw_slc("h_nuw_slc",";PFP presence relative to #nu slice;NuWro fake data", 2, 0.5, 2.5);
  TH1F h_nuw_vtx("h_nuw_vtx",";Vertex #nu Candidate Distance to Reco Vertex (cm);NuWro fake data", 22, 0.0, 220.0);
  TH1F h_nuw_len("h_nuw_len",";Candidate Track Length (cm);NuWro fake data", 14, 0.0, 70.0);
  TH1F h_nuw_pid("h_nuw_pid",";Candidate LLR PID Score; Candidates",36,-1.0,1.0);
  TH1F h_nuw_trkswr("h_nuw_trkswr",";Candidate Track/Shower Score; Candidates",20,0.0,1.0);
  TH1F h_nuw_prox("h_nuw_prox",";Candidate Proximity to Nearest PFP (cm); Number of Candidates",20,0.0,50.0);
  TH1F h_nuw_dir("h_nuw_dir",";cos(#chi);Number of Candidates",40,-1.0,1.0);

  h_nuw_slc.Sumw2(); h_nuw_vtx.Sumw2(); h_nuw_len.Sumw2(); h_nuw_pid.Sumw2(); h_nuw_trkswr.Sumw2(); h_nuw_prox.Sumw2(); h_nuw_dir.Sumw2();
  h_nuw_slc.SetDirectory(nullptr); h_nuw_vtx.SetDirectory(nullptr); h_nuw_len.SetDirectory(nullptr); h_nuw_pid.SetDirectory(nullptr); h_nuw_trkswr.SetDirectory(nullptr); h_nuw_prox.SetDirectory(nullptr); h_nuw_dir.SetDirectory(nullptr);
  h_nuw_slc.SetMarkerStyle(20); h_nuw_slc.SetLineColor(kBlack);
  h_nuw_slc.GetXaxis()->SetBinLabel(1, "Slice #nu Not In");
  h_nuw_slc.GetXaxis()->SetBinLabel(2, "Slice #nu In");
  h_nuw_vtx.SetMarkerStyle(20); h_nuw_vtx.SetLineColor(kBlack);
  h_nuw_len.SetMarkerStyle(20); h_nuw_len.SetLineColor(kBlack);
  h_nuw_pid.SetMarkerStyle(20); h_nuw_pid.SetLineColor(kBlack);
  h_nuw_trkswr.SetMarkerStyle(20); h_nuw_trkswr.SetLineColor(kBlack);
  h_nuw_prox.SetMarkerStyle(20); h_nuw_prox.SetLineColor(kBlack);
  h_nuw_dir.SetMarkerStyle(20); h_nuw_dir.SetLineColor(kBlack);

  for (const auto& rg : overlay_files) FillCandidates_PerFile(rg, Sample::kOverlay, get_bucket_hist_slc, get_bucket_hist_vtx, get_bucket_hist_len, get_bucket_hist_pid, get_bucket_hist_trkswr, get_bucket_hist_prox, get_bucket_hist_dir, cfg, max_events);
  for (const auto& rg : dirt_files) FillCandidates_PerFile(rg, Sample::kDirt,  get_bucket_hist_slc, get_bucket_hist_vtx, get_bucket_hist_len, get_bucket_hist_pid, get_bucket_hist_trkswr, get_bucket_hist_prox, get_bucket_hist_dir, cfg, max_events);
  for (const auto& rg : ext_files)  FillCandidates_PerFile(rg, Sample::kEXT,   get_bucket_hist_slc, get_bucket_hist_vtx, get_bucket_hist_len, get_bucket_hist_pid, get_bucket_hist_trkswr, get_bucket_hist_prox, get_bucket_hist_dir, cfg, max_events);

  for (const auto& rg : nuwro_files) {
    std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(),"READ"));
    if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; continue; }
    TTree* tr = (TTree*)f->Get(kTreeName);
    if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; continue; }
    tr->SetBranchStatus("*",0);
    Bool_t numucc = kFALSE;
    int slice_id = 0;
    float vx=0,vy=0,vz=0;
    std::vector<float> *sx_v=nullptr,*sy_v=nullptr,*sz_v=nullptr,*ex_v=nullptr,*ey_v=nullptr,*ez_v=nullptr,*len_v=nullptr,*pid_v=nullptr,*dist_v=nullptr,*trkscore_v=nullptr;
    std::vector<int> *slc_v=nullptr,*gen_v=nullptr;

    tr->SetBranchStatus("slice_id",1);                tr->SetBranchAddress("slice_id",&slice_id);
    tr->SetBranchStatus("CC1muNp0piNn_nu_mu_cc",1);	  tr->SetBranchAddress("CC1muNp0piNn_nu_mu_cc",&numucc);
    tr->SetBranchStatus("reco_nu_vtx_sce_x",1); 	  tr->SetBranchAddress("reco_nu_vtx_sce_x",&vx);
    tr->SetBranchStatus("reco_nu_vtx_sce_y",1); 	  tr->SetBranchAddress("reco_nu_vtx_sce_y",&vy);
    tr->SetBranchStatus("reco_nu_vtx_sce_z",1); 	  tr->SetBranchAddress("reco_nu_vtx_sce_z",&vz);
    tr->SetBranchStatus("nonprim_trk_sce_start_x_v",1);   tr->SetBranchAddress("nonprim_trk_sce_start_x_v",&sx_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_y_v",1);   tr->SetBranchAddress("nonprim_trk_sce_start_y_v",&sy_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_z_v",1);   tr->SetBranchAddress("nonprim_trk_sce_start_z_v",&sz_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_x_v",1);     tr->SetBranchAddress("nonprim_trk_sce_end_x_v",&ex_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_y_v",1);     tr->SetBranchAddress("nonprim_trk_sce_end_y_v",&ey_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_z_v",1);     tr->SetBranchAddress("nonprim_trk_sce_end_z_v",&ez_v);
    tr->SetBranchStatus("nonprim_slc_id_v",1);            tr->SetBranchAddress("nonprim_slc_id_v",&slc_v);
    tr->SetBranchStatus("nonprim_trk_len_v",1);           tr->SetBranchAddress("nonprim_trk_len_v",&len_v);
    tr->SetBranchStatus("nonprim_trk_llr_pid_score_v",1); tr->SetBranchAddress("nonprim_trk_llr_pid_score_v",&pid_v);
    tr->SetBranchStatus("nonprim_trk_distance_v",1);      tr->SetBranchAddress("nonprim_trk_distance_v",&dist_v);
    tr->SetBranchStatus("nonprim_trk_score_v",1);         tr->SetBranchAddress("nonprim_trk_score_v",&trkscore_v);
    tr->SetBranchStatus("nonprim_pfp_generation",1);      tr->SetBranchAddress("nonprim_pfp_generation",&gen_v);

    const Long64_t nent = tr->GetEntries();
    const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;
    for (Long64_t i=0;i<lim;++i){
      tr->GetEntry(i);
      if (!numucc) continue;
      if (!sx_v || !len_v || !pid_v || !slc_v ) continue;
      const int np = (int)sx_v->size();
      for (int ip=0; ip<np; ++ip){

	const double sx=sx_v->at(ip), sy=sy_v->at(ip), sz=sz_v->at(ip);
        const double ex=ex_v->at(ip), ey=ey_v->at(ip), ez=ez_v->at(ip);
	if (!FullyContained(sx,sy,sz,ex,ey,ez)) continue;

        const int gen = gen_v->at(ip);
        if (gen != 2) continue;

	const bool in_slice = (slc_v && ip < (int)slc_v->size() && slc_v->at(ip) == slice_id);
	h_nuw_slc.Fill(in_slice ? 2.0 : 1.0, 1.0);

        if (cfg.require_same_slice && slc_v->at(ip)!=slice_id) continue;

        const double L=len_v->at(ip);
        const double pid=pid_v->at(ip);
        const double trkscore=trkscore_v->at(ip);
        //const double dist=dist_v->at(ip);
	const double dist = std::sqrt(dist2(sx,sy,sz,vx,vy,vz));
        const double prox = ComputePFPProximity(sx_v,sy_v,sz_v,ex_v,ey_v,ez_v, sx, sy, sz, ex, ey, ez);
	const float dir = ComputeRelativeDirection(vx,vy,vz, sx,sy,sz, ex,ey,ez);
        if (L < cfg.min_track_len) continue;
        if (pid > cfg.max_pid_score) continue;
        if (trkscore < cfg.min_trk_score) continue;
	if (dist < cfg.min_distance) continue;
	if (dist > cfg.max_distance) continue;
        if (prox < cfg.min_prox) continue;
        if (dir < cfg.min_direction) continue;

        h_nuw_vtx.Fill(dist, 1.0);
        h_nuw_len.Fill(L, 1.0);
        h_nuw_pid.Fill(pid, 1.0);
        h_nuw_trkswr.Fill(trkscore,1.0);
        h_nuw_prox.Fill(prox,1.0);
        h_nuw_dir.Fill(dir,1.0);
      }
    }
    tr->SetBranchStatus("*",1);
    tr->ResetBranchAddresses();
  }

  /*if (add_ext_dirt_to_nuwro) {
    if (auto it = slc_by_bucket.find("EXT");  it != slc_by_bucket.end())  h_nuw_slc.Add(it->second.get());
    if (auto it = slc_by_bucket.find("Dirt"); it != slc_by_bucket.end())  h_nuw_slc.Add(it->second.get());
    if (auto it = vtx_by_bucket.find("EXT");  it != vtx_by_bucket.end())  h_nuw_vtx.Add(it->second.get());
    if (auto it = vtx_by_bucket.find("Dirt"); it != vtx_by_bucket.end())  h_nuw_vtx.Add(it->second.get());
    if (auto it = len_by_bucket.find("EXT");  it != len_by_bucket.end())  h_nuw_len.Add(it->second.get());
    if (auto it = len_by_bucket.find("Dirt"); it != len_by_bucket.end())  h_nuw_len.Add(it->second.get());
    if (auto it = pid_by_bucket.find("EXT");  it != pid_by_bucket.end())  h_nuw_pid.Add(it->second.get());
    if (auto it = pid_by_bucket.find("Dirt"); it != pid_by_bucket.end())  h_nuw_pid.Add(it->second.get());
    if (auto it = trkswr_by_bucket.find("EXT");  it != trkswr_by_bucket.end())  h_nuw_trkswr.Add(it->second.get());
    if (auto it = trkswr_by_bucket.find("Dirt"); it != trkswr_by_bucket.end())  h_nuw_trkswr.Add(it->second.get());
    if (auto it = prox_by_bucket.find("EXT");  it != prox_by_bucket.end())  h_nuw_prox.Add(it->second.get());
    if (auto it = prox_by_bucket.find("Dirt"); it != prox_by_bucket.end())  h_nuw_prox.Add(it->second.get());
    if (auto it = dir_by_bucket.find("EXT");  it != dir_by_bucket.end())  h_nuw_dir.Add(it->second.get());
    if (auto it = dir_by_bucket.find("Dirt"); it != dir_by_bucket.end())  h_nuw_dir.Add(it->second.get());
  }*/

  std::vector<std::string> order = {
    "Signal protons",
    "NI protons (non-sig.)",
    "Primary particles",
    "PI protons",
    "Other particles",
    "Dirt",
    "cosmic/EM",
    "EXT"
  };

  THStack st_vtx("st_vtx",";Candidate Distance to Reco Vertex (cm);Candidates");
  THStack st_len("st_len",";Candidate Track Length (cm);Candidates");
  THStack st_pid("st_pid",";Candidate LLR PID Score;Candidates");
  THStack st_trkswr("st_trkswr",";Candidate Track/Shower Score;Candidates");
  THStack st_prox("st_prox",";Candidate Proximity to Nearest PFP (cm);Candidates");
  THStack st_dir("st_dir",";cos(#chi);Candidates");

  for (const auto& name : order) { auto it = vtx_by_bucket.find(name);  if (it != vtx_by_bucket.end() && it->second->Integral()>0) st_vtx.Add(it->second.get()); }
  for (const auto& name : order) { auto it = len_by_bucket.find(name);  if (it != len_by_bucket.end() && it->second->Integral()>0) st_len.Add(it->second.get()); }
  for (const auto& name : order) { auto it=pid_by_bucket.find(name);    if (it!=pid_by_bucket.end()    && it->second->Integral()>0) st_pid.Add(it->second.get()); }
  for (const auto& name : order) { auto it=trkswr_by_bucket.find(name); if (it!=trkswr_by_bucket.end() && it->second->Integral()>0) st_trkswr.Add(it->second.get()); }
  for (const auto& name : order) { auto it=prox_by_bucket.find(name);   if (it!=prox_by_bucket.end()   && it->second->Integral()>0) st_prox.Add(it->second.get()); }
  for (const auto& name : order) { auto it=dir_by_bucket.find(name);    if (it!=dir_by_bucket.end()    && it->second->Integral()>0) st_dir.Add(it->second.get()); }

  auto draw_with_legend = [&](THStack& st, TH1F& h_nuw, const BucketMap& bm,
                            const std::string& suffix, bool logy=true){
    TCanvas c(("c_"+suffix).c_str(),"",1100,750);
    if (logy) { c.SetLogy(); st.SetMinimum(0.01); st.SetMaximum(100000000); }
    double ymax = std::max(h_nuw.GetMaximum(), st.GetMaximum()); 
    st.SetMaximum(1.25 * ymax);
    st.Draw("hist");
    h_nuw.Draw("E1 SAME");
    TLegend leg(0.55,0.6,0.9,0.9);
    leg.SetBorderSize(1); leg.SetLineColor(kBlack);
    leg.SetFillStyle(1001); leg.SetFillColor(kWhite);
    leg.SetTextFont(42); leg.SetTextSize(0.030);
    double tot=0.0; for (auto& kv : bm) tot += kv.second->Integral();
    for (const auto& name : order){
      auto it = bm.find(name); if (it==bm.end()) continue;
      const double v = it->second->Integral(); if (v<=0) continue;
      std::string lab = name; if (tot>0) lab += Form(" (%.1f%%)", 100.0*v/tot);
      leg.AddEntry(it->second.get(), lab.c_str(), "f");
    }
    leg.AddEntry(&h_nuw, "BND Data (1.15e20 POT)", "lep");
    leg.Draw();

    TLatex lat;
    lat.SetTextFont(42);
    lat.SetTextColor(kBlack);
    lat.SetTextSize(0.035);          // <-- explicit size!
    lat.SetTextAlign(13);           // left/top-ish
    lat.DrawLatexNDC(0.1, 0.95, "MicroBooNE Simulation");
    lat.DrawLatexNDC(0.4, 0.95, "Simulation events: 5539664");
    lat.DrawLatexNDC(0.7, 0.95, "Data Events: 4341548");

    c.SaveAs((out_prefix + "_" + suffix + ".pdf").c_str());
    c.SaveAs((out_prefix + "_" + suffix + ".png").c_str());
  };

  draw_with_legend(st_vtx,   h_nuw_vtx,   vtx_by_bucket,   "vtx",   false);
  draw_with_legend(st_len,   h_nuw_len,   len_by_bucket,   "len",   false);
  draw_with_legend(st_pid,   h_nuw_pid,   pid_by_bucket,   "pid",   false);
  draw_with_legend(st_trkswr,h_nuw_trkswr,trkswr_by_bucket,"trkswr",false);
  draw_with_legend(st_prox,  h_nuw_prox,  prox_by_bucket,  "prox",  false);
  draw_with_legend(st_dir,   h_nuw_dir,   dir_by_bucket,   "dir",   false);

  THStack st_slc("st_slc",";PFP presence relative to #nu slice;Candidates");
  for (const auto& name : order) {
    auto it = slc_by_bucket.find(name);
    if (it != slc_by_bucket.end() && it->second->Integral()>0) st_slc.Add(it->second.get());
  }

  TCanvas cslc("cslc","",1100,750);
  cslc.SetLogy();
  gStyle->SetOptStat(0);
  st_slc.SetMinimum(.01);
  st_slc.Draw("hist");
  h_nuw_slc.Draw("E1 SAME");

  st_slc.GetXaxis()->SetBinLabel(1, "Slice #nu Not In");
  st_slc.GetXaxis()->SetBinLabel(2, "Slice #nu In");

  TLegend legS(0.6,0.15,0.90,0.45);
  legS.SetBorderSize(1); legS.SetLineColor(kBlack);
  legS.SetFillStyle(1001); legS.SetFillColor(kWhite);
  legS.SetTextFont(42); legS.SetTextSize(0.030);
  double total_slc = 0.0; for (auto& kv : slc_by_bucket) total_slc += kv.second->Integral();
  for (const auto& name : order) {
    auto it = slc_by_bucket.find(name); if (it==slc_by_bucket.end()) continue;
    const double val = it->second->Integral(); if (val<=0) continue;
    std::string lab = name; if (total_slc>0) lab += Form(" (%.1f%%)", 100.0*val/total_slc);
    legS.AddEntry(it->second.get(), lab.c_str(), "f");
  }
  legS.AddEntry(&h_nuw_slc, "BND Data (1.15e20 POT)", "lep");
  legS.Draw();

  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.035);          // <-- explicit size!
  lat.SetTextAlign(13);           // left/top-ish
  lat.DrawLatexNDC(0.1, 0.95, "MicroBooNE Simulation");
  lat.DrawLatexNDC(0.4, 0.95, "Simulation events: 5539664");
  lat.DrawLatexNDC(0.7, 0.95, "Data Events: 4341548");

  cslc.SaveAs((out_prefix + "_slc.pdf").c_str());
  cslc.SaveAs((out_prefix + "_slc.png").c_str());
}

static std::vector<double> NeutronKEBins() {
  return {0.00,0.04,0.08,0.12,0.16,0.20,0.24,0.28,0.32,0.36,0.40,0.44,0.48,0.52,0.56,0.60};
}

static inline int FindMCIdxByTrkID(const std::vector<int>& trkid, int id) {
  for (int i=0;i<(int)trkid.size();++i) if (trkid[i]==id) return i;
  return -1;
}

static inline int AncestorFSNeutron(int mother_trkid,
                                    const std::vector<int>& trkid,
                                    const std::vector<int>& pdg,
                                    const std::vector<int>& mother)
{
  int cur = mother_trkid;
  int guard=0;
  while (cur!=0 && guard++<10000){
    int idx = FindMCIdxByTrkID(trkid, cur);
    if (idx<0) return 0;
    const int this_pdg = pdg[idx];
    const int mom      = mother[idx];
    if (this_pdg!=2112) return 0;
    if (this_pdg==2112 && mom==0) return trkid[idx]; // FS neutron ancestor trkid
    cur = mom;
  }
  return 0;
}

static inline bool IsNeutronInelastic(const std::string& s) {
  return s == "neutronInelastic";
}

static inline int MatchedFSNeutronForSignalProton(
  int backtracked_tid,
  const std::vector<int>&   mc_trkid,
  const std::vector<int>&   mc_pdg,
  const std::vector<int>&   mc_mother,
  const std::vector<std::string>& mc_process)
{
  if (backtracked_tid<=0) return 0;
  int pidx = FindMCIdxByTrkID(mc_trkid, backtracked_tid);
  if (pidx<0) return 0;
  const int apdg = std::abs(mc_pdg[pidx]);
  if (apdg != 2212) return 0;
  const std::string& proc = (pidx < (int)mc_process.size()) ? mc_process[pidx] : std::string();
  if (!IsNeutronInelastic(proc)) return 0;
  return AncestorFSNeutron(mc_mother[pidx], mc_trkid, mc_pdg, mc_mother);
}

struct ProtonMatch { int fsn_tid; double p_ke; };

static inline ProtonMatch MatchSignalProton_WithKE(
  int backtracked_tid,
  const std::vector<int>&   mc_trkid,
  const std::vector<int>&   mc_pdg,
  const std::vector<int>&   mc_mother,
  const std::vector<float>& mc_E,
  const std::vector<std::string>& mc_process)
{
  if (backtracked_tid<=0) return {0, 0.0};
  const int pidx = FindMCIdxByTrkID(mc_trkid, backtracked_tid);
  if (pidx<0) return {0, 0.0};
  if (std::abs(mc_pdg[pidx]) != 2212) return {0, 0.0};
  const std::string& proc = (pidx < (int)mc_process.size()) ? mc_process[pidx] : std::string();
  if (!IsNeutronInelastic(proc)) return {0, 0.0};
  const int fsn = AncestorFSNeutron(mc_mother[pidx], mc_trkid, mc_pdg, mc_mother);
  if (fsn==0) return {0, 0.0};
  const double p_ke = std::max(0.0, (double)mc_E[pidx] - 0.938272);
  return {fsn, p_ke};
}

struct FSNInfo {
  double ke;
  bool has_p50;              // at least one secondary proton with KE>=0.05
};

static inline void CollectTruth_FSN_and_SecProtons(
  const std::vector<int>&   trkid,
  const std::vector<int>&   pdg,
  const std::vector<int>&   mother,
  const std::vector<float>& E,
  const std::vector<std::string>& process,
  std::unordered_map<int,FSNInfo>& fsn_map,
  std::vector<std::pair<int,double>>& sec_protons)
{
  fsn_map.clear(); sec_protons.clear();

  for (int i=0;i<(int)trkid.size();++i){
    if (pdg[i]==2112 && mother[i]==0) {
      const double nke = std::max(0.0, (double)E[i] - 0.93957);
      fsn_map[trkid[i]] = {nke, false};
    }
  }
  for (int i=0;i<(int)trkid.size();++i){
    if (std::abs(pdg[i]) == 2212 && IsNeutronInelastic(process[i])) {
      const int fsn = AncestorFSNeutron(mother[i], trkid, pdg, mother);
      if (fsn != 0) {
        const double pke = std::max(0.0, (double)E[i] - 0.938272);
        sec_protons.emplace_back(fsn, pke);
        auto it = fsn_map.find(fsn);
        if (it != fsn_map.end() && pke >= 0.05) it->second.has_p50 = true;
      }
    }
  }
}

static void MakeNeutronEfficiencyPlots(const std::vector<RunGroup>& overlay_files,
                                       const SelCfg& cfg,
                                       const std::string& out_prefix,
                                       Long64_t max_events = -1)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  // === Integrated efficiency counters ===
  double N_true_fsn = 0.0;      // all FS neutrons (denominator)
  double N_sel_prx  = 0.0;      // proximity + direction (final)

  const auto edges = NeutronKEBins();
  const int nb = (int)edges.size()-1;

  TH1D h_den("h_den",";True Neutron Kinetic Energy (GeV);Candidates", nb, edges.data());
  TH1D h_num_pre  ("h_num_pre",  ";True Neutron Kinetic Energy (GeV);Tagged", nb, edges.data());
  TH1D h_num_score("h_num_score",";True Neutron Kinetic Energy (GeV);Tagged", nb, edges.data());
  TH1D h_num_prox ("h_num_prox", ";True Neutron Kinetic Energy (GeV);Tagged", nb, edges.data());
  h_den.Sumw2(); h_num_pre.Sumw2(); h_num_score.Sumw2(); h_num_prox.Sumw2();

  for (const auto& rg : overlay_files) {
    std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(),"READ"));
    if (!f || f->IsZombie()) { std::cerr << "[err] open " << rg.file << "\n"; continue; }
    TTree* tr = (TTree*)f->Get(kTreeName);
    if (!tr) { std::cerr << "[err] tree " << kTreeName << " in " << rg.file << "\n"; continue; }

    const double file_scale = PerFileScale(Sample::kOverlay, rg.run_index, f.get());

    tr->SetBranchStatus("*",0);

    Bool_t numucc = kFALSE;
    float vx=0,vy=0,vz=0;
    tr->SetBranchStatus("CC1muNp0piNn_nu_mu_cc",1); tr->SetBranchAddress("CC1muNp0piNn_nu_mu_cc",&numucc);
    tr->SetBranchStatus("reco_nu_vtx_sce_x",1); tr->SetBranchAddress("reco_nu_vtx_sce_x",&vx);
    tr->SetBranchStatus("reco_nu_vtx_sce_y",1); tr->SetBranchAddress("reco_nu_vtx_sce_y",&vy);
    tr->SetBranchStatus("reco_nu_vtx_sce_z",1); tr->SetBranchAddress("reco_nu_vtx_sce_z",&vz);

    std::vector<float> *sx_v=nullptr, *sy_v=nullptr, *sz_v=nullptr;
    std::vector<float> *ex_v=nullptr, *ey_v=nullptr, *ez_v=nullptr;
    std::vector<float> *len_v=nullptr, *pid_v=nullptr, *dist_v=nullptr, *trkscore_v=nullptr;
    std::vector<int>   *gen_v=nullptr, *bt_v=nullptr;

    tr->SetBranchStatus("nonprim_trk_sce_start_x_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_x_v",&sx_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_y_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_y_v",&sy_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_z_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_z_v",&sz_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_x_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_x_v",&ex_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_y_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_y_v",&ey_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_z_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_z_v",&ez_v);

    tr->SetBranchStatus("nonprim_trk_len_v",1);           tr->SetBranchAddress("nonprim_trk_len_v",&len_v);
    tr->SetBranchStatus("nonprim_trk_llr_pid_score_v",1); tr->SetBranchAddress("nonprim_trk_llr_pid_score_v",&pid_v);
    tr->SetBranchStatus("nonprim_trk_distance_v",1);      tr->SetBranchAddress("nonprim_trk_distance_v",&dist_v);
    tr->SetBranchStatus("nonprim_trk_score_v",1);         tr->SetBranchAddress("nonprim_trk_score_v",&trkscore_v);
    tr->SetBranchStatus("nonprim_pfp_generation",1);      tr->SetBranchAddress("nonprim_pfp_generation",&gen_v);
    tr->SetBranchStatus("nonprim_backtracked_tid",1);     tr->SetBranchAddress("nonprim_backtracked_tid",&bt_v);

    std::vector<int>*           mc_trkid=nullptr;
    std::vector<int>*           mc_pdg=nullptr;
    std::vector<int>*           mc_mother=nullptr;
    std::vector<float>*         mc_E=nullptr;
    std::vector<std::string>*   mc_process=nullptr;

    tr->SetBranchStatus("all_mc_trkid",1);   tr->SetBranchAddress("all_mc_trkid",&mc_trkid);
    tr->SetBranchStatus("all_mc_pdg",1);     tr->SetBranchAddress("all_mc_pdg",&mc_pdg);
    tr->SetBranchStatus("all_mc_mother",1);  tr->SetBranchAddress("all_mc_mother",&mc_mother);
    tr->SetBranchStatus("all_mc_E",1);       tr->SetBranchAddress("all_mc_E",&mc_E);
    tr->SetBranchStatus("all_mc_process",1); tr->SetBranchAddress("all_mc_process",&mc_process);

    std::vector<double>* w_spline=nullptr; std::vector<double>* w_genie=nullptr; std::vector<double>* w_nxsec=nullptr;
    tr->SetBranchStatus("weight_splines_general_Spline",1);
    tr->SetBranchStatus("weight_TunedCentralValue_UBGenie",1);
    tr->SetBranchStatus("weight_neutron_argon_xsec",1);
    tr->SetBranchAddress("weight_splines_general_Spline",&w_spline);
    tr->SetBranchAddress("weight_TunedCentralValue_UBGenie",&w_genie);
    tr->SetBranchAddress("weight_neutron_argon_xsec",&w_nxsec);

    const Long64_t nent = tr->GetEntries();
    const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

    for (Long64_t i=0;i<lim;++i){
      tr->GetEntry(i);
      if (!numucc) continue;
      if (!mc_trkid||!mc_pdg||!mc_mother||!mc_E) continue;

      double w = 1.0;
      if (!w_spline || !w_genie || !w_nxsec) continue;
      w *= (*w_spline)[0] * (*w_genie)[0] * (*w_nxsec)[0];
      if (!std::isfinite(w) || std::abs(w)>1e6) continue;
      w *= file_scale;

      std::unordered_map<int,double> fsn_ke;
      fsn_ke.reserve(mc_trkid->size()/4+1);
      for (int j=0;j<(int)mc_trkid->size();++j){
        if (mc_pdg->at(j)==2112 && mc_mother->at(j)==0){
          const double ke = std::max(0.0, (double)mc_E->at(j) - 0.93957);
          fsn_ke[mc_trkid->at(j)] = ke;
          h_den.Fill(ke, w);
	  N_true_fsn += w;
        }
      }
      if (fsn_ke.empty()) continue;

      if (!sx_v||!ex_v||!len_v||!pid_v||!dist_v||!trkscore_v||!gen_v||!bt_v) continue;
      const int np = (int)len_v->size();
      if (np==0) continue;

      std::unordered_set<int> tagged_pre, tagged_score, tagged_prox;
      tagged_pre.reserve(8); tagged_score.reserve(8); tagged_prox.reserve(8);

      bool pass_pre = false;

      for (int ip=0; ip<np; ++ip){
        if (gen_v->at(ip)!=2) continue;

        const double sx=sx_v->at(ip), sy=sy_v->at(ip), sz=sz_v->at(ip);
        const double ex=ex_v->at(ip), ey=ey_v->at(ip), ez=ez_v->at(ip);
        if (!FullyContained(sx,sy,sz,ex,ey,ez)) continue;

	pass_pre = true;

        const double L=len_v->at(ip);
        if (L < cfg.min_track_len) continue;

        const double dist = dist_v->at(ip);
        const bool pass_vtxsep = (dist >= cfg.min_distance && dist <= cfg.max_distance);

        const double pid = pid_v->at(ip);
        const double tscore = trkscore_v->at(ip);
        const bool pass_scores = (tscore >= cfg.min_trk_score && pid <= cfg.max_pid_score);

        const double prox = ComputePFPProximity(sx_v,sy_v,sz_v, ex_v,ey_v,ez_v, sx,sy,sz,ex,ey,ez);
        const float dir   = ComputeRelativeDirection(vx,vy,vz, sx,sy,sz, ex,ey,ez);

        const bool pass_proxdir = (prox >= cfg.min_prox && dir >= cfg.min_direction);

        const int bt_tid = (ip < (int)bt_v->size()) ? bt_v->at(ip) : 0;
        const int fsn_tid = MatchedFSNeutronForSignalProton(bt_tid, *mc_trkid, *mc_pdg, *mc_mother, *mc_process);
        if (fsn_tid==0) continue;

        auto it = fsn_ke.find(fsn_tid);
        if (it==fsn_ke.end()) continue;
        const double nke = it->second;

        if (pass_pre) {
          if (!tagged_pre.count(fsn_tid)) { h_num_pre.Fill(nke, w); tagged_pre.insert(fsn_tid); }
        }
        if (pass_vtxsep && pass_scores) {
          if (!tagged_score.count(fsn_tid)) { h_num_score.Fill(nke, w); tagged_score.insert(fsn_tid); }
        }
        if (pass_vtxsep && pass_scores && pass_proxdir) {
          if (!tagged_prox.count(fsn_tid)) { h_num_prox.Fill(nke, w); tagged_prox.insert(fsn_tid);  N_sel_prx += w;}
        }
      }
    }

    tr->SetBranchStatus("*",1);
    tr->ResetBranchAddresses();
  }

  TH1D h_eff_pre(h_num_pre);   h_eff_pre.SetName("h_eff_pre");   h_eff_pre.Divide(&h_den);
  TH1D h_eff_sco(h_num_score); h_eff_sco.SetName("h_eff_sco");   h_eff_sco.Divide(&h_den);
  TH1D h_eff_prx(h_num_prox);  h_eff_prx.SetName("h_eff_prx");   h_eff_prx.Divide(&h_den);

  TCanvas c("c_eff","",900,650);
  h_eff_pre.SetMarkerStyle(20); h_eff_pre.SetLineWidth(2); h_eff_pre.SetLineColor(1); h_eff_pre.SetMarkerColor(1);
  h_eff_sco.SetMarkerStyle(20); h_eff_sco.SetLineWidth(2); h_eff_sco.SetLineColor(2); h_eff_sco.SetMarkerColor(2);
  h_eff_prx.SetMarkerStyle(20); h_eff_prx.SetLineWidth(2); h_eff_prx.SetLineColor(4); h_eff_prx.SetMarkerColor(4);


  h_eff_prx.GetYaxis()->SetTitle("Efficiency");
  h_eff_prx.SetMinimum(0.0); h_eff_prx.SetMaximum(0.45);

  h_eff_prx.Draw("E1");
  h_eff_sco.Draw("E1 SAME");
  h_eff_pre.Draw("E1 SAME");

  TLegend leg(0.15,0.7,0.53,0.9);
  leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42);
  leg.AddEntry(&h_eff_pre, "Preselection", "lep");
  leg.AddEntry(&h_eff_sco, "Scores and Vertex Separation", "lep");
  leg.AddEntry(&h_eff_prx, "Proximity and Direction", "lep");
  leg.Draw();

  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.04);          // <-- explicit size!
  lat.SetTextAlign(13);           // left/top-ish
  lat.DrawLatexNDC(0.15, 0.95, "MicroBooNE Simulation");

  std::string pdf = out_prefix + "_neutron_efficiency_KE.pdf";
  std::string png = out_prefix + "_neutron_efficiency_KE.png";
  c.SaveAs(pdf.c_str());
  c.SaveAs(png.c_str());

  std::cout.setf(std::ios::fixed); std::cout<<std::setprecision(6);
  std::cout << "[Neutron Eff] weighted counts -- true=" << N_true_fsn
          << "  final=" << N_sel_prx  << "\n";
  if (N_true_fsn > 0.0) {
    std::cout << "[Neutron Eff] integrated " 
            << "  final=" << (N_sel_prx/N_true_fsn) << std::endl;
  }

}

static void MakeNeutronEfficiency_ProtGE50(
  const std::vector<RunGroup>& overlay_files,
  const SelCfg& cfg,
  const std::string& out_prefix,
  Long64_t max_events = -1)
{

  double N_true_fsn50 = 0.0, N_sel_prx50 = 0.0;

  const auto edges = NeutronKEBins(); const int nb = (int)edges.size()-1;
  TH1D h_den("h_den",";True Neutron KE (>50MeV Proton) (GeV);Efficiency", nb, edges.data());
  TH1D h_pre("h_pre",";True Neutron KE (>50MeV Proton) (GeV);Efficiency", nb, edges.data());
  TH1D h_sco("h_sco",";True Neutron KE (>50MeV Proton) (GeV);Efficiency", nb, edges.data());
  TH1D h_prx("h_prx",";True Neutron KE (>50MeV Proton) (GeV);Efficiency", nb, edges.data());
  h_den.Sumw2(); h_pre.Sumw2(); h_sco.Sumw2(); h_prx.Sumw2();

  for (const auto& rg : overlay_files) {
    std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(),"READ"));
    if (!f || f->IsZombie()) continue;
    TTree* tr = (TTree*)f->Get(kTreeName); if (!tr) continue;
    const double fscale = PerFileScale(Sample::kOverlay, rg.run_index, f.get());
    tr->SetBranchStatus("*",0);

    Bool_t numucc = kFALSE; float vx=0,vy=0,vz=0;
    tr->SetBranchStatus("CC1muNp0piNn_nu_mu_cc",1); tr->SetBranchAddress("CC1muNp0piNn_nu_mu_cc",&numucc);
    tr->SetBranchStatus("reco_nu_vtx_sce_x",1); tr->SetBranchAddress("reco_nu_vtx_sce_x",&vx);
    tr->SetBranchStatus("reco_nu_vtx_sce_y",1); tr->SetBranchAddress("reco_nu_vtx_sce_y",&vy);
    tr->SetBranchStatus("reco_nu_vtx_sce_z",1); tr->SetBranchAddress("reco_nu_vtx_sce_z",&vz);

    std::vector<float> *sx_v=nullptr,*sy_v=nullptr,*sz_v=nullptr,*ex_v=nullptr,*ey_v=nullptr,*ez_v=nullptr;
    std::vector<float> *len_v=nullptr,*pid_v=nullptr,*dist_v=nullptr,*trkscore_v=nullptr;
    std::vector<int>   *gen_v=nullptr,*bt_v=nullptr;

    tr->SetBranchStatus("nonprim_trk_sce_start_x_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_x_v",&sx_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_y_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_y_v",&sy_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_z_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_z_v",&sz_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_x_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_x_v",&ex_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_y_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_y_v",&ey_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_z_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_z_v",&ez_v);
    tr->SetBranchStatus("nonprim_trk_len_v",1);         tr->SetBranchAddress("nonprim_trk_len_v",&len_v);
    tr->SetBranchStatus("nonprim_trk_llr_pid_score_v",1); tr->SetBranchAddress("nonprim_trk_llr_pid_score_v",&pid_v);
    tr->SetBranchStatus("nonprim_trk_distance_v",1);    tr->SetBranchAddress("nonprim_trk_distance_v",&dist_v);
    tr->SetBranchStatus("nonprim_trk_score_v",1);       tr->SetBranchAddress("nonprim_trk_score_v",&trkscore_v);
    tr->SetBranchStatus("nonprim_pfp_generation",1);    tr->SetBranchAddress("nonprim_pfp_generation",&gen_v);
    tr->SetBranchStatus("nonprim_backtracked_tid",1);   tr->SetBranchAddress("nonprim_backtracked_tid",&bt_v);

    std::vector<int> *mc_trkid=nullptr,*mc_pdg=nullptr,*mc_mother=nullptr;
    std::vector<float> *mc_E=nullptr;
    std::vector<std::string>* mc_process=nullptr;
    tr->SetBranchStatus("all_mc_trkid",1);   tr->SetBranchAddress("all_mc_trkid",&mc_trkid);
    tr->SetBranchStatus("all_mc_pdg",1);     tr->SetBranchAddress("all_mc_pdg",&mc_pdg);
    tr->SetBranchStatus("all_mc_mother",1);  tr->SetBranchAddress("all_mc_mother",&mc_mother);
    tr->SetBranchStatus("all_mc_E",1);       tr->SetBranchAddress("all_mc_E",&mc_E);
    tr->SetBranchStatus("all_mc_process",1); tr->SetBranchAddress("all_mc_process",&mc_process);

    std::vector<double> *w1=nullptr,*w2=nullptr,*w3=nullptr;
    tr->SetBranchStatus("weight_splines_general_Spline",1);
    tr->SetBranchStatus("weight_TunedCentralValue_UBGenie",1);
    tr->SetBranchStatus("weight_neutron_argon_xsec",1);
    tr->SetBranchAddress("weight_splines_general_Spline",&w1);
    tr->SetBranchAddress("weight_TunedCentralValue_UBGenie",&w2);
    tr->SetBranchAddress("weight_neutron_argon_xsec",&w3);

    const Long64_t nent = tr->GetEntries();
    const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

    std::unordered_map<int,FSNInfo> fsn_map;
    std::vector<std::pair<int,double>> sec_protons;
    for (Long64_t ie=0; ie<lim; ++ie){
      tr->GetEntry(ie);
      if (!numucc) continue;
      if (!mc_trkid||!mc_pdg||!mc_mother||!mc_E) continue;
      double w = 1.0;
      if (!w1||!w2||!w3) continue;
      w *= (*w1)[0]*(*w2)[0]*(*w3)[0];
      if (!std::isfinite(w) || std::abs(w)>1e6) continue;
      w *= fscale;

      CollectTruth_FSN_and_SecProtons(*mc_trkid,*mc_pdg,*mc_mother,*mc_E,*mc_process, fsn_map, sec_protons);

      for (auto& kv : fsn_map) if (kv.second.has_p50) { h_den.Fill(kv.second.ke, w); N_true_fsn50 += w; }

      if (!sx_v||!ex_v||!len_v||!pid_v||!dist_v||!trkscore_v||!gen_v||!bt_v) continue;
      const int np = (int)len_v->size(); if (np==0) continue;

      std::unordered_set<int> seen_pre, seen_sco, seen_prx;
      bool pass_pre = false;
      for (int ip=0; ip<np; ++ip){
        if (gen_v->at(ip)!=2) continue;
        const double sx=sx_v->at(ip), sy=sy_v->at(ip), sz=sz_v->at(ip);
        const double ex=ex_v->at(ip), ey=ey_v->at(ip), ez=ez_v->at(ip);
        if (!FullyContained(sx,sy,sz,ex,ey,ez)) continue;

        pass_pre = true;

        const double L = len_v->at(ip);
        if (L < cfg.min_track_len) continue;

        const double pid = pid_v->at(ip);
        const double tsc = trkscore_v->at(ip);
        const double dist = dist_v->at(ip);
        const bool pass_vtxsep = (dist >= cfg.min_distance && dist <= cfg.max_distance);
        const bool pass_scores = (tsc >= cfg.min_trk_score && pid <= cfg.max_pid_score);
        const double prox = ComputePFPProximity(sx_v,sy_v,sz_v, ex_v,ey_v,ez_v, sx,sy,sz,ex,ey,ez);
        const float  dir  = ComputeRelativeDirection(vx,vy,vz, sx,sy,sz, ex,ey,ez);
        const bool pass_prxdir = (prox >= cfg.min_prox && dir >= cfg.min_direction);

        const int bt = (ip<(int)bt_v->size())? bt_v->at(ip) : 0;
        ProtonMatch m = MatchSignalProton_WithKE(bt, *mc_trkid,*mc_pdg,*mc_mother,*mc_E,*mc_process);
        if (m.fsn_tid==0) continue;
        if (m.p_ke < 0.05) continue;                        // require >¥50 MeV proton
        auto it = fsn_map.find(m.fsn_tid); if (it==fsn_map.end() || !it->second.has_p50) continue;

        if (pass_pre && !seen_pre.count(m.fsn_tid))   { h_pre.Fill(it->second.ke, w); seen_pre.insert(m.fsn_tid); }
        if (pass_vtxsep && pass_scores && !seen_sco.count(m.fsn_tid)) { h_sco.Fill(it->second.ke, w); seen_sco.insert(m.fsn_tid); }
        if (pass_vtxsep && pass_scores && pass_prxdir && !seen_prx.count(m.fsn_tid)) { h_prx.Fill(it->second.ke, w); seen_prx.insert(m.fsn_tid); N_sel_prx50 += w; }
      }
    }
  }

  TH1D e_pre(h_pre); e_pre.Divide(&h_den);
  TH1D e_sco(h_sco); e_sco.Divide(&h_den);
  TH1D e_prx(h_prx); e_prx.Divide(&h_den);

  e_prx.SetName("eff_proxdir"); e_sco.SetName("eff_scorevtx"); e_pre.SetName("eff_pre");
  e_prx.SetMinimum(0.0); e_prx.SetMaximum(0.55);
  e_prx.SetMarkerStyle(20); e_sco.SetMarkerStyle(20); e_pre.SetMarkerStyle(20);
  e_prx.SetMarkerColor(4); e_sco.SetMarkerColor(2); e_pre.SetMarkerColor(1);
  e_prx.SetLineWidth(2); e_sco.SetLineWidth(2); e_pre.SetLineWidth(2);
  e_prx.SetLineColor(4); e_sco.SetLineColor(2); e_pre.SetLineColor(1);

  TCanvas c("c_eff50","",900,650);
  e_prx.Draw("E1"); e_sco.Draw("E1 SAME"); e_pre.Draw("E1 SAME");
  TLegend leg(0.15,0.7,0.53,0.9); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42);
  leg.AddEntry(&e_pre, "Preselection", "lep");
  leg.AddEntry(&e_sco, "Scores and Vertex Separation", "lep");
  leg.AddEntry(&e_prx, "Proximity and Direction", "lep");
  leg.Draw();
  std::string base = out_prefix + "_neutron_efficiency_KE_50MeV_p_coarse_bin";

  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.04);          // <-- explicit size!
  lat.SetTextAlign(13);           // left/top-ish
  lat.DrawLatexNDC(0.15, 0.95, "MicroBooNE Simulation");

  c.SaveAs((base + ".pdf").c_str());
  c.SaveAs((base + ".png").c_str());

  std::cout<<std::fixed<<std::setprecision(6)
         << "[Neutron Eff, p>=50 MeV] true=" << N_true_fsn50
         << "  final=" << N_sel_prx50 << "\n";
  if (N_true_fsn50>0.0) {
    std::cout << "[Neutron Eff, p>=50 MeV final=" << (N_sel_prx50/N_true_fsn50) << std::endl;
  }
}

static std::vector<double> ProtonKEBins() {
  return {0.00,0.025,0.05,0.075,0.10,0.125,0.15,0.175,0.20,0.225,0.25};
}

static void MakeSecondaryProtonEfficiency(
  const std::vector<RunGroup>& overlay_files,
  const SelCfg& cfg,
  const std::string& out_prefix,
  Long64_t max_events = -1)
{
  const auto edges = ProtonKEBins(); const int nb = (int)edges.size()-1;
  TH1D h_den("h_den_p",";True Secondary Proton Kinetic Energy (GeV);Efficiency", nb, edges.data());
  TH1D h_pre("h_pre_p",";True Secondary Proton Kinetic Energy (GeV);Efficiency", nb, edges.data());
  TH1D h_sco("h_sco_p",";True Secondary Proton Kinetic Energy (GeV);Efficiency", nb, edges.data());
  TH1D h_prx("h_prx_p",";True Secondary Proton Kinetic Energy (GeV);Efficiency", nb, edges.data());
  h_den.Sumw2(); h_pre.Sumw2(); h_sco.Sumw2(); h_prx.Sumw2();

  for (const auto& rg : overlay_files) {
    std::unique_ptr<TFile> f(TFile::Open(rg.file.c_str(),"READ"));
    if (!f || f->IsZombie()) continue;
    TTree* tr = (TTree*)f->Get(kTreeName); if (!tr) continue;
    const double fscale = PerFileScale(Sample::kOverlay, rg.run_index, f.get());
    tr->SetBranchStatus("*",0);

    Bool_t numucc = kFALSE; float vx=0,vy=0,vz=0;
    tr->SetBranchStatus("CC1muNp0piNn_nu_mu_cc",1); tr->SetBranchAddress("CC1muNp0piNn_nu_mu_cc",&numucc);
    tr->SetBranchStatus("reco_nu_vtx_sce_x",1); tr->SetBranchAddress("reco_nu_vtx_sce_x",&vx);
    tr->SetBranchStatus("reco_nu_vtx_sce_y",1); tr->SetBranchAddress("reco_nu_vtx_sce_y",&vy);
    tr->SetBranchStatus("reco_nu_vtx_sce_z",1); tr->SetBranchAddress("reco_nu_vtx_sce_z",&vz);

    std::vector<float> *sx_v=nullptr,*sy_v=nullptr,*sz_v=nullptr,*ex_v=nullptr,*ey_v=nullptr,*ez_v=nullptr;
    std::vector<float> *len_v=nullptr,*pid_v=nullptr,*dist_v=nullptr,*trkscore_v=nullptr;
    std::vector<int>   *gen_v=nullptr,*bt_v=nullptr;

    tr->SetBranchStatus("nonprim_trk_sce_start_x_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_x_v",&sx_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_y_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_y_v",&sy_v);
    tr->SetBranchStatus("nonprim_trk_sce_start_z_v",1); tr->SetBranchAddress("nonprim_trk_sce_start_z_v",&sz_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_x_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_x_v",&ex_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_y_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_y_v",&ey_v);
    tr->SetBranchStatus("nonprim_trk_sce_end_z_v",1);   tr->SetBranchAddress("nonprim_trk_sce_end_z_v",&ez_v);
    tr->SetBranchStatus("nonprim_trk_len_v",1);         tr->SetBranchAddress("nonprim_trk_len_v",&len_v);
    tr->SetBranchStatus("nonprim_trk_llr_pid_score_v",1); tr->SetBranchAddress("nonprim_trk_llr_pid_score_v",&pid_v);
    tr->SetBranchStatus("nonprim_trk_distance_v",1);    tr->SetBranchAddress("nonprim_trk_distance_v",&dist_v);
    tr->SetBranchStatus("nonprim_trk_score_v",1);       tr->SetBranchAddress("nonprim_trk_score_v",&trkscore_v);
    tr->SetBranchStatus("nonprim_pfp_generation",1);    tr->SetBranchAddress("nonprim_pfp_generation",&gen_v);
    tr->SetBranchStatus("nonprim_backtracked_tid",1);   tr->SetBranchAddress("nonprim_backtracked_tid",&bt_v);

    std::vector<int> *mc_trkid=nullptr,*mc_pdg=nullptr,*mc_mother=nullptr;
    std::vector<float> *mc_E=nullptr;
    std::vector<std::string>* mc_process=nullptr;
    tr->SetBranchStatus("all_mc_trkid",1);   tr->SetBranchAddress("all_mc_trkid",&mc_trkid);
    tr->SetBranchStatus("all_mc_pdg",1);     tr->SetBranchAddress("all_mc_pdg",&mc_pdg);
    tr->SetBranchStatus("all_mc_mother",1);  tr->SetBranchAddress("all_mc_mother",&mc_mother);
    tr->SetBranchStatus("all_mc_E",1);       tr->SetBranchAddress("all_mc_E",&mc_E);
    tr->SetBranchStatus("all_mc_process",1); tr->SetBranchAddress("all_mc_process",&mc_process);

    std::vector<double> *w1=nullptr,*w2=nullptr,*w3=nullptr;
    tr->SetBranchStatus("weight_splines_general_Spline",1);
    tr->SetBranchStatus("weight_TunedCentralValue_UBGenie",1);
    tr->SetBranchStatus("weight_neutron_argon_xsec",1);
    tr->SetBranchAddress("weight_splines_general_Spline",&w1);
    tr->SetBranchAddress("weight_TunedCentralValue_UBGenie",&w2);
    tr->SetBranchAddress("weight_neutron_argon_xsec",&w3);

    const Long64_t nent = tr->GetEntries();
    const Long64_t lim  = (max_events>0 && max_events<nent)?max_events:nent;

    std::unordered_map<int,FSNInfo> fsn_dummy;
    std::vector<std::pair<int,double>> sec_protons;
    for (Long64_t ie=0; ie<lim; ++ie){
      tr->GetEntry(ie);
      if (!numucc) continue;
      if (!mc_trkid||!mc_pdg||!mc_mother||!mc_E) continue;
      double w = 1.0;
      if (!w1||!w2||!w3) continue;
      w *= (*w1)[0]*(*w2)[0]*(*w3)[0];
      if (!std::isfinite(w) || std::abs(w)>1e6) continue;
      w *= fscale;

      CollectTruth_FSN_and_SecProtons(*mc_trkid,*mc_pdg,*mc_mother,*mc_E,*mc_process, fsn_dummy, sec_protons);
      for (auto& pr : sec_protons) h_den.Fill(pr.second, w);

      if (!sx_v||!ex_v||!len_v||!pid_v||!dist_v||!trkscore_v||!gen_v||!bt_v) continue;
      const int np = (int)len_v->size(); if (np==0) continue;

      std::unordered_set<int> seen_pre, seen_sco, seen_prx;

      bool pass_pre = false;

      for (int ip=0; ip<np; ++ip){
        if (gen_v->at(ip)!=2) continue;
        const double sx=sx_v->at(ip), sy=sy_v->at(ip), sz=sz_v->at(ip);
        const double ex=ex_v->at(ip), ey=ey_v->at(ip), ez=ez_v->at(ip);
        if (!FullyContained(sx,sy,sz,ex,ey,ez)) continue;

	pass_pre = true;

        const double L = len_v->at(ip);
        if (L < cfg.min_track_len) continue;

        const double pid = pid_v->at(ip);
        const double tsc = trkscore_v->at(ip);
        const double dist = dist_v->at(ip);
        const bool pass_vtxsep = (dist >= cfg.min_distance && dist <= cfg.max_distance);
        const bool pass_scores = (tsc >= cfg.min_trk_score && pid <= cfg.max_pid_score);
        const double prox = ComputePFPProximity(sx_v,sy_v,sz_v, ex_v,ey_v,ez_v, sx,sy,sz,ex,ey,ez);
        const float  dir  = ComputeRelativeDirection(vx,vy,vz, sx,sy,sz, ex,ey,ez);
        const bool pass_prxdir = (prox >= cfg.min_prox && dir >= cfg.min_direction);

        const int bt = (ip<(int)bt_v->size())? bt_v->at(ip) : 0;
        ProtonMatch m = MatchSignalProton_WithKE(bt, *mc_trkid,*mc_pdg,*mc_mother,*mc_E,*mc_process);
        if (m.fsn_tid==0) continue;               // must be sec proton from FS neutron

        if (pass_pre && !seen_pre.count(bt))   { h_pre.Fill(m.p_ke, w); seen_pre.insert(bt); }
        if (pass_vtxsep && pass_scores && !seen_sco.count(bt)) { h_sco.Fill(m.p_ke, w); seen_sco.insert(bt); }
        if (pass_vtxsep && pass_scores && pass_prxdir && !seen_prx.count(bt)) { h_prx.Fill(m.p_ke, w); seen_prx.insert(bt); }
      }
    }
  }

  TH1D e_pre(h_pre); e_pre.Divide(&h_den);
  TH1D e_sco(h_sco); e_sco.Divide(&h_den);
  TH1D e_prx(h_prx); e_prx.Divide(&h_den);

  e_prx.SetName("effP_proxdir"); e_sco.SetName("effP_scorevtx"); e_pre.SetName("effP_pre");
  e_prx.SetMinimum(0.0); e_prx.SetMaximum(0.6);
  e_prx.SetMarkerStyle(20); e_sco.SetMarkerStyle(20); e_pre.SetMarkerStyle(20);
  e_prx.SetMarkerColor(4); e_sco.SetMarkerColor(2); e_pre.SetMarkerColor(1);
  e_prx.SetLineWidth(2); e_sco.SetLineWidth(2); e_pre.SetLineWidth(2);
  e_prx.SetLineColor(4); e_sco.SetLineColor(2); e_pre.SetLineColor(1);

  TCanvas c("c_effP","",900,650);
  e_prx.Draw("E1"); e_sco.Draw("E1 SAME"); e_pre.Draw("E1 SAME");
  TLegend leg(0.15,0.7,0.53,0.9); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextFont(42);
  leg.AddEntry(&e_pre, "Preselection", "lep");
  leg.AddEntry(&e_sco, "Scores and Vertex Separation", "lep");
  leg.AddEntry(&e_prx, "Proximity and Direction", "lep");
  leg.Draw();
  std::string base = out_prefix + "_secondary_proton_efficiency_KE_coarse_bin";

  TLatex lat;
  lat.SetTextFont(42);
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.04);          // <-- explicit size!
  lat.SetTextAlign(13);           // left/top-ish
  lat.DrawLatexNDC(0.15, 0.95, "MicroBooNE Simulation");

  c.SaveAs((base + ".pdf").c_str());
  c.SaveAs((base + ".png").c_str());
}

static void MakeAllEfficiencyAddons(const std::vector<RunGroup>& overlay_files,
                                    const SelCfg& cfg,
                                    const std::string& out_prefix,
                                    Long64_t max_events = -1)
{
  MakeNeutronEfficiency_ProtGE50(overlay_files, cfg, out_prefix, max_events);
  MakeSecondaryProtonEfficiency(overlay_files, cfg, out_prefix, max_events);
}

int main() {
  std::vector<RunGroup> overlay_files = {
  /*  {"/pnfs/uboone/persistent/users/mhernan/processed/run1_CC1muNp0piNn_overlay_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/mhernan/processed/run2_CC1muNp0piNn_overlay_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_overlay_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_overlay_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_overlay_POST_multisim.root", 4},
     };
*/  {"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root", 0},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_overlay_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_overlay_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_overlay_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_overlay_POST_multisim.root", 4},
   };
/*  std::vector<RunGroup> dirt_files = {
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_dirt_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_dirt_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_dirt_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_dirt_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_dirt_POST_multisim.root", 4},
  };
  std::vector<RunGroup> ext_files = {
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_EXT_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2_CC1muNp0piNn_EXT_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_EXT_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_EXT_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_EXT_POST_multisim.root", 4},
  };
*/  
/*  std::vector<RunGroup> nuwro_files = {
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run1_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run2b_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run3_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run4_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/run5_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 4},
  };
*/
/*  std::vector<RunGroup> overlay_files = {
    {"/pnfs/uboone/persistent/users/mhernan/processed/run1_CC1muNp0piNn_overlay_POST_multisim.root", 0},
    {"/pnfs/uboone/persistent/users/mhernan/processed/run2_CC1muNp0piNn_overlay_POST_multisim.root", 1},
    {"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_overlay_POST_multisim.root", 2},
    {"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_overlay_POST_multisim.root", 3},
    {"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_overlay_POST_multisim.root", 4},
  };
*/
  std::vector<RunGroup> dirt_files = {
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run1_CC1muNp0piNn_dirt_POST_multisim.root", 0},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run2_CC1muNp0piNn_dirt_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_dirt_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_dirt_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_dirt_POST_multisim.root", 4},
  };

  std::vector<RunGroup> ext_files = {
    {"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_beamoffEXT_surprise_blipPlus_20260522.root", 0},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run2_CC1muNp0piNn_EXT_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_EXT_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_EXT_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_EXT_POST_multisim.root", 4},
  };

  std::vector<RunGroup> nuwro_files = {
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run1_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 0},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run2b_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run3_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run4_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/mhernan/processed/run5_CC1muNp0piNn_nuwro_fake_data_POST_multisim.root", 4},
  };
  std::vector<RunGroup> data_files = {
    {"/Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_11_BNB_opendata_surprise_blipPlus_20260522.root", 0},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run2_Beam_On_POST_multisim.root", 1},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run3_Beam_On_POST_multisim.root", 2},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run4_Beam_On_POST_multisim.root", 3},
    //{"/pnfs/uboone/persistent/users/birwin/PostProcessedNTuples/temp_PostProcess/Run5_Beam_On_POST_multisim.root", 4},
  };

  SelCfg cfg;
  cfg.min_track_len = 0.0;
  cfg.max_pid_score = -0.15;
  cfg.require_same_slice = true;
  cfg.min_distance = 10.0;
  cfg.max_distance = 200.0;
  cfg.min_trk_score = 0.5;
  cfg.min_prox = 21.0;
  cfg.min_direction = 0.2;


  MakeSecondaryProtonPlots(overlay_files, dirt_files, ext_files, data_files, cfg, "Run1-5_BNB_Data_candidates", false, -1);

  //MakeNeutronEfficiencyPlots(overlay_files, cfg, "Run1-5", -1);
  //MakeAllEfficiencyAddons(overlay_files, cfg, "Run1-5", -1);

  return 0;
}

