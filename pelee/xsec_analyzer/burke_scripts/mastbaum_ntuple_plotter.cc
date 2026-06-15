#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

bool passes_signal_definition(int nfsps, int pdg[], float px[], float py[], float pz[], float energy[], bool flagCC0pi) {
  if (!flagCC0pi) return false;

  double muon_p = -1.0;
  std::vector<double> proton_ps;
  bool has_high_ke_neutron = false;

  for (int i = 0; i < nfsps; ++i) {
    int id = pdg[i];
    double p = std::sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
    double mass = 0.0;
    if (id == 13) {
      muon_p = p;
    } else if (id == 2212) {
      proton_ps.push_back(p);
    } else if (id == 2112) {
      mass = 0.9395654; // neutron mass
      double ke = energy[i] - mass;
      if (ke > 0.1) has_high_ke_neutron = true;
    }
  }

  if (muon_p < 0.1 || muon_p > 1.2) return false;
  if (proton_ps.empty()) return false;

  double lead_proton_p = *std::max_element(proton_ps.begin(), proton_ps.end());
  if (lead_proton_p < 0.25 || lead_proton_p > 1.0) return false;
  if (!has_high_ke_neutron) return false;

  return true;
}

int main() {
  /*std::vector<std::string> generators = {
    "genie_v2_12_10",
    "genie_v2_12_10.mec",
    "genie_v3_00_06",
    "neut_5_4_0_1",
    "nuwro_19_02_1"
  };*/

  std::vector<std::string> generators = {
    "neut_5_4_0_1",
    "AR23_20i_00_000",
    "gibuu_2025",
    "NuWro_numu_21_09_2",
    "genie_v3_00_06",
  };
  //std::string folder = "/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/";
  std::vector<std::string> paths_prefix = {
    "/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.",
    "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/",
    "/exp/uboone/data/users/cnguyen/CC0Pi_Selection/modelcomparson/",
    "/pnfs/uboone/persistent/users/apapadop/NuWro_Samples/",
    "/pnfs/uboone/persistent/users/mastbaum/tuning2022/mc/bnb_ub/flat/bnb.ub.num.",
  };

  std::vector<std::string> paths_post = {
    ".flat.root",
    "/14_1000180400_CC_v3_4_0_AR23_20i_00_000.flat.root",
    "/GiBUU_2025.flat_250files.root",
    "_2000000.flat.root",
    ".flat.root",
  };

  TFile* fout = new TFile("ClosureTest_signal_counts.root", "RECREATE");

  for (size_t i = 0; i < generators.size(); ++i) {
    std::string gen = generators[i];
    std::string path_pre = paths_prefix[i];
    std::string path_post = paths_post[i];
    //std::string path = folder + "bnb.ub.num." + gen + ".flat.root";
    std::string path = path_pre+gen+path_post;

    TFile* file = TFile::Open(path.c_str(), "READ");
    if (!file || file->IsZombie()) {
      std::cerr << "Failed to open file: " << path << std::endl;
      continue;
    }

    TTree* tree = dynamic_cast<TTree*>(file->Get("FlatTree_VARS"));
    if (!tree) {
      std::cerr << "Failed to get tree in file: " << path << std::endl;
      file->Close();
      continue;
    }

    int nfsps, pdg[100];
    float px[100], py[100], pz[100], energy[100];
    bool flagCC0pi;
    double scaleFactor = 0.0;
    float weight = 0.0;

    tree->SetBranchAddress("nfsp", &nfsps);
    tree->SetBranchAddress("pdg", pdg);
    tree->SetBranchAddress("px", px);
    tree->SetBranchAddress("py", py);
    tree->SetBranchAddress("pz", pz);
    tree->SetBranchAddress("E", energy);
    tree->SetBranchAddress("flagCC0pi", &flagCC0pi);
    tree->SetBranchAddress("fScaleFactor", &scaleFactor);
    tree->SetBranchAddress("Weight",&weight);

    Long64_t nentries = tree->GetEntries();
    int count = 0;
    for (Long64_t j = 0; j < nentries; ++j) {
      tree->GetEntry(j);
      if (passes_signal_definition(nfsps, pdg, px, py, pz, energy, flagCC0pi)) 
      {
	if (gen == "gibuu_2025" ) { count += weight; }
	else { ++count; }
      }
    }
    std::cout << "gen " << gen << " and count: " << count << std::endl;
    std::cout << "scale factor: " << scaleFactor << std::endl;

    double xsec = count * scaleFactor * 40; //times 40 to make it per nucleus

    std::string hist_name = "h_" + gen;
    TH1D* h = new TH1D(hist_name.c_str(), ("Signal count for " + gen).c_str(), 1, 0, 1);
    h->SetDirectory(nullptr);
    h->SetBinContent(1, xsec);
    fout->cd();
    h->Write();
    file->Close();
  }

  fout->Close();

  return 0;
}

