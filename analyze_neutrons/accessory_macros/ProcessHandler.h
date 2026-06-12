#ifndef PROCESS_HANDLER_H
#define PROCESS_HANDLER_H

#include <unordered_map>
#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include "TH1D.h"

class ProcessHandler {
public:
    using HandlerFunc = std::function<void(double cos, int p_mom)>;

    // Structure to hold energy statistics for each process
    struct EnergyStats {
        double total_energy = 0.0; int count = 0;
        void addEnergy(double energy) { total_energy += energy; count++; }
        double getAverage() const { return count > 0 ? total_energy / count : 0.0; }
    };

    // Constructor takes references to your histograms and counters
    ProcessHandler(//TH1D* h_cos_nIne_gamma,
                   //TH1D* h_cos_nCap_gamma,
                   //TH1D* h_cos_eBrem_gamma,
                   //TH1D* h_cos_nOther_gamma,
                   int&  total_nCap_gammas)
        : //h_cos_nIne_gamma_(h_cos_nIne_gamma),
          //h_cos_nCap_gamma_(h_cos_nCap_gamma),
          //h_cos_eBrem_gamma_(h_cos_eBrem_gamma),
          //h_cos_nOther_gamma_(h_cos_nOther_gamma),
          total_nCap_gammas_(total_nCap_gammas) {
        setupDefaultHandlers();
    }

    // Main processing function
    void handleProcess(const std::string& process_name, double cos, int p_mom, float energy) {
        // Track energy for this process
        energy_stats_[process_name].addEnergy(energy);
        
        auto it = process_handlers_.find(process_name);
        if (it != process_handlers_.end()) {
            it->second(cos, p_mom);
        } else {
            h_cos_nOther_gamma_->Fill(cos);
            std::cout << "Unknown process: " << process_name << std::endl;
        }
    }

    double getAverageEnergy(const std::string& process_name) const {
        auto it = energy_stats_.find(process_name);
        return it != energy_stats_.end() ? it->second.getAverage() : 0.0;
    }

    std::vector<std::string> getProcessNames() const {
        std::vector<std::string> names;
        for (const auto& pair : energy_stats_) {
            names.push_back(pair.first);
        }
        return names;
    }

    // Print energy statistics for all processes
    void printEnergyStatistics() const {
        std::cout<<"\n=== Gamma Energy Stat by Process ===";
        std::cout<<"\n---------------------------------------\n";

        for (const auto& pair : energy_stats_) {
          const std::string& name = pair.first;
          const EnergyStats& stats = pair.second;
	  printf("Name: %-25s Count: %-10i  Energy (Avg) =%6.2f MeV\n", name.c_str(), stats.count, stats.getAverage() );
        }
    }

    // Get energy statistics map (for further analysis)
    const std::unordered_map<std::string, EnergyStats>& getEnergyStats() const {
        return energy_stats_;
    }

    // Add new process handlers dynamically (if needed)
    void addHandler(const std::string& process_name, HandlerFunc handler) {
        process_handlers_[process_name] = handler;
    }

// feels like this should be on top
private:
    std::unordered_map<std::string, HandlerFunc> process_handlers_;
    std::unordered_map<std::string, EnergyStats> energy_stats_;

    // References to your analysis objects
    //TH1D* h_cos_nIne_gamma_;
    //TH1D* h_cos_nCap_gamma_;
    //TH1D* h_cos_eBrem_gamma_;
    //TH1D* h_cos_nOther_gamma_;
    int&  total_nCap_gammas_;

    void setupDefaultHandlers() {
      process_handlers_ = {
        {"neutronInelastic", [this](double cos, int p_mom) //{ h_cos_nIne_gamma_->Fill(cos); }},
        {"nCapture", [this](double cos, int p_mom) { h_cos_nCap_gamma_->Fill(cos);
                                                	     total_nCap_gammas_++; }},
        {"eBrem", [this](double cos, int p_mom) { h_cos_eBrem_gamma_->Fill(cos); }},
        {"annihil", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"hBertiniCaptureAtRest", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"pi-Inelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"pi+Inelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"protonInelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"dInelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"sigma+Inelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"kaon0LInelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"kaon+Inelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"kaon-Inelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"tInelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"lambdaInelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"alphaInelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"He3Inelastic", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"Decay", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"photonNuclear", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"electronNuclear", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"positronNuclear", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"muMinusCaptureAtRest", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }},
        {"primary", [this](double cos, int p_mom) { h_cos_nOther_gamma_->Fill(cos); }}
      };
    }
};

#endif
