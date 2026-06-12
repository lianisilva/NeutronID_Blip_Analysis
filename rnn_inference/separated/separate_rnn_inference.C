// ROOT macro for RNN inference
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "network.h"
#include <iostream>
#include <vector>
#include <algorithm>

void rnn_inference(const char* input_file = "input.root", 
                   const char* output_file = "output.root", 
                   const char* model_path = "model_traced.pt") {
    
    std::cout << "=== ROOT Macro RNN Inference ===" << std::endl;
    
    // Initialize predictor
    RNNPredictor predictor;
    if (!predictor.loadModel(model_path)) {
        std::cerr << "Cannot load model, exiting." << std::endl;
        return;
    }
    
    // Open input file
    TFile* inFile = TFile::Open(input_file, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Cannot open input file: " << input_file << std::endl;
        return;
    }
    
    TTree* inTree = (TTree*)inFile->Get("blipana/anatree");
    if (!inTree) {
        std::cerr << "Cannot find tree" << std::endl;
        inFile->Close();
        return;
    }
    
    // Setup branches
    std::vector<float>* feature1 = nullptr;
    std::vector<float>* feature2 = nullptr;
    std::vector<float>* feature3 = nullptr;
    Int_t cluster_size = 0;
    
    inTree->SetBranchAddress("feature1", &feature1);
    inTree->SetBranchAddress("feature2", &feature2);
    inTree->SetBranchAddress("feature3", &feature3);
    inTree->SetBranchAddress("cluster_size", &cluster_size);
    
    // Create output
    TFile* outFile = TFile::Open(output_file, "RECREATE");
    TTree* outTree = new TTree("rnn_results", "RNN Predictions");
    
    std::vector<float> rnn_output;
    Float_t prediction_score = 0;
    Int_t sequence_length = 0;
    
    outTree->Branch("rnn_output", &rnn_output);
    outTree->Branch("prediction_score", &prediction_score);
    outTree->Branch("sequence_length", &sequence_length);
    
    // Process events
    Long64_t nEntries = inTree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    
    for (Long64_t i = 0; i < nEntries; ++i) {
        if (i % 100 == 0) {
            std::cout << "Processing event " << i << "/" << nEntries << std::endl;
        }
        
        inTree->GetEntry(i);
        
        rnn_output.clear();
        prediction_score = -1.0;
        sequence_length = 0;
        
        if (!feature1 || !feature2 || !feature3 || 
            feature1->empty() || feature2->empty() || feature3->empty()) {
            outTree->Fill();
            continue;
        }
        
        int actual_seq_length = std::min({
            static_cast<int>(feature1->size()), 
            static_cast<int>(feature2->size()), 
            static_cast<int>(feature3->size()),
            cluster_size
        });
        
        if (actual_seq_length <= 0) {
            outTree->Fill();
            continue;
        }
        
        sequence_length = actual_seq_length;
        
        std::vector<std::vector<float>> sequence;
        for (int t = 0; t < actual_seq_length; ++t) {
            sequence.push_back({
                (*feature1)[t],
                (*feature2)[t], 
                (*feature3)[t]
            });
        }
        
        auto prediction = predictor.predict(sequence, actual_seq_length);
        
        if (!prediction.empty()) {
            rnn_output = prediction;
            prediction_score = prediction[0];
        }
        
        outTree->Fill();
    }
    
    outTree->Write();
    outFile->Close();
    inFile->Close();
    
    std::cout << "✓ Processing complete!" << std::endl;
    std::cout << "Results saved to: " << output_file << std::endl;
}
