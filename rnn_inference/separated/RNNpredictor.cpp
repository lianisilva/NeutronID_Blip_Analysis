#include "network.h"
#include <iostream>

RNNPredictor::RNNPredictor() : device(torch::kCPU), model_loaded(false) {}

bool RNNPredictor::loadModel(const std::string& model_path) {
    try {
        model = torch::jit::load(model_path);
        model.to(device);
        model.eval();
        model_loaded = true;
        std::cout << "Model loaded: " << model_path << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to load model: " << e.what() << std::endl;
        return false;
    }
}

std::vector<float> RNNPredictor::predict(
    const std::vector<std::vector<float>>& sequence, 
    int actual_length) {
    
    if (!model_loaded) {
        std::cerr << "Model not loaded!" << std::endl;
        return {};
    }
    
    try {
        int seq_len = sequence.size();
        int n_features = sequence[0].size();
        
        auto tensor = torch::zeros({1, seq_len, n_features});
        for (int i = 0; i < seq_len; ++i) {
            for (int j = 0; j < n_features; ++j) {
                tensor[0][i][j] = sequence[i][j];
            }
        }
        
        auto lengths_tensor = torch::tensor({actual_length}, torch::kLong);
        
        torch::NoGradGuard no_grad;
        std::vector<torch::jit::IValue> inputs;
        inputs.push_back(tensor);
        inputs.push_back(lengths_tensor);
        
        at::Tensor output = model.forward(inputs).toTensor();
        
        std::vector<float> result;
        auto accessor = output.accessor<float, 2>();
        for (int i = 0; i < output.size(1); ++i) {
            result.push_back(accessor[0][i]);
        }
        return result;
    } catch (const std::exception& e) {
        std::cerr << "Prediction failed: " << e.what() << std::endl;
        return {};
    }
}
