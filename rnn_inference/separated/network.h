#pragma once
#include <torch/torch.h>
#include <torch/script.h>
#include <string>
#include <vector>

class RNNPredictor {
private:
    torch::jit::script::Module model;
    torch::Device device;
    bool model_loaded;
 
public:
    RNNPredictor();
    bool loadModel(const std::string& model_path);
    std::vector<float> predict(const std::vector<std::vector<float>>& sequence, 
                                int actual_length);
};
