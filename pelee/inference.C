#include <xgboost/c_api.h>
#include <iostream>
#include <vector>
#include <cmath>

void inference() {

    BoosterHandle booster;
    XGBoosterCreate(nullptr, 0, &booster);
    XGBoosterLoadModel(booster, "best_models/n_p_select_noShow_100k.model");

    // 1 event x 9 features, row-major
    std::vector<float> data = {1.5f, -0.3f, 2.1f, 1.5f, -0.3f, 2.1f, 1.5f, -0.3f, 2.1f};

    DMatrixHandle dmat;
    XGDMatrixCreateFromMat(data.data(), 1, 9, NAN, &dmat);

    bst_ulong outLen;
    const float* outResult;
    // output_margin=1 gives raw logitraw score, matching predict(output_margin=True)
    XGBoosterPredict(booster, dmat, /*output_margin=*/1, 0, 0, &outLen, &outResult);

    std::cout << "Raw margin: " << outResult[0] << std::endl;

    XGDMatrixFree(dmat);
    XGBoosterFree(booster);
}

