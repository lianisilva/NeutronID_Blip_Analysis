#!/bin/bash
export TORCH_ROOT=/Users/liani/Library/libtorch
export ROOT_INCLUDE_PATH="/Users/liani/Library/libtorch/include:/Users/liani/Library/libtorch/include/torch/csrc/api/include:$ROOT_INCLUDE_PATH"
export DYLD_LIBRARY_PATH=$TORCH_ROOT/lib:$DYLD_LIBRARY_PATH
