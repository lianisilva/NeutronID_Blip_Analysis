#!/bin/bash

echo "Compiling blip_utils shared library..."

g++ -shared -fPIC -o libbliputils.so ../rnn_inference/blip_utils.cpp \
    -std=c++17 \
    $(root-config --cflags --libs) \
    -I$HOME/Library/libtorch/include \
    -I$HOME/Library/libtorch/include/torch/csrc/api/include \
    -L$HOME/Library/libtorch/lib \
    -Wl,-rpath,$HOME/Library/libtorch/lib \
    -ltorch -ltorch_cpu -lc10

if [ $? -eq 0 ]; then
    echo "Library compiled successfully: libbliputils.so"
else
    echo "Compilation failed!"
    exit 1
fi
