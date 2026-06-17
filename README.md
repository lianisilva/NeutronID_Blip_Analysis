WELCOME!

This is the ENTIRETY of my "neutron identification with proton blips" analysis code excluding the prepared/output .root files, .txt files and other types of files you can find in the .gitignore.
Most of the ROOT dataset files that are referenced within the code can be found at /exp/uboone/data/users/wforeman/uboonecode-mcc10/anatrees/, with access to the MicroBooNE GPVMs.

The main/most recently updated folder is the "pelee" folder which contains all the code needed to make our final primary (n,p) selection.
It also contains "xsec_analyzer" which is mostly code from others with tweaks/additions to make it work for our analysis.
(The code from others is that from Steven's framework as well as from Burke and Miguel's work on proton tracks for neutron ID.)

The "analyze_neutrons" folder is from before we switched to the PeLEE NTuples and instead use a specifically curated run 3 overlay root file.
This folder also includes our blip validation work ("validation") and the "hit_info.C" macro which was used to create the .txt file input for the RNN before we used Libtorch/had RNN variables directly in our NTuples.

Speaking of, we have an "rnn_inference" folder here which contains all of the scripts needed to run the RNN model within ROOT/C++ utilizing the LibTorch library.
More information on how to do this is in this guide document -> https://docs.google.com/document/d/1PF5z-vx-rHSHfrryYS5lpSa-mMfgve9uZSoTvvZcdKU/edit?pli=1&tab=t.0.
This document also explains how to input blips into the BDT which gives you a score based on how signal or background-like the blip is.
