# WELCOME!

## This is the ENTIRETY of my "neutron identification with proton blips" analysis code excluding the prepared/output .root files, .txt files and other types of files you can find in the .gitignore.

  * Most of the ROOT dataset files that are referenced within the code can be found (with access to the MicroBooNE GPVMs) at:

    * `/exp/uboone/data/users/wforeman/uboonecode-mcc10/anatrees/`

## * The `"pelee"` folder:
  * The main/most recently updated folder which contains all the code needed to make our final primary (n,p) selection.

  * It also contains `"xsec_analyzer"` which is mostly code from others with tweaks/additions to make it work for our analysis.
    * The code from "others" is that from Steven's framework as well as from Burke and Miguel's work on proton tracks for neutron ID.

## * The `"analyze_neutrons"` folder:
  * OLD; This is from before we switched to the PeLEE NTuples and instead use a specifically curated run 3 overlay root file.

  * The file `"neutrino_analysis.C"` here does practically the same thing as the more updated/recent main analysis macro `"pelee/surprise_nu_analysis.C"` but there are some useful calculations and functions that give some of the statistics we cite and that haven't been transferred to the new macro.

  * This folder also includes:
    * Our blip validation work (`"validation"`)
    * The `"hit_info.C"` macro which was used to create the .txt file input for the RNN before we used Libtorch; before we had RNN variables directly in our NTuples.

## * For more information on how to optimally use the standalone Python RNN see this document:
  * https://docs.google.com/document/d/1GcmniWuQfjSWTkSgCd5PSh_i_p5iKdFHEq_olMJ1syc/edit?usp=sharing

  * Though we no longer use this workflow procedure since we added RNN variables to the NTuples, this old procedure requires:
    * a local ROOT file (passed to)-> our (Python notebook) model
    * the model created ROOT file (passed back to)-> the local ROOT macro to use in the analysis.
  * Therefore, the only remaining use case for this workflow is dedicated RNN improvement testing, i.e. testing several different model tweaks before choosing the best model to integrate into the selection.

## * The best way to use the RNN model:
  * Locally run RNN inference within ROOT/C++ utilizing the LibTorch library.

  * The `"rnn_inference"` folder contains all of the scripts needed to do this.

  * More information on how exactly to set this up is in this guide document:
    * https://docs.google.com/document/d/1PF5z-vx-rHSHfrryYS5lpSa-mMfgve9uZSoTvvZcdKU/edit?usp=sharing
  * This document also explains how to input blips into the BDT which gives you a score based on how signal or background-like the blip is.

* P.S.
  * In `"analyze_neutrons/accessory_macros"` I started making a `"Visualizer3D.C"` script that uses true particle variables to try to create EVD-like 3D vizualizations of certain interactions.
  * It would be very nice to have this working so that you can easily investigate anomalous interaction cases, but I never had the time to finish it, you who inherits this analysis might want to give it a try!

# GOOD LUCK !
