#!/bin/bash

ProcessNTuples /Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_blipPlus_20260304.root Blip_n_Np overlay_out.root > overlay.log 2>&1 &
ProcessNTuples /Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_09_BNB_beamoffEXT_surprise_blipPlus_20260522.root Blip_n_Np EXT_out.root     > EXT.log     2>&1 &
ProcessNTuples /Users/liani/nu_work/pelee/files/MCC9p10_Run4b_v10_04_07_11_BNB_opendata_surprise_blipPlus_20260522.root Blip_n_Np openData_out.root > openData.log 2>&1 &

wait
echo "All done"
