#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-ubuntu2404-gcc13-opt/setup.sh
python3 bdt.py --signal_root "/work/ambient_generation/examples/c++/mixed_signal.root" \
       --bkg_root "/work/ambient_generation/examples/c++/ambient_out.root" \
       --n_signal 0 \
       --n_bkg 0 \
       --scale_pos_weight 10 \
       --max_len 90 \
       --exclude_branches "t_min" "t_mean" \
       --batch_size 10000
