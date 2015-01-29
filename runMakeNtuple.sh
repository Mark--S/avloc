#!/bin/bash
source /mnt/lustre/epp_scratch/neutrino/rat/snoing/install/env_rat-5.0.0.sh
for file in outputs/Fibre_FT0* 
do
    echo $file
    bin/make_ntuple $file

done
