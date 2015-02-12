#!/bin/bash
for file in outputs/Fibre_FT0* 
do
    qsub -cwd -q mps.q ./submit.sh $file 
done
