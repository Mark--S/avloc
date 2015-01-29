#!/bin/bash

for file in outputs/Fibre_FT0* 
do
    echo $file
    bin/make_ntuple $file

done
