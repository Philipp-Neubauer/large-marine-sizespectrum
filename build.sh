#!/bin/bash
set -ex 

# add scripts to run here
cd analysis/Barents\ Sea && Rscript BarentsSea_minimal.R
# copy outputs to Gorbachev outputs folder
cp -r plots/*.pdf /output/