#!/bin/bash

#This script builds the conda envs and fixes in error within the rmats_nextflow env

conda env create -f ~/nextflow-rnaseq/envs/R_nextflow_rnaseq.yml

conda env create -f ~/nextflow-rnaseq/envs/zlib_nextflow_rnaseq.yml

conda env create -f ~/nextflow-rnaseq/envs/rmats_nextflow.yml

mv ~/.conda/envs/rmats_nextflow/lib/libncurses.so.6 ~/.conda/envs/rmats_nextflow/lib/libncurses.so.5
