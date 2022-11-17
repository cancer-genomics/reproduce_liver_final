#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_fsize=10G
#$ -l mem_free=10G
#$ -l h_vmem=10G
#$ -l h_rt=144:00:00

module load conda_R
#Rscript run_screen.r
Rscript run_screen_delfi.r