#!/bin/bash

# Before running this, read step0.sh, step1.sh, step2.sh and step3.sh so that all details can be checked, especially if running the pipeline for a new experiment set or on a new system

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=03:00:00
#SBATCH --job-name=LV_30
#SBATCH --output=outlog
#SBATCH --error=errlog


# module purge
# eval "$(conda shell.bash hook)"
# conda activate Barcode_extraction
#Code to run the overall python script for the pipeline: This is the file where you modified all the inputs.
python3 /home/keerthana/Goyal_Lab/FateMapPipeline_Ubuntu/gDNA_Barcode_Extraction/experimentIDInputRun.py

