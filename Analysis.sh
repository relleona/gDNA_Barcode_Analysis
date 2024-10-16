#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=03:00:00
#SBATCH --job-name=LV_30
#SBATCH --output=outlog
#SBATCH --error=errlog


module purge
eval "$(conda shell.bash hook)"
conda activate Barcode_extraction

#Code to run the overall python script for the pipeline: This is the file where you modified all the inputs.
python3 /projects/b1042/GoyalLab/aleona/gDNA_Barcode_Analysis/Inputs.py

