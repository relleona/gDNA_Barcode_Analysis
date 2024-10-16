'''
FASTQ File Reorganization Script when SampleSheet and Primers containing Stagger Length are provided

This script reorganizes sequencing data files (FASTQ) based on sample names,
grouping them by their base name (e.g., SA2-039, SA2-040), creates a comprehensive StagerFileAll.csv,
and individual StaggerFile.csv files for each base sample group.

Required Packages: os, re, shutil, argparse, csv, pandas

Expected Input:
- A directory containing nested folders with FASTQ files (.fastq.gz).
- File naming pattern: SampleName_Number_Description_S\d+_L00\d_R\d_001.fastq.gz
- A sample sheet Excel file (e.g., SAUpdated_GoyalPool3_SampleSheet_Apr2023.xlsx)
- A primers Excel file (e.g., 20200811_10X_BCSeq_Primers.xlsx)

Expected Output:
- New directories created for each base sample name (e.g., SA2-039, SA2-040).
- Subdirectories for each full sample name within the base sample directories.
- FASTQ files copied into their respective sample subdirectories.
- A comprehensive StagerFileAll.csv in the main output directory.
- Individual StaggerFile.csv files in each base sample group directory.
- Console output indicating which files were copied and any errors encountered.
- New directories to create extData or extracted data 

Example:
input_path/
├── Subfolder1/
│   ├── Folder1/
│   │   └── SA2-039_3_GemFull_S2_L001_R1_001.fastq.gz
│   ├── Folder2/
│   │   └── SA2-039_2_GemFull_S2_L001_R1_001.fastq.gz
│   ├── Folder3/
│   │   └── SA2-040_3_SotFull_S4_L001_R1_001.fastq.gz
│   └── Folder4/
│       └── SA2-040_2_SotFull_S4_L001_R1_001.fastq.gz
└── Subfolder2/
    ├── Folder5/
    │   └── SA2-041_1_TerFull_S3_L001_R1_001.fastq.gz
    ├── Folder6/
    │   └── SA2-041_2_TerFull_S3_L001_R1_001.fastq.gz
    ├── Folder7/
    │   └── SA2-042_4_QuaFull_S1_L001_R1_001.fastq.gz
    └── Folder8/
        └── SA2-042_3_QuaFull_S1_L001_R1_001.fastq.gz

output_path/
├── StagerFileAll.csv
├── SA2-039/
│   ├── StaggerFile.csv
│   └── rawData/
│       ├── SA2-039_3_GemFull/
│       │   └── SA2-039_3_GemFull_S2_L001_R1_001.fastq.gz
│       └── SA2-039_2_GemFull/
│           └── SA2-039_2_GemFull_S2_L001_R1_001.fastq.gz
├── SA2-040/
│   ├── StaggerFile.csv
│   └── rawData/
│       ├── SA2-040_3_SotFull/
│       │   └── SA2-040_3_SotFull_S4_L001_R1_001.fastq.gz
│       └── SA2-040_2_SotFull/
│           └── SA2-040_2_SotFull_S4_L001_R1_001.fastq.gz
├── SA2-041/
│   ├── StaggerFile.csv
│   └── rawData/
│       ├── SA2-041_1_TerFull/
│       │   └── SA2-041_1_TerFull_S3_L001_R1_001.fastq.gz
│       └── SA2-041_2_TerFull/
│           └── SA2-041_2_TerFull_S3_L001_R1_001.fastq.gz
└── SA2-042/
    ├── StaggerFile.csv
    └── rawData/
        ├── SA2-042_4_QuaFull/
        │   └── SA2-042_4_QuaFull_S1_L001_R1_001.fastq.gz
        └── SA2-042_3_QuaFull/
            └── SA2-042_3_QuaFull_S1_L001_R1_001.fastq.gz

Command to run this file:
python3 <path to BaseSpaceReorganization.py> <path/to/nested/folders> --sample_sheet <path/to/sample_sheet.xlsx> --primers_sheet <path/to/primers_sheet.xlsx> [-exp <path/to/custom/output>] 

Note: If -raw is not specified, output will be in a 'reorganized_data' folder within the input directory.

StagerFileAll.csv Content:
Sample,Path,Index2,Stagger_Length
SA2-039_3_GemFull,/path/to/output/SA2-039/SA2-039_3_GemFull,TAGATCGC,5
SA2-039_2_GemFull,/path/to/output/SA2-039/SA2-039_2_GemFull,TAGATCGC,4
...

StaggerFile.csv Content (e.g., in SA2-039 folder):
SA2-039_3_GemFull,5
SA2-039_2_GemFull,4
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''

# Import necessary libraries
import os  # For file and directory operations
import re  # For regular expressions
import shutil  # For high-level file operations
import csv  # For reading and writing CSV files
import pandas as pd  # For handling Excel files
import openpyxl  # Required for pandas to read Excel files
from argparse import ArgumentParser  # For parsing command-line arguments
from Reorganizationfunction import process_excel_files, extract_primer_data, \
    extract_sample_data, create_stagger_file, \
    create_stager_file_all, extract_sample_info, \
    sanitize_dirname, group_samples, find_common_prefix, write_main_folders_dict, \
    get_sample_name, combine_fastq_files

def main(input_path, output_path, sample_sheet, primers_sheet, script_path):
    """Main function to reorganize FASTQ files and create necessary CSV files."""
    # Process Excel files to get combined sample data
    sample_data = process_excel_files(sample_sheet, primers_sheet)
    # Create output directory if it doesn't exist
    os.makedirs(output_path, exist_ok=True)

    # Create StagerFileAll.csv first
    stager_file_path = create_stager_file_all(output_path, sample_data)
    
    # Read StagerFileAll.csv to get all sample IDs
    with open(stager_file_path, 'r') as f:
        reader = csv.DictReader(f)
        sample_ids = [row['Sample'] for row in reader]

    # Dictionary to store main folders and their sample IDs
    main_folders = {}

    # Iterate through the directory structure
    for root, _, files in os.walk(input_path):
        for file in files:
            # Process only .fastq.gz files
            if file.endswith('.fastq.gz'):
                file_path = os.path.join(root, file)
                for sample_id in sample_ids:
                    if sample_id in file:
                        # Extract the main folder name (e.g., SA2-037)
                        main_folder = '-'.join(sample_id.split('_')[0].split('-')[:2])
                        
                        # Create full path for the main folder
                        main_folder_path = os.path.join(output_path, main_folder)
                        
                        # Add main folder path and sample ID to the dictionary
                        if main_folder_path not in main_folders:
                            main_folders[main_folder_path] = []
                        if sample_id not in main_folders[main_folder_path]:
                            main_folders[main_folder_path].append(sample_id)

                        # Create new directory structure and copy file
                        new_sample_dir = os.path.join(main_folder_path, "raw", sample_id)
                        os.makedirs(new_sample_dir, exist_ok=True)
                        
                        destination = os.path.join(new_sample_dir, file)
                        try:
                            shutil.copy2(file_path, destination)
                            print(f"Copied {file_path} to {destination}")
                        except shutil.Error as e:
                            print(f"Error copying {file_path}: {e}")
                        break  # Stop searching for matching sample_ids once found

    for main_folder_path, sample_id_list in main_folders.items():
        # Create StaggerFile.csv for each main folder
        create_stagger_file(os.path.basename(main_folder_path), stager_file_path, output_path)
        staggerfilepath=os.path.join(main_folder_path, "StaggerFile.csv")
        # destinationfile = os.path.join(main_folder_path, "templateInputFullRun.py")
        # modify_config(sourcefile, destinationfile, main_folder_path,staggerfilepath, sample_id_list)

    print("File reorganization and CSV file creation complete.")
    print(f"Reorganized files and CSV files are located in: {output_path}")
    
    # Add on to the txt file that contains all the necessary inputs for the bash script in the different experiments 
    os.chdir(output_path)
    # Create txt file in the Experiment folder 
    output_file = "Templateinput.txt"
    write_main_folders_dict(output_file, main_folders)


if __name__ == "__main__":
    # Set up command-line argument parser
    parser = ArgumentParser(description="Reorganize FASTQ files and create StagerFileAll.csv and StaggerFile.csv files.")
    parser.add_argument("path", help="Specify the path containing the nested sample folders.", type=str)
    parser.add_argument("--sample_sheet", help="Path to the sample sheet Excel file", required=True)
    parser.add_argument("--primers_sheet", help="Path to the primers Excel file", required=True)
    parser.add_argument("-spath", "--scriptpath", help="Path to all script files.",required=True)
    parser.add_argument("-exp", "--experimentfile", help="Specify the output path. If not provided, a 'Experiment' folder will be created in the input path.", type=str, required=False)
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Determine the output directory
    output_dir = args.experimentfile if args.experimentfile else os.path.join(args.path, "Experiment")
    
    # Run the main function with the provided arguments
    main(args.path, output_dir, args.sample_sheet, args.primers_sheet, args.script_path)
