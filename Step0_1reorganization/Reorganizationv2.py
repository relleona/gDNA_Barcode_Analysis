'''
FASTQ File Reorganization Script when only StaggerFileAll.csv is provided

This script reorganizes sequencing data files (FASTQ) based on sample names,
grouping them by their base name (e.g., SA2-039, SA2-040), creates a comprehensive StagerFileAll.csv,
and individual StaggerFile.csv files for each base sample group.

Required Packages: os, re, shutil, argparse, csv, pandas

'''
# Import necessary libraries
import os  # For file and directory operations
import re  # For regular expressions
import shutil  # For high-level file operations
import csv  # For reading and writing CSV files
import pandas as pd  # For handling Excel files
from collections import defaultdict
import openpyxl  # Required for pandas to read Excel files
from argparse import ArgumentParser  # For parsing command-line arguments
from Reorganizationfunction import process_excel_files, extract_primer_data, \
    extract_sample_data, create_stagger_file, \
    create_stager_file_all, extract_sample_info, \
    sanitize_dirname, group_samples, find_common_prefix, write_main_folders_dict, \
    get_sample_name, combine_fastq_files


def main(input_path, output_path, stagger_file_path, script_path):
    """Main function to reorganize FASTQ files and create necessary CSV files."""

    # Print starting message with input parameters
    print(f"Starting main function with inputs: {input_path}, {output_path}, {stagger_file_path}, {script_path}")

    ## Open and read the StaggerFileAll.csv file
    with open(stagger_file_path, 'r') as f:
        reader = csv.DictReader(f)
        # Trim whitespace from fieldnames
        reader.fieldnames = [field.strip() for field in reader.fieldnames]

        # Find the correct column name for "Sample" (case-insensitive)
        sample_column = next(
            (field for field in reader.fieldnames if field.lower() == 'sample'), None)

        # Raise an error if the 'Sample' column is not found
        if sample_column is None:
            raise ValueError(
                f"'Sample' column not found in {stagger_file_path}. Available columns: {', '.join(reader.fieldnames)}")

        # Extract sample IDs from the 'Sample' column
        sample_ids = [row[sample_column] for row in reader]


    # Group the samples and get the reverse mapping function
    sample_groups, sample_to_group = group_samples(sample_ids)

    # Print the number of groups
    print(f"Number of groups: {len(sample_groups)}")

    # Print the first 5 groups and their samples
    for group, samples in list(sample_groups.items())[:5]:
        print(f"Group {group}: {samples}")

    # Initialize variables for file counting and matching
    file_count = 0
    match_count = 0
    main_folders = defaultdict(list)  # Use defaultdict instead of regular dict


    # Walk through the input directory
    for root, _, files in os.walk(input_path):
        # Filter for .fastq.gz files
        # Find if in that folder sets there are files of the format
        fastq_files = [f for f in files if f.endswith('.fastq.gz')]
        
        # Process each fastq file
        for file in fastq_files:
            # Get the full file path
            file_path = os.path.join(root, file)
            
            # Find the matching sample and group
            # sample_to_group.keys() is the subfolder under raw 
            # sample_to_group.values() is the mainfolder
            sample_name = get_sample_name(file, sample_to_group.keys())
            
            # If a matching sample is found, process the file
            if sample_name:
                # Get the group for this sample
                group = sample_to_group[sample_name]
                
                # Update counters
                match_count += 1
                
                # Create main folder path and add to dictionary
                main_folder_path = os.path.join(output_path, group)
                main_folders[main_folder_path].append(sample_name)
                
                # Create new directory structure for the sample
                sample_dir = os.path.join(main_folder_path, "raw", sample_name)
                
                # Create the sample directory
                os.makedirs(sample_dir, exist_ok=True)
                
                # Set the destination path for the file
                destination = os.path.join(sample_dir, file)
                
                # Copy the file to its new location
                try:
                    shutil.copy2(file_path, destination)
                    print(f"Copied {file_path} to {destination}")
                    file_count += 1
                except Exception as e:
                    print(f"Error copying {file_path}: {e}")
            else:
                # Print a message if no matching sample is found
                print(f"Could not determine sample for file: {file}")

    # Print summary of processed files
    print(f"Total files processed: {file_count}")
    print(f"Found {match_count} matching files")
    print(f"Processed {len(main_folders)} main folders")

    # Define the source file path for the config template
    sourcefile = os.path.join(script_path, "templateInputFullRun.py")

    # Iterate through main folders and their sample IDs
    for main_folder_path, sample_id_list in main_folders.items():
        print(f"Modifying config for {main_folder_path} with {len(sample_id_list)} samples")
        
        # Create StaggerFile.csv for each main folder
        create_stagger_file(os.path.basename(main_folder_path), stagger_file_path, output_path)
        staggerfilepath = os.path.join(main_folder_path, "StaggerFile.csv")

        # # Modify the config file using the imported function
        # modify_config(
        #     sourcefile,
        #     destinationfile,
        #     main_folder_path,
        #     staggerfilepath,
        #     sample_id_list)

    # Print completion message
    print("File reorganization and CSV file creation complete.")
    print(f"Reorganized files and CSV files are located in: {output_path}")
    
    # Add on to the txt file that contains all the necessary inputs for the bash script in the different experiments 
    os.chdir(output_path)
    # Create txt file in the Experiment folder 
    output_file = "Templateinput.txt"
    write_main_folders_dict(output_file, main_folders)



# Check if the script is being run directly
if __name__ == "__main__":
    # Set up command-line argument parser
    parser = ArgumentParser(
        description="Reorganize FASTQ files and create StagerFileAll.csv and StaggerFile.csv files.")
    parser.add_argument(
        "path",
        help="Specify the path containing the nested sample folders.",
        type=str)
    parser.add_argument(
        "-stagger",
        "--staggerfileall",
        help="Path to StaggerFileAll.csv .",
        required=True)
    parser.add_argument(
        "-spath",
        "--scriptpath",
        help="Path to all script files.",
        required=True)
    parser.add_argument(
        "-exp",
        "--experimentfile",
        help="Specify the output path. If not provided, a 'Experiment' folder will be created in the input path.",
        type=str,
        required=False)

    # Parse the command-line arguments
    args = parser.parse_args()

    # Determine the output directory
    output_dir = args.experimentfile if args.experimentfile else os.path.join(
        args.path, "Experiment")

    # Run the main function with the provided arguments
    main(args.path, output_dir, args.staggerfileall, args.scriptpath)

