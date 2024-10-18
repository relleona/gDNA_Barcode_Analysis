
'''
FASTQ File Reorganization Script with Excel Data Integration and CSV Generation

This script contains functions needed to reorganize sequencing data files (FASTQ) based on sample names,
grouping them by their base name (e.g., SA2-039, SA2-040), creates a comprehensive StagerFileAll.csv, 
individual StaggerFile.csv files for each base sample group and 
create an TemplateInput copies for all Experiment folders

Required Packages: os, re, shutil, argparse, csv, pandas

Key Features:
1. Handles deeply nested folder structures (two levels of subfolders).
2. Extracts sample names and base names from FASTQ files.
3. Processes sample sheet and primers Excel files for additional data.
4. Creates a two-level directory structure: base name -> full sample name.
5. Copies files to maintain the original data integrity.
6. Creates a comprehensive StagerFileAll.csv in the main output directory.
7. Creates individual StaggerFile.csv files for each base sample group.
8. Provides flexibility in output location (custom path or within input directory).
9. Handles potential file system errors and provides informative output.

'''
# Import necessary libraries
import os  # For file and directory operations
import ast #For .py file config
import re  # For regular expressions
import shutil  # For high-level file operations
import csv  # For reading and writing CSV files
import pandas as pd  # For handling Excel files
import openpyxl  # Required for pandas to read Excel files
import numpy as np 
from argparse import ArgumentParser  # For parsing command-line arguments
from collections import defaultdict

def sanitize_dirname(name):
    """Sanitize directory names by removing or replacing invalid characters."""
    # Remove or replace invalid characters with underscores
    sanitized = re.sub(r'[^\w\-_\. ]', '_', name)
    # Remove leading/trailing whitespace and carriage returns
    return sanitized.strip().rstrip('\r')

def extract_sample_info(filename, sample_ids):
    """
    Extract sample info based on the filename and known sample IDs.
    
    Args:
    filename (str): The filename to extract information from.
    sample_ids (list): List of known sample IDs from StagerFileAll.csv.
    
    Returns:
    tuple: A tuple containing (main_folder, sample_id) or (None, None) if no match.
    """
    # Iterate through all known sample IDs
    for sample_id in sample_ids:
        # If the sample ID is found in the filename
        if sample_id in filename:
            # Create the main folder name from the first two sections of the sample ID
            main_folder = '_'.join(sample_id.split('_')[:2])
            return main_folder, sample_id
    # If no match is found, return None for both values
    return None, None


def create_stager_file_all(output_dir, sample_data):
    """Create a comprehensive StagerFileAll.csv containing information for all samples."""
    # Define the path for StagerFileAll.csv
    stager_path = os.path.join(output_dir, "StagerFileAll.csv")
    # Open the file for writing
    with open(stager_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write the header row
        writer.writerow(['Sample', 'Path', 'Index2', 'Stagger_Length'])
        # Iterate through all samples and write their data
        for sample_id, (index2, stagger_length) in sample_data.items():
            # Convert sample_id to string
            str_sample_id = str(sample_id)
            # Create main_folder from sample_id
            if '_' in str_sample_id:
                main_folder = '_'.join(str_sample_id.split('_')[:2])
            else:
                # Fallback for sample IDs without underscores
                main_folder = str_sample_id
            # Write the row for this sample
            writer.writerow([str_sample_id, os.path.join(output_dir, main_folder, "raw", str_sample_id), index2, stagger_length])
    # Print confirmation message
    print(f"Created StagerFileAll.csv in {output_dir}")
    # Return the path of the created file
    return stager_path

def create_stagger_file(main_folder, stager_file_all_path, output_dir):
    """
    Create a StaggerFile.csv for each main folder by reading from StagerFileAll.csv
    and checking the actual subfolders present in the main folder.
    """
    # Path for the new StaggerFile.csv
    stagger_path = os.path.join(output_dir, main_folder, "StaggerFile.csv")
    
    # Path to the main folder
    raw_data_path = os.path.join(output_dir, main_folder, "raw")
    
    # Read StagerFileAll.csv into a dictionary for quick lookup
    stager_data = {}
    with open(stager_file_all_path, 'r') as stager_file_all:
        reader = csv.DictReader(stager_file_all)
        for row in reader:
            stager_data[row['Sample']] = row['Stagger_Length']
    
    # Get the list of actual subfolders in the main folder
    actual_subfolders = [f for f in os.listdir(raw_data_path) if os.path.isdir(os.path.join(raw_data_path, f))]
    
    # Write the filtered and verified data to StaggerFile.csv
    with open(stagger_path, 'w', newline='') as stagger_file:
        writer = csv.writer(stagger_file)
        for subfolder in actual_subfolders:
            if subfolder in stager_data:
                writer.writerow([subfolder, stager_data[subfolder]])
            else:
                print(f"Warning: {subfolder} found in directory but not in StagerFileAll.csv")
                writer.writerow([subfolder, 'N/A'])

    print(f"Created StaggerFile.csv for {main_folder}")

def extract_sample_data(file_path):
    """Extract sample data from the sample sheet Excel file."""
    try:
        # Read the entire Excel file
        df = pd.read_excel(file_path, header=None)
        
        # Find the row indices for '[BCLConvert_Data]' and '[Cloud_Settings]'
        start_row = df[df[0] == '[BCLConvert_Data]'].index[0]
        end_row = df[df[0] == '[Cloud_Settings]'].index[0]
        
        # Extract data between '[BCLConvert_Data]' and '[Cloud_Settings]'
        data = df.iloc[start_row+1:end_row].reset_index(drop=True)
        
        # Set the first row as column names
        data.columns = data.iloc[0]
        data = data.iloc[1:].reset_index(drop=True)

        # Print all column names for debugging
        print("Columns in the sample sheet Excel file:")
        for col in data.columns:
            print(f"  - '{col}' (type: {type(col)})")

        # Check if we have the required columns
        required_columns = ['Sample_ID', 'Index2']
        missing_columns = [col for col in required_columns if col not in data.columns]
        if missing_columns:
            print(f"Error: The following required columns are missing: {', '.join(missing_columns)}")
            print("Please ensure the Excel file contains 'Sample_ID' and 'Index2' columns.")
            raise ValueError("Missing required columns in Excel file.")

        # Convert Sample_ID to string and create dictionary
        return data.astype({'Sample_ID': str}).set_index('Sample_ID')['Index2'].to_dict()
    
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        print("Please ensure the Excel file is not empty and is correctly formatted.")
        raise

def extract_primer_data(file_path):
    """Extract primer data from the primers Excel file."""
    try:
        # Read the first sheet of the Excel file, using the first row as header
        df = pd.read_excel(file_path, sheet_name=0, header=0)
        
        # Print all column names for debugging
        print("Columns in the primers Excel file:")
        for col in df.columns:
            print(f"  - '{col}' (type: {type(col)})")

        # Check if we have the required columns
        required_columns = ['Index', 'Stagger_Length']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            print(f"Error: The following required columns are missing: {', '.join(missing_columns)}")
            print("Please ensure the Excel file contains 'Index' and 'Stagger_Length' columns.")
            
            # Try to find columns with similar names
            for missing_col in missing_columns:
                similar_cols = [col for col in df.columns if missing_col.lower() in str(col).lower()]
                if similar_cols:
                    print(f"Similar column names found for '{missing_col}': {', '.join(map(str, similar_cols))}")
            
            raise ValueError("Missing required columns in Excel file.")

        # Convert Stagger_Length to integer, replacing NaN with -1
        df['Stagger_Length'] = pd.to_numeric(df['Stagger_Length'], errors='coerce').fillna(-1).astype(int)

        # Print unique values in Stagger_Length for verification
        print("Unique values in Stagger_Length after conversion:")
        print(df['Stagger_Length'].unique())

        # Return dictionary of Index and Stagger_Length pairs
        return df.set_index('Index')['Stagger_Length'].to_dict()
    
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        print("Please ensure the Excel file is not empty and is correctly formatted.")
        raise

def process_excel_files(sample_sheet, primers_sheet):
    """Process both Excel files and combine their data."""
    # Extract data from both Excel files
    sample_data = extract_sample_data(sample_sheet)
    primer_data = extract_primer_data(primers_sheet)
    
    # Combine data from both sources
    combined_data = {}
    for sample_id, index2 in sample_data.items():
        # Get the stagger length for this index (default to 'N/A' if not found)
        stagger_length = primer_data.get(index2, 'N/A')
        combined_data[sample_id] = (index2, stagger_length)
    
    return combined_data

def find_common_prefix(sample_ids):
    if not sample_ids:
        return ""
    
    # Start with the first sample ID as the prefix
    prefix = sample_ids[0]
    
    for sample_id in sample_ids[1:]:
        # Find the common part between the current prefix and the new sample ID
        common = ""
        for i, (c1, c2) in enumerate(zip(prefix, sample_id)):
            if c1 != c2:
                break
            common += c1
        
        # Update the prefix to be the common part
        prefix = common
        
        # If we've reduced the prefix to nothing, stop
        if not prefix:
            break
    
    # Remove any trailing non-alphanumeric characters
    prefix = re.sub(r'[^a-zA-Z0-9]+$', '', prefix)
    
    return prefix

def combine_fastq_files(files, output_file):
    with open(output_file, 'wb') as outfile:
        for fname in sorted(files):
            with open(fname, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

def group_samples(sample_ids):
    # Initialize dictionaries to store groups and reverse mapping
    groups = defaultdict(list)
    sample_to_group = {}
    
    # Find the common prefix among all sample IDs
    common_prefix = find_common_prefix(sample_ids)
    prefix_length = len(common_prefix)
    
    # If common prefix is very short, try to extend it
    if prefix_length < 3:
        pattern = re.compile(rf'^({re.escape(common_prefix)}[A-Za-z0-9-]+)')
    else:
        pattern = re.compile(rf'^({re.escape(common_prefix)}[A-Za-z0-9]*)')
    
    # Group each sample ID and create reverse mapping
    for sample_id in sample_ids:
        match = pattern.match(sample_id)
        if match:
            group = match.group(1)
            groups[group].append(sample_id)
            sample_to_group[sample_id] = group
        else:
            # Fallback: use the whole sample ID as its own group
            groups[sample_id].append(sample_id)
            sample_to_group[sample_id] = sample_id
    
    return groups, sample_to_group

# Function to get sample name from filename
def get_sample_name(filename, group_samples):
    # filename is the name of the file we are trying to match to the samplenames/subfolder
    # group_samples is the name of the samplenames/subfolder

    # Return the first sample name that matches the filename
    # Next returns the second argument, which is None if the first value is empty, or there is no match
    return next((sample for sample in group_samples if sample in filename), None)

# Printing this to a txt file 
def write_main_folders_dict(output_file, main_folder_dict):
    # print(main_folder_dict)
    with open(output_file, "a") as file:
        for folder_path, samples in main_folder_dict.items():
            # Join samples with comma and space, ensure that there are no duplicates 
            samples_string = ",".join(np.unique(samples))
            # Write the folder path and samples to the file
            file.write(f"{folder_path} = {samples_string}\n")
