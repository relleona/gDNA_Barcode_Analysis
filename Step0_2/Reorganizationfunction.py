
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
            writer.writerow([str_sample_id, os.path.join(output_dir, main_folder, str_sample_id), index2, stagger_length])
    # Print confirmation message
    print(f"Created StagerFileAll.csv in {output_dir}")
    # Return the path of the created file
    return stager_path

def create_stagger_file(main_folder, output_dir):
    """
    Create a StaggerFile.csv for each main folder by reading from StagerFileAll.csv
    and checking the actual subfolders present in the main folder.
    """
    # Path to StagerFileAll.csv
    stager_file_all_path = os.path.join(output_dir, "StaggerFileAll.csv")
    
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
    

def modify_config(source_file, destination_file, new_path_experiment_folder, 
                  new_stagger_file, new_sample_array):
    """
    Modify specific parameters in a Python configuration file.
    
    :param source_file: Path to the original configuration file
    :param destination_file: Path where the modified configuration will be saved
    :param new_path_experiment_folder: New value for pathExperimentFolder
    :param new_stagger_file: New value for staggerFile
    :param new_sample_array: New list for sampleArray and sampleArrayStarcode
    """
    
    try:
        # Check if the source file exists
        if not os.path.exists(source_file):
            raise FileNotFoundError(f"Source file not found: {source_file}")

        # Read the content of the source file
        with open(source_file, 'r') as file:
            content = file.read()

        print(f"Original content length: {len(content)}")

        # Perform string replacements
        new_content = content

        new_content = re.sub(r'pathExperimentFolder\s*=\s*["\'].*?["\']', 
                             f'pathExperimentFolder = "{new_path_experiment_folder}"', new_content)
        
        new_content = re.sub(r'staggerFile\s*=\s*["\'].*?["\']', 
                             f'staggerFile = "{new_stagger_file}"', new_content)
        
        sample_array_str = str(new_sample_array).replace("'", '"')
        new_content = re.sub(r'sampleArray\s*=\s*\[.*?\]', 
                             f'sampleArray = {sample_array_str}', new_content, flags=re.DOTALL)
        new_content = re.sub(r'sampleArrayStarcode\s*=\s*\[.*?\]', 
                             f'sampleArrayStarcode = {sample_array_str}', new_content, flags=re.DOTALL)

        print(f"New content length: {len(new_content)}")

        if content == new_content:
            print("Warning: No changes were made to the content.")
        else:
            print("Changes detected in the content.")

        # Ensure the directory for the destination file exists
        os.makedirs(os.path.dirname(destination_file), exist_ok=True)

        # Write the updated content to the destination file
        with open(destination_file, 'w') as file:
            file.write(new_content)

        print(f"Configuration successfully updated and saved to: {destination_file}")
        
        # Verify the file was written
        if os.path.exists(destination_file):
            with open(destination_file, 'r') as file:
                final_content = file.read()
            print(f"Final file content length: {len(final_content)}")
            if final_content == new_content:
                print("File content matches expected content.")
            else:
                print("Warning: File content does not match expected content.")
        else:
            print(f"Warning: Destination file {destination_file} does not exist after write operation.")

    except FileNotFoundError as e:
        print(f"Error: {e}")
    except PermissionError:
        print(f"Error: Permission denied when trying to write to {destination_file}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")