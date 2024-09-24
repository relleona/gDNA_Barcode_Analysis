'''
FASTQ File Reorganization Script with Excel Data Integration and CSV Generation

This script reorganizes sequencing data files (FASTQ) based on sample names,
grouping them by their base name (e.g., SA2-039, SA2-040), creates a comprehensive StagerFileAll.csv,
and individual StaggerFile.csv files for each base sample group.

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
    stager_file_all_path = os.path.join(output_dir, "StagerFileAll.csv")
    
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

def main(input_path, output_path, sample_sheet, primers_sheet):
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

    # Set to store main folders
    main_folders = set()

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
                        main_folders.add(main_folder)

                        # Create new directory structure and copy file
                        new_sample_dir = os.path.join(output_path, main_folder, "raw", sample_id)
                        os.makedirs(new_sample_dir, exist_ok=True)
                        
                        destination = os.path.join(new_sample_dir, file)
                        try:
                            shutil.copy2(file_path, destination)
                            print(f"Copied {file_path} to {destination}")
                        except shutil.Error as e:
                            print(f"Error copying {file_path}: {e}")
                        break  # Stop searching for matching sample_ids once found

                    
   
    for main_folder in main_folders:
        # Create StaggerFile.csv for each main folder
        create_stagger_file(main_folder, output_path)
        

    print("File reorganization and CSV file creation complete.")
    print(f"Reorganized files and CSV files are located in: {output_path}")

if __name__ == "__main__":
    # Set up command-line argument parser
    parser = ArgumentParser(description="Reorganize FASTQ files and create StagerFileAll.csv and StaggerFile.csv files.")
    parser.add_argument("path", help="Specify the path containing the nested sample folders.", type=str)
    parser.add_argument("--sample_sheet", help="Path to the sample sheet Excel file", required=True)
    parser.add_argument("--primers_sheet", help="Path to the primers Excel file", required=True)
    parser.add_argument("-exp", "--experimentfile", help="Specify the output path. If not provided, a 'Experiment' folder will be created in the input path.", type=str, required=False)

    # Parse the command-line arguments
    args = parser.parse_args()

    # Determine the output directory
    output_dir = args.experimentfile if args.experimentfile else os.path.join(args.path, "Experiment")
    
    # Run the main function with the provided arguments
    main(args.path, output_dir, args.sample_sheet, args.primers_sheet)