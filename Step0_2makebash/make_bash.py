
''' 
Genome Index Generation Script for STAR Aligner

Note about the script:
This script generates a job submission file for creating a genome index using STAR aligner.
It's designed to work with high-performance computing environments, particularly those using SLURM
for job scheduling.

Required Packages: datetime, time, os, sys, numpy, subprocess, csv, argparse, glob and multiprocessing

Key Features:
1. Generates a SLURM job submission script for STAR genome indexing.
2. Supports both FASTA and GFF genome formats.
3. Configurable for different HPC environments (Quest-specific options included).
4. Creates necessary directory structure for job files, logs, and genome index.
5. Handles large genomes with many scaffolds by adjusting genomeChrBinNbits.

Script Output:
- A SLURM job submission script for STAR genome indexing
- Directory structure for jobs, logs, and genome index
- Console output indicating job submission status

STAR Genome Generation Output:
After the STAR genome generation process completes, the following files will be created in the specified genome index directory:
1. chrName.txt: Names of the chromosomes
2. chrLength.txt: Lengths of the chromosomes
3. chrStart.txt: Start positions of the chromosomes in the SAindex
4. chrNameLength.txt: Names and lengths of the chromosomes
5. genomeParameters.txt: Parameters of genome generation
6. Genome: Suffix array (SA) of the genome
7. SA: Suffix array offsets
8. SAindex: Index for the suffix array file
9. sjdbInfo.txt: Information about splice junctions
10. sjdbList.out.tab: List of splice junctions
11. sjdbList.fromGTF.out.tab: List of splice junctions extracted from the GTF file (if GTF was provided)


Note: This script is designed for use in a SLURM-managed HPC environment. Modify as needed for
different job scheduling systems. Always test on a small scale before running on full datasets.

Command to run this file: 
python3 <path to make_bash.py> <path/to/Experiment_folder> -d depth -q quest -a account -p partition -n nodes -nt ntasks -t time -m memory -e email -mt mail_type

- created by Aurelia on 2024-08-26 and last modified on 2024-08-27 at 11:21
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
from datetime import date, datetime, timedelta
import os
import argparse
import subprocess
import regex as re

def file_to_dict(file_path):
    """
    Read a text file and convert its contents to a dictionary.
    Each line of the form KEY=VALUE becomes a dictionary entry.
    
    :param file_path: Path to the text file
    :return: Dictionary with keys and values from the file
    """
    result_dict = {}

    try:
        with open(file_path, 'r') as file:
            for line in file:
                line = line.strip()
                # Match lines of the form KEY=VALUE or KEY="VALUE"
                match = re.match(r'^(.+?)\s*=\s*(.+)$', line)
                if match:
                    key, value = match.groups()
                    result_dict[key] = value

        if not result_dict:
            print("No valid key-value pairs found in the file.")
        
        return result_dict
    
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return None
    except IOError:
        print(f"Error: Could not read file '{file_path}'.")
        return None


def create_and_submit_slurm_script(output_path, subfolder, args):
    # Get current date and time
    today = date.today().strftime('%Y-%m-%d')
    current_time = datetime.now().strftime("%H.%M")
    
    # Initialize path 
    pathtosubfolder= os.path.join(output_path, subfolder)
    # Create a directory to store jobs and logs 
    job_directory = os.path.join(output_path, subfolder, "jobs")
    logs_directory = os.path.join(output_path, subfolder, "logs")
    os.makedirs(job_directory, exist_ok=True)
    os.makedirs(logs_directory, exist_ok=True)
    
    job_file = os.path.join(job_directory, f"job_{subfolder}_{today}_{current_time}.sh")
    
    # txt file that contains the inputs to the commands 
    input_file = os.path.join(args.output_path, "Templateinput.txt")

    # convert it back into a dictionary 
    result_dict = file_to_dict(input_file)
    
    print(result_dict)

    # Create job file name
    job_file = os.path.join(job_directory, f"run_barcode_analysis_{today}_{current_time}.sh")
    print(f"Creating job file: {job_file}")

    with open(job_file, 'w') as fh:
        # SLURM directives
        fh.write("#!/bin/bash\n")
        slurm_directives = [
            f"#SBATCH --account={args.account}",
            f"#SBATCH --nodes={args.nodes}",
            f"#SBATCH --partition={args.partition}",
            f"#SBATCH --ntasks={args.ntasks}",
            f"#SBATCH --time={args.job_time}",
            f"#SBATCH --mem={args.memory}",
            f"#SBATCH --mail-user={args.email}",
            f"#SBATCH --mail-type={args.mail_type}",
            f"#SBATCH --job-name=gDNApipeline_{subfolder}",
            f"#SBATCH --output={logs_directory}/%x_%j.out",
            f"#SBATCH --error={logs_directory}/%x_%j.err",
        ]
        fh.write('\n'.join(slurm_directives) + '\n\n')

        # Module and conda setup
        fh.write("module purge\n")
        fh.write("eval \"$(conda shell.bash hook)\"\n")
        fh.write("conda activate Barcode_extraction\n\n")

        # Step1 Variables
        fh.write(f"STEP1PATH={result_dict['STEP1PATH']}\n")
        fh.write(f"EXPERIMENT_DIR={pathtosubfolder}/\n")
        fh.write(f"SCRIPT_DIR={result_dict['SCRIPT_DIR']}/\n")
        fh.write(f"STAGGERFILEDIR={os.path.join(pathtosubfolder, 'StaggerFile.csv')}\n")
        fh.write(f"CHECKVECTOR={result_dict['CHECKVECTOR']}\n")
        fh.write(f"BARCODELENGTH={result_dict['BARCODELENGTH']}\n")
        fh.write(f"MINPHRED={result_dict['MINPHRED']}\n")
        fh.write(f"EXCLUDEREADS={result_dict['EXCLUDEREADS']}\n")
        fh.write(f"ASCIIOFFSET={result_dict['ASCIIOFFSET']}\n\n")

        # Step 1
        fh.write("# Step 1\n")
        step1_cmd = [
            'python3 "${STEP1PATH}" \\',
            '    "${EXPERIMENT_DIR}" \\',
            '    "${SCRIPT_DIR}" \\',
            '    --pathStaggerFile "${STAGGERFILEDIR}" \\',
            '    -r \\',
            '    -checkVector "${CHECKVECTOR}" \\',
            '    -barcodeLength "${BARCODELENGTH}" \\',
            '    -Q "${MINPHRED}" \\',
            '    -a "${ASCIIOFFSET}" \\',
            '    -e "${EXCLUDEREADS}"'
        ]
        fh.write('\n'.join(step1_cmd) + '\n\n')


        # Step2 Variables
        fh.write(f"STEP2PATH={result_dict['STEP2PATH']}\n")
        # The key to the sample array is the path to the subfolder 
        fh.write(f'SAMPLEARRAY="{result_dict[str(pathtosubfolder)]}"\n')
        fh.write(f"FRACTION={result_dict['FRACTION']}\n")
        fh.write(f"INPUTLENGTH={result_dict['INPUTLENGTH']}\n\n")

        # Step 2
        fh.write("# Step 2\n")
        step2_cmd = [
            'python3 "${STEP2PATH}" \\',
            '    "${SCRIPT_DIR}" \\',
            '    "${EXPERIMENT_DIR}" \\',
            '    "${SAMPLEARRAY}" \\',
            '    "${FRACTION}" \\',
            '    "${INPUTLENGTH}"'
        ]
        fh.write('\n'.join(step2_cmd) + '\n\n')
        
        # Step3 Variables
        fh.write(f"STEP3PATH={result_dict['STEP3PATH']}\n")
        # The key to the sample array is the path to the subfolder 
        fh.write(f"COMBINESAMPLE={result_dict['COMBINESAMPLE']}\n")
        fh.write(f"DISTANCE={result_dict['DISTANCE']}\n")
        if result_dict['COMBINESAMPLE'] == "yes":
            fh.write(f'STARCODEARRAY="{"Multiple_Samples"}"\n\n')
        else:
            fh.write(f'STARCODEARRAY="{result_dict[str(pathtosubfolder)]}"\n\n')

        # Step 3
        fh.write("# Step 3\n")
        step3_cmd = [
            'python3 "${STEP3PATH}" \\',
            '    "${SCRIPT_DIR}" \\',
            '    "${EXPERIMENT_DIR}" \\',
            '    "${COMBINESAMPLE}" \\',
            '    "${INPUTLENGTH}" \\',
            '    "${DISTANCE}" \\',
            f'    {args.ntasks} \\',
            '    "${STARCODEARRAY}" \\',
            '    "${FRACTION}"',
        ]
        fh.write('\n'.join(step3_cmd) + '\n')

    print(f"Job file created: {job_file}")
    # return job_file
    
    print(f"Created job script for {subfolder}: {job_file}")

    # Submit the job
    result = subprocess.run(["sbatch", job_file], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Job for {subfolder} submitted successfully")
        print("Output:", result.stdout)
    else:
        print(f"Job submission for {subfolder} failed")
        print("Error:", result.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate and submit SLURM job scripts for subfolders")
    parser.add_argument("output_path", help="Path to the experiment directory")
    parser.add_argument("-d", "--depth", type=int, default=1, help="Depth of subfolders to process (default: 1)")
    parser.add_argument("-q", "--quest", help="If True, Quest is used", type=bool, default=True)
    parser.add_argument("-a", "--account", help="Account for computational resources")
    parser.add_argument("-p", "--partition", help="Partition for job scheduling")
    parser.add_argument("-n", "--nodes", help="Number of nodes to use", type=int, default=1)
    parser.add_argument("-nt", "--ntasks", help="Number of tasks to run", type=int)
    parser.add_argument("-t", "--job_time", help="Time limit for the job")
    parser.add_argument("-m", "--memory", help="Memory needed to run the job", default="60GB")
    parser.add_argument("-e", "--email", help="Email for job status")
    parser.add_argument("-mt", "--mail_type", help="Mail type for job status")
    
    args = parser.parse_args()
    
    # os.walk() is a generator that yields a 3-tuple for each directory it visits:
    # root: The current directory path
    # dirs: A list of subdirectories in the current directory
    # files: A list of files in the current directory
    for root, dirs, files in os.walk(args.output_path):
        # This line calculates the depth of the current directory relative to the starting path.
        depth = root[len(args.output_path):].count(os.sep)
        # If the current depth match the desired depth 
        # ensures we ony process directories at the specified depth level
        if depth == args.depth:
            subfolder = os.path.relpath(root, args.output_path)
            create_and_submit_slurm_script(args.output_path, subfolder, args)

