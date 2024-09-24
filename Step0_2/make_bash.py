
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

def create_and_submit_slurm_script(output_path, subfolder, args):
    today = date.today().strftime('%Y-%m-%d')
    current_time = datetime.now().strftime("%H.%M")
    
    job_directory = os.path.join(output_path, subfolder, "jobs")
    logs_directory = os.path.join(output_path, subfolder, "logs")
    os.makedirs(job_directory, exist_ok=True)
    os.makedirs(logs_directory, exist_ok=True)
    
    job_file = os.path.join(job_directory, f"job_{subfolder}_{today}_{current_time}.sh")
    
    with open(job_file, 'w') as fh:
        fh.write("#!/bin/bash\n")
        if args.quest:
            fh.write(f"#SBATCH --account={args.account}\n")
            fh.write(f"#SBATCH --nodes={args.nodes}\n")
        fh.write(f"#SBATCH --partition={args.partition}\n")
        fh.write(f"#SBATCH --ntasks={args.ntasks}\n")
        fh.write(f"#SBATCH --time={args.time}\n")
        fh.write(f"#SBATCH --mem={args.memory}\n")
        fh.write(f"#SBATCH --mail-user={args.email}\n")
        fh.write(f"#SBATCH --mail-type={args.mail_type}\n")
        fh.write(f"#SBATCH --job-name={subfolder}_gDNApipeline\n")
        fh.write(f"#SBATCH --output={logs_directory}/%x_%j.out\n")
        fh.write(f"#SBATCH --error={logs_directory}/%x_%j.err\n\n")
        
        fh.write("module purge\n")
        fh.write('eval "$(conda shell.bash hook)"\n\n')
        fh.write("conda activate Barcode_extraction\n\n")
        
        # Fixed last line as specified
        path = os.path.join(output_path,subfolder,"templateInputFullRun.py")
        fh.write(f"python3 {path}\n")
    
    print(f"Created job script for {subfolder}: {job_file}")

    # Submit the job
    result = subprocess.run(["sbatch", job_file], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"Job for {subfolder} submitted successfully")
        print("Output:", result.stdout)
    else:
        print(f"Job submission for {subfolder} failed")
        print("Error:", result.stderr)

def main():
    parser = argparse.ArgumentParser(description="Generate and submit SLURM job scripts for subfolders")
    parser.add_argument("output_path", help="Path to the experiment directory")
    parser.add_argument("-d", "--depth", type=int, default=1, help="Depth of subfolders to process (default: 1)")
    parser.add_argument("-q", "--quest", help="If True, Quest is used", type=bool, default=True)
    parser.add_argument("-a", "--account", help="Account for computational resources")
    parser.add_argument("-p", "--partition", help="Partition for job scheduling")
    parser.add_argument("-n", "--nodes", help="Number of nodes to use", type=int, default=1)
    parser.add_argument("-nt", "--ntasks", help="Number of tasks to run", type=int)
    parser.add_argument("-t", "--time", help="Time limit for the job")
    parser.add_argument("-m", "--memory", help="Memory needed to run the job", default="60GB")
    parser.add_argument("-e", "--email", help="Email for job status")
    parser.add_argument("-mt", "--mail_type", help="Mail type for job status")
    
    args = parser.parse_args()
    
    for root, dirs, files in os.walk(args.output_path):
        depth = root[len(args.output_path):].count(os.sep)
        if depth == args.depth:
            subfolder = os.path.relpath(root, args.output_path)
            create_and_submit_slurm_script(args.output_path, subfolder, args)

if __name__ == "__main__":
    main()
