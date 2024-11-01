import os
import sys
import logging
import time
import pandas as pd

# Import functions from pipeline_functions.py
from pipeline_functions import (
    check_and_download_reference,
    download_samples,
    convert_sra_to_fq,
    run_fastqc,
    trim_sequences,
    check_and_create_bwa_index,
    align_and_sort_reads,
    mark_duplicates,
    index_bam,
    call_variants
)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def main():
    # Define the number of threads
    num_threads = 4  # Adjust as needed based on your system capabilities

    # Define the base directory (current working directory)
    base_dir = os.getcwd()

    # Define directories for the workflow
    directories = {
        'raw_sequences': os.path.join(base_dir, 'raw_sequences'),
        'qc_sequences': os.path.join(base_dir, 'qc_sequences'),
        'quality_reports': os.path.join(base_dir, 'quality_reports'),
        'alignment': os.path.join(base_dir, 'alignment'),
        'marked_duplicates': os.path.join(base_dir, 'marked_duplicates'),
        'vcf_files': os.path.join(base_dir, 'vcf_files'),
        'plink_files': os.path.join(base_dir, 'plink_files'),
        'kinship_files': os.path.join(base_dir, 'kinship_files')
    }

    # Create necessary directories if they don't exist
    for dir_name, dir_path in directories.items():
        os.makedirs(dir_path, exist_ok=True)
        logging.info(f"Directory ready: {dir_path}")

    # Define pre and post-trim directories for FastQC reports
    quality_reports_dir = directories['quality_reports']
    pre_trim_dir = os.path.join(quality_reports_dir, 'pre_trim')
    post_trim_dir = os.path.join(quality_reports_dir, 'post_trim')

    for dir_path in [pre_trim_dir, post_trim_dir]:
        os.makedirs(dir_path, exist_ok=True)
        logging.info(f"Directory ready: {dir_path}")

    # Define the reference genome paths
    reference_genome_dir = os.path.join(base_dir, 'reference_genome')
    reference_genome_file = os.path.join(reference_genome_dir, 'Homo_sapiens.GRCh38.dna.primary_assembly.fa')

    # Create reference genome directory
    os.makedirs(reference_genome_dir, exist_ok=True)

    # Define sample information directory and file path
    sample_info_dir = os.path.join(base_dir, 'sample_information')
    srr_file_path = os.path.join(sample_info_dir, 'SRR_Acc_List.txt')

    # Create sample information directory
    os.makedirs(sample_info_dir, exist_ok=True)

    # Check if the sample_information directory exists
    if not os.path.isdir(sample_info_dir):
        logging.error(f"Directory not found: {sample_info_dir}. No file list available.")
        sys.exit(1)
    else:
        # Check if the SRR_Acc_List.txt file exists
        if not os.path.isfile(srr_file_path):
            logging.error(f"File not found: {srr_file_path}. No file list available.")
            sys.exit(1)
        else:
            logging.info(f"File list available: {srr_file_path}")

    # Load sample list
    sample_file = srr_file_path
    try:
        samples = pd.read_csv(sample_file, header=None)[0].tolist()
        logging.info(f"Loaded samples: {samples}")
    except Exception as e:
        logging.error(f"Failed to load samples from {sample_file}: {e}")
        sys.exit(1)

    # Define path to Trimmomatic JAR file
    trimmomatic_jar_path = os.path.join(base_dir, 'trimmomatic.jar')
    if not os.path.isfile(trimmomatic_jar_path):
        logging.error(f"Trimmomatic JAR file not found at {trimmomatic_jar_path}. Please ensure it is available.")
        sys.exit(1)

    # Define path to Picard (assumed to be available in PATH)
    picard_cmd = 'picard'  # If Picard is not in PATH, provide the full path

    # Check and download reference genome
    if not check_and_download_reference(reference_genome_dir, reference_genome_file):
        logging.error("Failed to prepare reference genome. Exiting.")
        sys.exit(1)

    # Download samples from SRA
    download_samples(samples)

    # Convert SRA files to FASTQ format
    convert_sra_to_fq(samples, directories)

    # Run FastQC before trimming
    for sample in samples:
        input_files = [
            os.path.join(directories['raw_sequences'], f"{sample}_1.fastq"),
            os.path.join(directories['raw_sequences'], f"{sample}_2.fastq")
        ]
        run_fastqc(input_files, output_dir=pre_trim_dir, num_threads=num_threads)


if __name__ == "__main__":
    main()
