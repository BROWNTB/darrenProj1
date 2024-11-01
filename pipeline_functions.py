import os
import subprocess
import time
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def check_and_download_reference(reference_genome_dir, reference_genome_file):
    """
    Check if the reference genome and its index exist.
    If not, download and prepare them.

    Args:
        reference_genome_dir (str): Directory where the reference genome should be stored.
        reference_genome_file (str): Path to the reference genome FASTA file.
    """
    # Create the reference genome directory if it doesn't exist
    if not os.path.exists(reference_genome_dir):
        os.makedirs(reference_genome_dir)
        logging.info(f"Created directory: {reference_genome_dir}")

    # Define the index file path
    index_file = reference_genome_file + '.fai'

    # Check if the reference genome file exists
    if not os.path.isfile(reference_genome_file):
        logging.info(f"{reference_genome_file} not found. Downloading...")
        # Download the reference genome
        ref_url = 'ftp://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
        gz_file = reference_genome_file + '.gz'
        try:
            subprocess.run(['wget', '-O', gz_file, ref_url], check=True)
            logging.info("Download completed.")
            # Unzip the downloaded file
            subprocess.run(['gunzip', '-f', gz_file], check=True)
            logging.info("Unzipping completed.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to download or unzip the reference genome: {e}")
            return False
    else:
        logging.info(f"Reference genome file already exists: {reference_genome_file}")

    # Check if the index file exists
    if not os.path.isfile(index_file):
        logging.info(f"{index_file} not found. Creating index...")
        # Index the reference genome
        try:
            subprocess.run(['samtools', 'faidx', reference_genome_file], check=True)
            logging.info("Indexing completed.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to index the reference genome: {e}")
            return False
    else:
        logging.info(f"Index file already exists: {index_file}")
    return True

def download_samples(samples):
    """
    Download samples from the SRA using the prefetch command.

    Args:
        samples (list): A list of sample identifiers.
    """
    for sample in samples:
        start_time = time.time()  # Start the timer
        logging.info(f"Downloading sample: {sample}")
        try:
            subprocess.run(['prefetch', sample], check=True)
            elapsed_time = time.time() - start_time  # Calculate elapsed time
            logging.info(f"Finished downloading sample: {sample} in {elapsed_time / 60:.2f} minutes.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to download sample {sample}: {e}")
            continue

def validate_file_exists(file_path):
    """
    Validate that a file exists at the specified path.

    Args:
        file_path (str): The path to the file to check.

    Returns:
        bool: True if the file exists, False otherwise.
    """
    if not os.path.isfile(file_path):
        logging.error(f"Expected output file not found: {file_path}")
        return False
    return True

def convert_sra_to_fq(samples, directories):
    """
    Convert downloaded .sra files to .fq (FASTQ) format.

    Args:
        samples (list): A list of sample identifiers.
        directories (dict): A dictionary containing paths to directories.
    """
    for sample in samples:
        start_time = time.time()  # Start the timer
        logging.info(f"Converting {sample}.sra to FASTQ files...")
        sra_file = os.path.join(sample, f"{sample}.sra")
        output_dir = directories['raw_sequences']
        cmd_fastq_dump = [
            'fastq-dump',
            '--split-files',
            '--outdir',
            output_dir,
            sra_file
        ]
        try:
            subprocess.run(cmd_fastq_dump, check=True)
            elapsed_time = time.time() - start_time  # Calculate elapsed time
            logging.info(f"Finished converting {sample}.sra to FASTQ files in {elapsed_time / 60:.2f} minutes.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to convert {sample}.sra to FASTQ: {e}")
            continue

def run_fastqc(input_files, output_dir, num_threads=2):
    """
    Run FastQC on the specified FASTQ files to perform quality control.

    Args:
        input_files (list): A list of paths to the FASTQ files to be analyzed.
        output_dir (str): Directory to save FastQC reports.
        num_threads (int): Number of threads to use.
    """
    logging.info("Running FastQC...")
    for file in input_files:
        start_time = time.time()
        # Execute FastQC
        cmd_fastqc = [
            'fastqc',
            file,
            '-o',
            output_dir,
            '--quiet',
            '--threads',
            str(num_threads)
        ]
        logging.info(f"Executing command: {' '.join(cmd_fastqc)}")
        try:
            subprocess.run(cmd_fastqc, check=True)
            # Check for expected output files
            base_filename = os.path.basename(file).replace('.fastq', '')
            html_report = os.path.join(output_dir, f"{base_filename}_fastqc.html")
            zip_report = os.path.join(output_dir, f"{base_filename}_fastqc.zip")

            if os.path.isfile(html_report) and os.path.isfile(zip_report):
                logging.info(f"FastQC reports generated for {file}")
            else:
                logging.warning(f"Expected FastQC output files not found for {file}.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error occurred during FastQC for {file}: {e}")
            continue
        elapsed_time = time.time() - start_time
        logging.info(f"Finished FastQC for {file} in {elapsed_time / 60:.2f} minutes.")
    logging.info("Finished running FastQC.")

def trim_sequences(sample, trimmomatic_jar_path, directories, num_threads=2):
    """
    Trim sequences using Trimmomatic to remove low-quality reads.

    Args:
        sample (str): The identifier for the sample to be trimmed.
        trimmomatic_jar_path (str): Path to the Trimmomatic JAR file.
        directories (dict): Dictionary of directory paths.
        num_threads (int): Number of threads to use.
    """
    start_time = time.time()  # Start the timer
    # Define paths for trimming
    fq1 = os.path.join(directories['raw_sequences'], f"{sample}_1.fastq")
    fq2 = os.path.join(directories['raw_sequences'], f"{sample}_2.fastq")
    trimmed_fq1 = os.path.join(directories['qc_sequences'], f"{sample}_1_paired.fastq")
    unpaired_fq1 = os.path.join(directories['qc_sequences'], f"{sample}_1_unpaired.fastq")
    trimmed_fq2 = os.path.join(directories['qc_sequences'], f"{sample}_2_paired.fastq")
    unpaired_fq2 = os.path.join(directories['qc_sequences'], f"{sample}_2_unpaired.fastq")
    log_file = os.path.join(directories['quality_reports'], f"{sample}_trim_report.log")

    # Trimming command
    cmd_trim = [
        'java',
        '-jar',
        trimmomatic_jar_path,
        'PE',
        '-threads',
        str(num_threads),
        '-phred33',
        fq1,
        fq2,
        trimmed_fq1,
        unpaired_fq1,
        trimmed_fq2,
        unpaired_fq2,
        'LEADING:3',
        'TRAILING:3',
        'SLIDINGWINDOW:4:20',
        'MINLEN:36'
    ]

    logging.info(f"Trimming sequences for sample {sample}")
    try:
        with open(log_file, 'w') as logfile:
            subprocess.run(cmd_trim, stdout=logfile, stderr=subprocess.STDOUT, check=True)
        elapsed_time = time.time() - start_time
        logging.info(f"Finished trimming for sample {sample} in {elapsed_time / 60:.2f} minutes.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred during trimming for sample {sample}: {e}")

def check_and_create_bwa_index(reference_genome_file):
    """
    Check for the presence of BWA index files for the reference genome.

    If the index files are not present, this function will create them using
    the BWA index command.

    Args:
        reference_genome_file (str): Path to the reference genome FASTA file.
    """
    index_files = [
        f"{reference_genome_file}.amb",
        f"{reference_genome_file}.ann",
        f"{reference_genome_file}.bwt",
        f"{reference_genome_file}.pac",
        f"{reference_genome_file}.sa"
    ]

    # Check if all index files exist
    if all(os.path.isfile(index_file) for index_file in index_files):
        logging.info("BWA index files are present.")
    else:
        logging.info("BWA index files not found. Creating BWA index...")
        cmd_index = ['bwa', 'index', reference_genome_file]
        try:
            subprocess.run(cmd_index, check=True)
            logging.info("BWA index created.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error occurred while creating BWA index: {e}")
            return False
    return True


import subprocess
import logging
import os
import time

def align_and_sort_reads(sample, reference_genome_file, directories, num_threads=4):
    """
    Aligns trimmed FASTQ files to the reference genome using BWA and sorts the output BAM file.
    Fixes the header to ensure only one @HD line is present.

    Args:
        sample (str): The sample identifier.
        reference_genome_file (str): Path to the reference genome FASTA file.
        directories (dict): Dictionary of directory paths.
        num_threads (int): Number of threads to use.
    """
    try:
        logging.info(f"Aligning and sorting reads for sample {sample}")
        start_time = time.time()

        # Input files
        fq1 = os.path.join(directories['qc_sequences'], f"{sample}_1_paired.fastq")
        fq2 = os.path.join(directories['qc_sequences'], f"{sample}_2_paired.fastq")

        # Output file
        sorted_bam_output = os.path.join(directories['alignment'], f"{sample}_sorted.bam")

        # Read group information
        read_group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"

        # BWA command with read group
        bwa_cmd = [
            'bwa', 'mem',
            '-t', str(num_threads),
            '-R', read_group,
            reference_genome_file,
            fq1, fq2
        ]

        # Samtools sort command
        samtools_sort_cmd = [
            'samtools', 'sort',
            '-@', str(num_threads),
            '-o', sorted_bam_output,
            '-'
        ]

        # Execute BWA and pipe into Samtools for sorting
        with subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE) as bwa_proc:
            with subprocess.Popen(samtools_sort_cmd, stdin=bwa_proc.stdout, stdout=subprocess.PIPE) as sort_proc:
                bwa_proc.stdout.close()  # Close BWA stdout when piping
                sort_proc.communicate()  # Execute samtools sort

        # Check if output file exists and is not empty
        if os.path.isfile(sorted_bam_output) and os.path.getsize(sorted_bam_output) > 0:
            elapsed_time = time.time() - start_time
            logging.info(f"Finished alignment and sorting for sample {sample} in {elapsed_time / 60:.2f} minutes.")
        else:
            logging.error(f"Sorted BAM file not created for sample {sample}")

    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred during alignment and sorting for sample {sample}: {e}")
    except Exception as e:
        logging.error(f"Unexpected error during alignment and sorting for sample {sample}: {e}")

        # Run BWA and pipe to samtools sort
        with subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE) as bwa_process:
            with subprocess.Popen(samtools_sort_cmd, stdin=bwa_process.stdout) as sort_process:
                bwa_process.stdout.close()  # Allow bwa_process to receive a SIGPIPE if sort_process exits
                sort_process.communicate()

        # After sorting, fix the header to remove duplicate @HD lines
        fix_bam_header(sorted_bam_output)

        end_time = time.time()
        time_elapsed = (end_time - start_time) / 60
        logging.info(f"Finished alignment and sorting for sample {sample} in {time_elapsed:.2f} minutes.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred during alignment and sorting for sample {sample}: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred during alignment and sorting for sample {sample}: {e}")

def fix_bam_header(bam_file):
    """
    Fixes the BAM file header to ensure only one @HD line is present.

    Args:
        bam_file (str): Path to the BAM file to be fixed.
    """
    try:
        # Extract the header
        header_output = subprocess.check_output(['samtools', 'view', '-H', bam_file], text=True)
        header_lines = header_output.strip().split('\n')

        # Remove duplicate @HD lines
        new_header_lines = []
        hd_found = False
        for line in header_lines:
            if line.startswith('@HD'):
                if not hd_found:
                    # Keep the first @HD line
                    new_header_lines.append(line)
                    hd_found = True
                else:
                    # Skip additional @HD lines
                    continue
            else:
                new_header_lines.append(line)

        # Write the new header to a temporary file
        temp_header_file = bam_file + '.header'
        with open(temp_header_file, 'w') as hf:
            hf.write('\n'.join(new_header_lines) + '\n')

        # Create a temporary BAM file with the fixed header
        temp_bam_file = bam_file + '.tmp'
        with open(temp_bam_file, 'wb') as temp_bam:
            subprocess.run(['samtools', 'reheader', temp_header_file, bam_file], stdout=temp_bam, check=True)

        # Replace the original BAM file with the fixed one
        os.replace(temp_bam_file, bam_file)
        os.remove(temp_header_file)
        logging.info(f"Fixed BAM header for {bam_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred while fixing BAM header for {bam_file}: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred while fixing BAM header for {bam_file}: {e}")
    """
    Aligns trimmed FASTQ files to the reference genome using BWA and sorts the output BAM file.

    Args:
        sample (str): The sample identifier.
        reference_genome_file (str): Path to the reference genome FASTA file.
        directories (dict): Dictionary of directory paths.
        num_threads (int): Number of threads to use.
    """
    try:
        logging.info(f"Aligning and sorting reads for sample {sample}")
        start_time = time.time()

        # Input files
        fq1 = os.path.join(directories['qc_sequences'], f"{sample}_1_paired.fastq")
        fq2 = os.path.join(directories['qc_sequences'], f"{sample}_2_paired.fastq")

        # Output file
        sorted_bam_output = os.path.join(directories['alignment'], f"{sample}_sorted.bam")

        # Read group information
        read_group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"

        # BWA command with read group
        bwa_cmd = [
            'bwa', 'mem',
            '-t', str(num_threads),
            '-R', read_group,
            reference_genome_file,
            fq1, fq2
        ]

        # Samtools sort command
        samtools_sort_cmd = [
            'samtools', 'sort',
            '-@', str(num_threads),
            '-o', sorted_bam_output
        ]

        # Run BWA and pipe to samtools sort
        with subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as bwa_process:
            with subprocess.Popen(
                samtools_sort_cmd,
                stdin=bwa_process.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            ) as sort_process:
                # Allow bwa_process to receive a SIGPIPE if sort_process exits
                bwa_process.stdout.close()

                # Wait for sort_process to complete and capture stderr
                _, stderr = sort_process.communicate()

                # Check for errors in samtools sort
                if sort_process.returncode != 0:
                    logging.error(f"Samtools sort failed for sample {sample}: {stderr.decode('utf-8')}")
                    raise subprocess.CalledProcessError(sort_process.returncode, samtools_sort_cmd)

        end_time = time.time()
        time_elapsed = (end_time - start_time) / 60
        logging.info(f"Finished alignment and sorting for sample {sample} in {time_elapsed:.2f} minutes.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred during alignment and sorting for sample {sample}: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred during alignment and sorting for sample {sample}: {e}")


def mark_duplicates(sample, picard_cmd, directories):
    """
    Marks duplicates in BAM files using Picard.

    Args:
        sample (str): The sample identifier.
        picard_cmd (str): Command or path to run Picard.
        directories (dict): Dictionary of directory paths.
    """
    input_bam = os.path.join(directories['alignment'], f"{sample}_sorted.bam")
    output_bam = os.path.join(directories['marked_duplicates'], f"{sample}_marked.bam")
    metrics_file = os.path.join(directories['quality_reports'], f"{sample}_dup_metrics.txt")

    if os.path.isfile(input_bam) and os.path.getsize(input_bam) > 0:
        cmd_mark_duplicates = [
            picard_cmd,
            'MarkDuplicates',
            f'I={input_bam}',
            f'O={output_bam}',
            f'M={metrics_file}',
            'REMOVE_DUPLICATES=true'
        ]
        logging.info(f"Marking duplicates for sample {sample}")
        try:
            subprocess.run(cmd_mark_duplicates, check=True)
            logging.info(f"Finished marking duplicates for sample {sample}")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error occurred during marking duplicates for sample {sample}: {e}")
    else:
        logging.warning(f"Input BAM file not found or empty for sample {sample}. Skipping marking duplicates.")

def index_bam(sample, directories):
    """
    Indexes BAM files using Samtools.

    Args:
        sample (str): The sample identifier.
        directories (dict): Dictionary of directory paths.
    """
    marked_bam_file = os.path.join(directories['marked_duplicates'], f"{sample}_marked.bam")

    start_time = time.time()
    cmd_index = [
        'samtools', 'index',
        marked_bam_file
    ]
    logging.info(f"Indexing BAM file for sample {sample}")
    try:
        subprocess.run(cmd_index, check=True)
        elapsed_time = time.time() - start_time
        logging.info(f"Finished indexing BAM for sample {sample} in {elapsed_time / 60:.2f} minutes.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred during BAM indexing for sample {sample}: {e}")

def call_variants(sample, reference_genome_file, directories):
    """
    Call variants from BAM files using BCFtools and create a VCF file.

    Args:
        sample (str): The sample identifier.
        reference_genome_file (str): Path to the reference genome FASTA file.
        directories (dict): Dictionary of directory paths.
    """
    bam_file = os.path.join(directories['marked_duplicates'], f"{sample}_marked.bam")
    vcf_file = os.path.join(directories['vcf_files'], f"{sample}.vcf")

    start_time = time.time()
    cmd_bcftools_mpileup = [
        'bcftools', 'mpileup',
        '-Ou',
        '-f', reference_genome_file,
        bam_file
    ]
    cmd_bcftools_call = [
        'bcftools', 'call',
        '-mv',
        '-o', vcf_file
    ]
    logging.info(f"Calling variants for sample {sample}")
    try:
        mpileup = subprocess.Popen(cmd_bcftools_mpileup, stdout=subprocess.PIPE)
        subprocess.run(cmd_bcftools_call, stdin=mpileup.stdout, check=True)
        mpileup.stdout.close()
        elapsed_time = time.time() - start_time
        logging.info(f"Finished variant calling for sample {sample} in {elapsed_time / 60:.2f} minutes.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error occurred during variant calling for sample {sample}: {e}")
