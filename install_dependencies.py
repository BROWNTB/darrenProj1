import subprocess
import sys
import logging
import time
import os
import shutil

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def configure_conda_channels():
    """Configure Conda channels for package installation."""
    channels = ['defaults', 'bioconda', 'conda-forge']
    for channel in channels:
        subprocess.run(['conda', 'config', '--add', 'channels', channel], check=True)


# Function to install a package using Conda
def install_package(package_name):
    """Install a package using Conda if it's not already installed."""
    logging.info(f"Checking if {package_name} is installed.")
    try:
        # Check if the package is installed
        result = subprocess.run(['conda', 'list', package_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = result.stdout.decode()
        if package_name.lower() not in output.lower():
            logging.info(f"{package_name} is not installed. Installing via Conda...")
            start_time = time.time()  # Start timer for this installation
            subprocess.run(['conda', 'install', '-y', package_name], check=True)
            elapsed_time = time.time() - start_time  # Calculate elapsed time
            logging.info(f"{package_name} installed successfully in {elapsed_time / 60:.2f} minutes.")
        else:
            logging.info(f"{package_name} is already installed.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to install {package_name} via Conda: {e}")
        sys.exit(1)

# Function to download files with wget (using conda-installed wget)
def download_file(url, dest_path):
    """Download a file from a URL to a destination path."""
    if not os.path.exists(dest_path):
        logging.info(f"Downloading {url} to {dest_path}")
        start_time = time.time()  # Start timer for this download
        try:
            subprocess.run(['wget', '-O', dest_path, url], check=True)
            elapsed_time = time.time() - start_time  # Calculate elapsed time
            logging.info(f"Downloaded {dest_path} in {elapsed_time / 60:.2f} minutes.")
        except subprocess.CalledProcessError as e:
            logging.error(f"Failed to download {url}: {e}")
            sys.exit(1)
    else:
        logging.info(f"File {dest_path} already exists.")

# Function to unzip files
def unzip_file(zip_path, extract_to='.'):
    """Unzip a file to the specified directory."""
    logging.info(f"Unzipping {zip_path} to {extract_to}")
    start_time = time.time()  # Start timer for unzipping
    try:
        subprocess.run(['unzip', '-o', zip_path, '-d', extract_to], check=True)
        elapsed_time = time.time() - start_time  # Calculate elapsed time
        logging.info(f"Unzipped {zip_path} in {elapsed_time:.2f} seconds.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to unzip {zip_path}: {e}")
        sys.exit(1)

# Function to check for and download Trimmomatic JAR file
def check_trimmomatic():
    """Ensure that Trimmomatic JAR file is available."""
    trimmomatic_jar_path = os.path.join(os.getcwd(), 'trimmomatic.jar')
    if not os.path.exists(trimmomatic_jar_path):
        logging.info("Trimmomatic JAR file not found. Downloading...")
        trimmomatic_url = 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip'
        trimmomatic_zip = os.path.join(os.getcwd(), 'Trimmomatic-0.39.zip')
        download_file(trimmomatic_url, trimmomatic_zip)
        unzip_file(trimmomatic_zip, os.getcwd())
        # Move the JAR file to the current directory
        shutil.move(os.path.join(os.getcwd(), 'Trimmomatic-0.39', 'trimmomatic-0.39.jar'), trimmomatic_jar_path)
        # Clean up
        os.remove(trimmomatic_zip)
        shutil.rmtree(os.path.join(os.getcwd(), 'Trimmomatic-0.39'))
        logging.info(f"Trimmomatic JAR file is ready at {trimmomatic_jar_path}")
    else:
        logging.info(f"Trimmomatic JAR file already exists at {trimmomatic_jar_path}")
    return trimmomatic_jar_path

def main():
    # Start the overall timer
    total_start_time = time.time()

    # Configure Conda channels
    configure_conda_channels()

    # List of required packages
    required_packages = [
        'samtools',
        'sra-tools',   # Use 'sra-tools' instead of 'sratoolkit'
        'fastqc',
        'bwa',
        'vcftools',
        'bcftools',
        'wget',
        'unzip',
        'picard',
        'openjdk',     # For Java dependencies
        'plink',
        'pandas'
    ]

    # Install required packages
    for package in required_packages:
        install_package(package)

    # Check and prepare Trimmomatic
    trimmomatic_jar_path = check_trimmomatic()

    # End of the overall timer
    total_end_time = time.time()
    total_time_elapsed = total_end_time - total_start_time

    # Convert total time into minutes and seconds
    minutes = int(total_time_elapsed // 60)
    seconds = total_time_elapsed % 60

    logging.info(f"All installations and setups are complete in {minutes} minutes and {seconds:.2f} seconds.")

if __name__ == "__main__":
    main()
