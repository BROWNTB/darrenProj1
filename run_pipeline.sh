#!/bin/bash

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "Conda is not installed or not in the PATH. Please install Conda before running this script."
    exit 1
fi

# Usage: ./run_pipeline.sh [environment_name]
# Default environment name: bioinfo_env

# Set default environment name if not provided
ENV_NAME=${1:-bioinfo_env}

# Specify the desired Python version
PYTHON_VERSION=3.9

# Function to check if a Conda environment exists
function conda_env_exists() {
    conda info --envs | grep -w "$1" > /dev/null 2>&1
}

# Function to create a Conda environment
function create_conda_env() {
    echo "Creating Conda environment '$1' with Python $PYTHON_VERSION..."
    conda create -y -n "$1" python="$PYTHON_VERSION"
}

# Start of the script
echo "Checking for Conda environment: $ENV_NAME"

if conda_env_exists "$ENV_NAME"; then
    echo "Environment '$ENV_NAME' exists."
else
    echo "Environment '$ENV_NAME' does not exist."
    create_conda_env "$ENV_NAME"
fi

# Activate the Conda environment
echo "Activating environment '$ENV_NAME'"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Run the installation script
echo "Running install_dependencies.py..."
python install_dependencies.py
if [ $? -ne 0 ]; then
    echo "Installation failed. Exiting."
    exit 1
fi

# Run the pipeline script
echo "Running pipeline_script.py..."
python pipeline_script.py
if [ $? -ne 0 ]; then
    echo "Pipeline execution failed. Exiting."
    exit 1
fi

echo "Pipeline completed successfully."
