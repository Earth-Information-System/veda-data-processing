#!/usr/bin/env bash
#SBATCH --account=s3673
#SBATCH --time=11:59:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=64G
#SBATCH --constraint='sky|cas|hasw'
#SBATCH --output=logs/slurm-%J.log

source ~/.bash_functions
mod_python

echo "Python: $(which python)"
echo "Working directory: $(pwd)"

echo "Starting script"
python -u 02-create-zarr.py
echo "Done!"

exit
