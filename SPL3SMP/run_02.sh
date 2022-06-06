#!/usr/bin/env bash
#SBATCH --account=s3673
#SBATCH --time=11:59:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=128G
#SBATCH --constraint='sky|cas|hasw'
#SBATCH --output=logs/create-zarr-%j.log

source ~/.bash_functions
mod_py39

echo "Conda path: $(which conda)"
echo "Python path: $(conda run -n smap-download which python)"
echo "Python version: $(conda run -n smap-download python --version)"
echo "Working directory: $(pwd)"

echo "Starting script"
conda run -n smap-download python -u 02-create-zarr.py > logs/02-create-zarr.log
RESULT=$?
if [[ $RESULT -eq 0 ]]; then
  echo "Done!"
else
  echo "Error executing Python script"
  exit $RESULT
fi
