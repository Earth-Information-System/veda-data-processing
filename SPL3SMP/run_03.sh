#!/usr/bin/env bash
#SBATCH --account=s3673
#SBATCH --time=03:00:00
#SBATCH --mem=128G
#SBATCH --constraint='sky|cas|hasw'
#SBATCH --output=logs/rechunk-%j.log

source ~/.bashrc
mod_py39

echo "Python path: $(conda run -n smap-download which python)"
echo "Python version: $(conda run -n smap-download python --version)"
echo "Working directory: $(pwd)"

echo "Starting script"
conda run -n smap-download python -u 03-rechunk.py > logs/03-rechunk.log
RESULT=$?
if [[ $RESULT -eq 0 ]]; then
  echo "Done!"
else
  echo "Error executing Python script"
  exit $RESULT
fi
exit
