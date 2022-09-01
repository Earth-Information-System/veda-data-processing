#!/usr/bin/env bash
#SBATCH --account=s3673
#SBATCH --time=03:00:00
#SBATCH --mem=128G
#SBATCH --constraint='sky|cas|hasw'
#SBATCH --output=logs/zarr2cog-%j.log

source ~/.bashrc
mod_py39

echo "Starting script"
conda run -n vdp-common python -u ../common/zarr2cog.py \
  --outdir=SPL3SMP-cog \
  --dataset=SPL3SMP.zarr \
  --prefix=SPL3SMP \
  --crs=EPSG:6933 \
  --x_dim=easting_m \
  --y_dim=northing_m \
  --timevar=datetime \
  --timeunit=D \
  > logs/zarr2cog.log
