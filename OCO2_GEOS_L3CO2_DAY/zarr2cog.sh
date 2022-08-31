#!/usr/bin/env bash

mkdir -p logs

conda run -n vdp-common python -u ../common/zarr2cog.py \
  --dataset="OCO2_GEOS_L3CO2_day.zarr" \
  --outdir="OCO2_GEOS_L3CO2_day.COG" \
  --prefix="OCO2_GEOS_L3CO2_day" \
  --crs="EPSG:4326" \
  --x_dim="lon" \
  --y_dim="lat" \
  --timevar="time" \
  --timeunit="D" \
  > logs/zarr2cog.log
