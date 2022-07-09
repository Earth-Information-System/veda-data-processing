#!/usr/bin/env bash
aws s3 sync OCO2_GEOS_L3CO2_day.zarr/ s3://veda-data-store-staging/EIS/zarr/OCO2_GEOS_L3CO2_day.zarr/
aws s3 sync OCO2_GEOS_L3CO2_day.COG/ s3://veda-data-store-staging/EIS/COG/OCO2_GEOS_L3CO2_day/
