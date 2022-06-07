#!/usr/bin/env bash

aws s3 sync SPL3SMP.zarr/ s3://veda-data-store-staging/EIS/zarr/SPL3SMP.zarr/
aws s3 sync SPL3SMP-cog/ s3://veda-data-store-staging/EIS/cog/SPL3SMP/
