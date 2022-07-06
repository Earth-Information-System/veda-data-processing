#!/usr/bin/env bash

env | grep AWS &> /dev/null || (echo 'Error: No AWS credentials detected in environment.' && exit 1)

aws s3 sync processed-output/by-date-version/ s3://veda-data-store-staging/EIS/geoparquet/FIRMS-active-fires/
