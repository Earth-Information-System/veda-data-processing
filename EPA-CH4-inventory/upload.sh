#!/usr/bin/env bash
aws s3 sync processed-output/ s3://veda-data-store-staging/EIS/cog/EPA-inventory-2012
