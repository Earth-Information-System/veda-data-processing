# MERRA-2-IceSheets
This directory contains scripts that preprocess MERRA-2 variables for firn and ice sheet modeling. The high-level steps for both scripts are:

1. Open the necessary files
1. Subset MERRA-2 variables geographically by latitude/longitude bounding box
1. Resample to daily resolution
1. Save to zarr format

There are two different Python scripts that perform the preprocessing:
1. MERRA_concat_discover.py
1. MERRA_to_Zarr.py

They differ slightly in functionality and MERRA_to_Zarr.py outputs one zarr file for each day, whereas MERRA_concat_discover.py outputs one large zarr file with all daily data.

There is also an sbatch script (MERRA_concat_discover.j) that launches MERRA_concat_discover.py on Discover.

Note that these scripts are not meant to be generic and are very much tailored to the needs of the firn and ice sheet modeling work that is being done as part of the Earth Information System: Sea-Level Change project.
