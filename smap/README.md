# Soil Moisture Active/Passive (SMAP) data processing

SMAP data originate from the National Snow and Ice Data Center (NSIDC).
This specific workflow is for [SMAP L3 Radiometer Global Daily 36 km EASE-Grid Soil moisture][smap].

- `01-nsidc-download.py` -- Download script for SMAP data, for a given temporal subset (set near the top of the file). 99% of this script originates from the NSIDC data download Python script. Here, it has been adapted specifically for SMAP.

- `02-create-zarr.py` -- Convert downloaded SMAP data to Zarr. A key component of this script is to add EASE grid coordinates, based on the EASE grid transformation specified in the SMAP algorithm document.

[smap]: https://nsidc.org/data/SPL3SMP
