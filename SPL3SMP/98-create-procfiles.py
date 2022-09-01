import xarray as xr
import pandas as pd
import glob

dat = xr.open_zarr("SPL3SMP.zarr")
dates = pd.to_datetime(dat.datetime.values)

fpattern = "SPL3SMP/SMAP_L3_SM_P_{y:04}{m:02d}{d:02d}{s}.h5"

flist = []
for d in dates:
    # d = dates[0]
    search = fpattern.format(y = d.year, m = d.month, d = d.day, s = "*") 
    found = glob.glob(search)
    assert len(found) == 1
    flist.append(found[0])

with open("processed-files.txt", "w") as f:
    f.writelines(fs + "\n" for fs in flist)
