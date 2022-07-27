#!/usr/bin/env bash

for f in $(find FWI_GPM_LATE_v5_Daily.zarr/ -type f -user rfield1); do
  newname="$f-rf"
  mv "$f" "$newname"
  cp "$newname" "$f"
  rm -f "$newname"
done
