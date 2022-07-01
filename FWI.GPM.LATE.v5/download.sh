#!/usr/bin/env bash

wget --recursive --no-parent \
  -nH --cut-dirs=6 \
  -e robots=off \
  -R "index.html*" \
  -P "raw-input/" \
  -nc \
  https://portal.nccs.nasa.gov/datashare/GlobalFWI/v2.0/fwiCalcs.GEOS-5/Default/GPM.LATE.v5/2022/
