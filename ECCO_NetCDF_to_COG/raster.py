#!/bin/env python

from osgeo import gdal, osr
from osgeo import gdalconst
import numpy

def writeArrayAsRasterBand(filename,geoTransform,array,noDataValue,metadataDict=None,wktProj=None,epsg=None,dataType='Float32'):
   cols = array.shape[1]
   rows = array.shape[0]

   dt = eval('gdal.GDT_' + dataType)
   
   driver = gdal.GetDriverByName('GTiff')
   outRaster = driver.Create(filename, cols, rows, 1, dt)
   outRaster.SetGeoTransform(geoTransform)
   if metadataDict is not None:
      outRaster.SetMetadata( metadataDict )
   outBand = outRaster.GetRasterBand(1)

   arrayout = numpy.where(~numpy.isnan(array), array, noDataValue)

   outBand.WriteArray(arrayout)
   outBand.SetNoDataValue(noDataValue)
   outRasterSRS = osr.SpatialReference()
   if wktProj:
       outRasterSRS.ImportFromWkt(wktProj)
   elif epsg:
       outRasterSRS.ImportFromEPSG(epsg)
   else:
       outRasterSRS.ImportFromEPSG(3413)
   outRaster.SetProjection(outRasterSRS.ExportToWkt())
   outBand.FlushCache()

