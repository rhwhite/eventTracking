# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os, errno
import numpy as np
import numpy.ma as ma
import netCDF4
from netCDF4 import Dataset
import datetime as dt
import re
import sys
import Ngl
import xray
import math

Data = "TRMMrg"

sumlats = 5
sumlons = 5

startyr = 1998 # Don't change - tied to file names!
endyr = 2014

anstartyr = 1998 #year for analysis start
anendyr = 2014 #year for analysis end
nyears = anendyr - anstartyr + 1
nmonths = 12

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

filetimespan = '3hrly'

if Data == "TRMM":
	DirP = '/home/disk/eos4/rachel/Obs/TRMM/'
	FileP = 'SeasAnn_TRMM_' + str(startyr) + '-' + str(endyr) + '_3B42_3hrly_nonan.nc'
elif Data == "TRMMrg":
        DirP = '/home/disk/eos4/rachel/Obs/TRMM/'
        FileP = 'SeasAnn_regrid2ERAI_TRMM_3B42_' + str(startyr) + '-' + str(endyr) + '.nc'

else:
	print "Wrong data type"
	exit()

DirO = DirP 
FileO = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_' + FileP
print DirP + FileP
#pcp, time, ongitude, latitude
precipdata = xray.open_dataset(DirP + FileP)

lons = precipdata['Longitude']
lats = precipdata['Latitude']
precipSeas = precipdata['PrecipSeas']
precipAnn = precipdata['PrecipAnn']
precipSeasC = precipdata['PrecipClimSeas']
precipAnnC = precipdata['PrecipClimAnn']

nseas = precipSeasC.shape[0]
nyears = precipAnn.shape[0]

Seas = precipdata['Season']

nlons = len(lons)
nlats = len(lats)

nlonsnew = nlons/sumlons
nlatsnew = nlats/sumlats

Latsnew = np.zeros(nlatsnew,np.float)
Lonsnew = np.zeros(nlonsnew,np.float)

# Deal with Seasonal and Annual separately
NewS = np.zeros((nyears,nseas,nlatsnew,nlonsnew),np.float)
NewA = np.zeros((nyears,nlatsnew,nlonsnew),np.float)

NewSC = np.zeros((nseas,nlatsnew,nlonsnew),np.float)
NewAC = np.zeros((nlatsnew,nlonsnew),np.float)

#'regrid' by averaging over large boxes (averaging lons and lats)
inlat = 0
for ilats in range(0,nlats,sumlats):
	inlon = 0
	Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlats])
#	print np.mean(lats[ilats:ilats+sumlats])
	for ilons in range(0,nlons,sumlons):
		Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlons])
#		print np.mean(lons[ilons:ilons+sumlons])
		NewS[:,:,inlat,inlon] = np.mean(precipSeas[:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(2,3),dtype=np.float)
                NewA[:,inlat,inlon] = np.mean(precipAnn[:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(1,2),dtype=np.float)

                NewSC[:,inlat,inlon] = np.mean(precipSeasC[:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(1,2),dtype=np.float)
                NewAC[inlat,inlon] = np.mean(precipAnnC[ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(0,1),dtype=np.float)

		inlon += 1
	inlat += 1

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('seas', 4)
ncfile.createDimension('year',nyears)
ncfile.createDimension('lon', nlonsnew)
ncfile.createDimension('lat', nlatsnew)

ODensityMapC = ncfile.createVariable('PrecipAnnClim','f4',('lat','lon'),fill_value=-9999)
ODensityMapSeasC = ncfile.createVariable('PrecipSeasClim','f4',('seas','lat','lon'),fill_value=-9999)
ODensityMap = ncfile.createVariable('PrecipAnn','f4',('year','lat','lon'),fill_value=-9999)
ODensityMapSeas = ncfile.createVariable('PrecipSeas','f4',('year','seas','lat','lon'),fill_value=-9999)


setattr(ODensityMap,'Extra Info','average precip in mm/hr (regridding not weighted by surface area)')
setattr(ODensityMapSeas,'Extra Info','average precip in mm/hr (regridding not weighted by surface area)')
setattr(ODensityMap,'units','mm/hr')
setattr(ODensityMapSeas,'units','mm/hr')

setattr(ODensityMapC,'Extra Info','average precip in mm/hr (regridding not weighted by surface area)')
setattr(ODensityMapSeasC,'Extra Info','average precip in mm/hr (regridding not weighted by surface area)')
setattr(ODensityMapC,'units','mm/hr')
setattr(ODensityMapSeasC,'units','mm/hr')


OLongitude = ncfile.createVariable('lon','f4',('lon'),fill_value=-9999)
OLatitude = ncfile.createVariable('lat','f4',('lat'),fill_value=-9999)
OSeas = ncfile.createVariable('seas',str,('seas'))

OLatitude[:] = Latsnew
OLongitude[:] = Lonsnew
ODensityMapC[...] = NewAC
ODensityMapSeasC[...] = NewSC

print NewA
ODensityMap[...] = NewA
ODensityMapSeas[...] = NewS

OSeas[...] = Seas.values

Ngl.end()









