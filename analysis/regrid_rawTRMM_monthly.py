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


sumlats = 400
sumlons = 1440

startyr = 1998 # Don't change - tied to file names!
endyr = 2014

anstartyr = 1998 #year for analysis start
anendyr = 2014 #year for analysis end
nyears = anendyr - anstartyr + 1
nmonths = 12

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

starttsteps = 0
endtsteps = 46752 # 16 years
anntsteps = 2920 # timesteps per year

filetimespan = '3hrly'

DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
FileP = 'TRMM_' + str(startyr) + '-' + str(endyr) + '_seasmean_nonan.nc'

DirO = DirP 
FileO = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_' + FileP

#pcp, time, ongitude, latitude
precipdata = xray.open_dataset(DirP + FileP)
lons = precipdata['longitude']
lats = precipdata['latitude']
nlons = len(lons)
nlats = len(lats)

precipSeas = precipdata['pcp']
ntimespre = len(precipdata['time'])

ntimes = precipSeas.shape[0]
print ntimes

nlonsnew = nlons/sumlons
nlatsnew = nlats/sumlats

Latsnew = np.zeros(nlatsnew,np.float)
Lonsnew = np.zeros(nlonsnew,np.float)

RDenMapnew = np.zeros((ntimes,nlatsnew,nlonsnew),np.float)

#'regrid' by averaging over large boxes (averaging lons and lats)
inlat = 0
for ilats in range(0,nlats,sumlats):
	inlon = 0
	Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlats])
#	print np.mean(lats[ilats:ilats+sumlats])
	for ilons in range(0,nlons,sumlons):
		Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlons])
#		print np.mean(lons[ilons:ilons+sumlons])
		RDenMapnew[:,inlat,inlon] = np.mean(precipSeas[:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(1,2),dtype=np.float)
		inlon += 1
	inlat += 1

#Now put into different seasons
RDenMapSeas = np.zeros((nyears-1,4,nlatsnew,nlonsnew),np.float)

RDenMapSeas[:,0,:,:] = RDenMapnew[range(4,ntimes-3,4),:,:]	# do DJF separately to skip first one - as only JF.
for iseas in range(1,4):
	print range(iseas,ntimes-4,4)
	RDenMapSeas[:,iseas,:,:] = RDenMapnew[range(iseas,ntimes-4,4),:,:]	# all others include first year, skpi last

# Create Annual sum

RDenMapAnn = np.sum(RDenMapSeas,axis=1)
print RDenMapAnn.shape

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('seas', 4)
ncfile.createDimension('year',nyears-1)
ncfile.createDimension('lon', nlonsnew)
ncfile.createDimension('lat', nlatsnew)

ODensityMap = ncfile.createVariable('PrecipAnn','f4',('year','lat','lon'),fill_value=-9999)
ODensityMapSeas = ncfile.createVariable('PrecipSeas','f4',('year','seas','lat','lon'),fill_value=-9999)

OLongitude = ncfile.createVariable('Longitude','f4',('lon'),fill_value=-9999)
OLatitude = ncfile.createVariable('Latitude','f4',('lat'),fill_value=-9999)
OSeas = ncfile.createVariable('seas',str,('seas'))

OLatitude[:] = Latsnew
OLongitude[:] = Lonsnew
ODensityMap[...] = RDenMapAnn
ODensityMapSeas[...] = RDenMapSeas

Ngl.end()









