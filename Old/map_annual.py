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

Data = "TRMM"
Version = "Standard"
startyr = 1998 # Don't change - tied to file names!
endyr = 2014

Pfilestartyr = 1998
Pfileendyr = 2015

filetimespan = "3hrly"

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

minevent = 100000

if Data == "TRMM":
	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
	DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
	DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
	FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'

        FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
	FileI2 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

elif Data == "ERAI":
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/ERAI_output/Standard/Precip/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/ERAI_output/Standard/Precip/'
        DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'

	FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + 'preprocess.nc'
	FileI1 = 'Precip_Sizes_ERAI' + str(startyr) + '-' + str(endyr) + '_Standard.nc'
        FileI2 = 'Precip_Sizes_ERAI' + str(startyr) + '-' + str(endyr) + '_Standard.nc'

FileO = 'test2_DensityMap_' + str(startyr) + '-' + str(endyr) + '.nc'


filetimespan = "3hrly"
print DirP + FileP
tempfile = xray.open_dataset(DirP + FileP)
lons = tempfile['longitude'].values
lats = tempfile['latitude'].values

nlons = len(lons)
nlats = len(lats)


#Create new datasets
RDenMap = np.zeros((nyears,4,nlats,nlons),np.int)

print DirO + FileO
ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('years', nyears)
ncfile.createDimension('months',nmonths)
ncfile.createDimension('size', 4)
ncfile.createDimension('lon', nlons)
ncfile.createDimension('lat', nlats)


DensityMap = ncfile.createVariable('DensityMap','f4',('years','size','lat','lon'),fill_value=-9999)
Longitude = ncfile.createVariable('Longitude','f4',('lon'),fill_value=-9999)
Latitude = ncfile.createVariable('Latitude','f4',('lat'),fill_value=-9999)
Years = ncfile.createVariable('years','f4',('years'),fill_value=-9999)
Months = ncfile.createVariable('months','f4',('months'),fill_value=-9999)


setattr(DensityMap,'Extra Info','Size based on timespan: 0 is < 8 (1 day), 1: < 16 (2 days) < 48 (6 days), 2: > 48)')
 
Longitude[:] = lons
Latitude[:] = lats

print DirI + FileI1
datain = xray.open_dataset(DirI + FileI1)
datain2 = xray.open_dataset(DirI + FileI2)

print(datain.coords)

tspan=datain['timespan'].values
tmean = datain['tmean']
tstart = datain2['tstart'].values
xmin1 = datain['xmin']
xmin2 = datain2['xmin']

xcenter = datain['xcentermean'].values
ycenter = datain['ycentermean'].values

if (xmin1 != xmin2).any():
	sys.exit('xmin from two files not equal')

nevents = tspan.shape[0]
print nevents

ievent = 1

curtime = tstart[starttsteps]
n = 0
iyear = 0
for n in range(0,nevents):
	curtime = tstart[n]
	if curtime == starttsteps:
		mints[iyear] = n
		break

for iyear in range(0,nyears):
	Years[iyear] = anstartyr + iyear
	for n in range(int(mints[iyear]),nevents):
		curtime = tstart[n]
		if curtime == (iyear + 1) * anntsteps:
			if iyear < nyears-1:
				mints[iyear+1] = n
			
			maxts[iyear] = n  #switched to n instead of n-1 because python doesn't use the last index in a range
			break

print mints
print maxts

daycount01 = 0
daycount12 = 0
daycount25 = 0
daycount5p = 0

for iyear in range(0,nyears):
	print mints[iyear]
	print maxts[iyear]
	for ievent in range(int(mints[iyear]),int(maxts[iyear])):

		if (ievent % 100000 == 0):
			print "ievent: " + str(ievent)	
		if tspan[ievent] <= 8:
			RDenMap[iyear,0,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
			daycount01 += 1
		elif tspan[ievent] <= 16:
                        RDenMap[iyear,1,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
			daycount12 += 1
                elif tspan[ievent] <= 40:
                        RDenMap[iyear,2,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
			daycount25 +=1
		else:
                        RDenMap[iyear,3,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
			daycount5p +=1

print '0-1day events: ', daycount01
print '1-2day events: ', daycount12
print '2-5day events: ', daycount25
print '5+ day events: ', daycount5p

DensityMap[:,:,:,:] = RDenMap

datain.close()
datain2.close()

Ngl.end()









