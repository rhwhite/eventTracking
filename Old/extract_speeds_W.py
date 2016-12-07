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
from scipy import stats

#Version = 'standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
Version = 'Standard'
#Version = '5th_from48'

Data = "TRMM"
filetimespan = "3hrly"

# In m/s
tbound = np.array([10,6,3])
nbounds = len(tbound)

R = 6371000     # radius of Earth in m


if Data == "CESM":
	Fstartyr = 1990
	Fendyr = 2014

	startyr = 1990
	endyr = 2011 
elif Data == "TRMM":
	Fstartyr = 1998
	Fendyr = 2014
	startyr = 1998
	endyr = 2014
elif Data == "ERAI":
	Fstartyr = 1980
	Fendyr = 2014

	startyr = 1980
	endyr = 2014

nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

starttsteps = 0
endtsteps = 46752 # 16 years
anntsteps = 2920 # timesteps per year

minevent = 100000

if Data == "TRMM":
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
        	FileI1 = 'Precip_Sizes_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
		DirO = DirI
        FileInLats = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc'

elif Data == "ERAI":
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
        FileI1 = 'Precip_Sizes_ERAI' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	DirO = DirI

elif Data == "CESM":
	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	FileI1 = 'Precip_Sizes_CESM' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	DirO = DirI

else:
	print("unexpected data type")
	exit()

#Get lons and lats
iday = 0
print FileInLats
FileIn = xray.open_dataset(FileInLats)

lats = FileIn['Latitude'].values
lons = FileIn['Longitude'].values

nlats = len(lats)
nlons = len(lons)


print DirI + FileI1
datain = xray.open_dataset(DirI + FileI1)

print(datain.coords)

gridboxspanSA = datain['gridboxspanSA'].values
totalprecipSA = datain['totalprecipSA'].values
uniquegridboxspanSA = datain['uniquegridboxspanSA'].values
gridboxspan = datain['gridboxspan'].values
totalprecip = datain['totalprecip'].values
uniquegridboxspan = datain['uniquegridboxspan'].values
timespan = datain['timespan'].values
tstart = datain['tstart'].values
tmean = datain['tmean'].values
xcenterstart = datain['xcenterstart'].values
xcenterend = datain['xcenterend'].values
ycenterstart = datain['ycenterstart'].values
ycenterend = datain['ycenterend'].values
xcentermean = datain['xcentermean'].values
ycentermean = datain['ycentermean'].values
xmin = datain['xmin'].values
xmax = datain['xmax'].values
ymin = datain['ymin'].values
ymax = datain['ymax'].values

timespan = datain['timespan'].values

nevents = timespan.shape[0]

# Calculate zonal speeds
lats1 = np.radians(lats[datain['ycenterstart']])
lats2 = np.radians(lats[datain['ycenterend']])
lons1 = np.radians(lons[datain['xcenterstart']])
lons2 = np.radians(lons[datain['xcenterend']])

latsmean = (lats1 + lats2) / 2.0

az = (np.cos(latsmean) * np.cos(latsmean) * np.power(np.sin((lons2-lons1)/2),2))
cz = 2.0 * np.arctan(np.sqrt(az),np.sqrt(1-az))

distancez = R * cz

pi = 3.14
# Calculate speed
angle = np.arctan2((lats2-lats1),(lons2-lons1))
ones = np.zeros(lons1.shape)
ones[...] = 1.0
negones = ones * -1.0

direction = np.where(lons2>=lons1,ones,negones) # True where 
zonalspeed = direction * distancez/(timespan*3.0*60.0*60.0)


print nevents
for ibound in range(0,nbounds-1):
        FileO = 'Precip_Speeds_' + str(int(tbound[ibound+1])) + '-' + str(int(tbound[ibound])) + 'ms_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

	# Create new file 
	ncfile = Dataset(DirO + FileO, 'w')
	ncfile.createDimension('events', None)

	OgridboxspanSA = ncfile.createVariable('gridboxspanSA','f4',('events'),fill_value=-9999)
	OtotalprecipSA = ncfile.createVariable('totalprecipSA','f4',('events'),fill_value=-9999)
	OuniquegridboxspanSA = ncfile.createVariable('uniquegridboxspanSA','f4',('events'),fill_value=-9999)
	Ogridboxspan = ncfile.createVariable('gridboxspan','f4',('events'),fill_value=-9999)
	Ototalprecip = ncfile.createVariable('totalprecip','f4',('events'),fill_value=-9999)
	Ouniquegridboxspan = ncfile.createVariable('uniquegridboxspan','f4',('events'),fill_value=-9999)
	Otimespan = ncfile.createVariable('timespan','f4',('events'),fill_value=-9999)
	Otstart = ncfile.createVariable('tstart','f4',('events'),fill_value=-9999)
	Otmean = ncfile.createVariable('tmean','f4',('events'),fill_value=-9999)
	Oxcenterstart = ncfile.createVariable('xcenterstart','f4',('events'),fill_value=-9999)
	Oxcenterend = ncfile.createVariable('xcenterend','f4',('events'),fill_value=-9999)
	Oycenterstart = ncfile.createVariable('ycenterstart','f4',('events'),fill_value=-9999)
	Oycenterend = ncfile.createVariable('ycenterend','f4',('events'),fill_value=-9999)
	Oxcentermean = ncfile.createVariable('xcentermean','f4',('events'),fill_value=-9999)
	Oycentermean = ncfile.createVariable('ycentermean','f4',('events'),fill_value=-9999)
	
	Ozonspeed = ncfile.createVariable('zonalspeed','f4',('events'),fill_value=-9999)

	Oxmin = ncfile.createVariable('xmin','f4',('events'),fill_value=-9999)
	Oxmax = ncfile.createVariable('xmax','f4',('events'),fill_value=-9999)
	Oymin = ncfile.createVariable('ymin','f4',('events'),fill_value=-9999)
	Oymax = ncfile.createVariable('ymax','f4',('events'),fill_value=-9999)

	Oeventid = ncfile.createVariable('eventid','f8',('events'),fill_value=-9999)

	writeidx = 0
	#Extract events of startlen-endlen days
	for ievent in range(0,nevents):
		if (zonalspeed[ievent] <= tbound[ibound]) and (zonalspeed[ievent] > tbound[ibound+1]):
			OgridboxspanSA[writeidx] = gridboxspanSA[ievent]
			OtotalprecipSA[writeidx] = totalprecipSA[ievent]
			OuniquegridboxspanSA[writeidx]  = uniquegridboxspanSA[ievent]
			Ogridboxspan[writeidx] = gridboxspan[ievent]
			Ototalprecip[writeidx] = totalprecip[ievent]
			Ouniquegridboxspan[writeidx]  = uniquegridboxspan[ievent] 
			Otimespan[writeidx] = timespan[ievent]
			Otstart[writeidx] = tstart[ievent]
			Otmean[writeidx] = tmean[ievent]
			Oxcenterstart[writeidx] = xcenterstart[ievent]
			Oxcenterend[writeidx] = xcenterend[ievent]
			Oycenterstart[writeidx] = ycenterstart[ievent]
			Oycenterend[writeidx] = ycenterend[ievent]
			Oxcentermean[writeidx] = xcentermean[ievent]
			Oycentermean[writeidx] = ycentermean[ievent]
			Ozonspeed[writeidx] = zonalspeed[ievent]
			Oxmin[writeidx] = xmin[ievent]
			Oxmax[writeidx] = xmax[ievent]
			Oymin[writeidx] = ymin[ievent]
			Oymax[writeidx] = ymax[ievent]
			Oeventid[writeidx] = ievent
			writeidx = writeidx + 1	
		
		if ievent % 10000000 == 0:
			print ievent 
			
	print writeidx

	datain.close()
	ncfile.close()








