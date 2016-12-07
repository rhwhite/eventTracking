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
maxspeed = 1

direction = "westerly"

# In m/s
tbound = np.array([1000,30,10,6,3])
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


f4vars = np.array(['xmaxtravel','gridboxspanSA','totalprecipSA','uniquegridboxspanSA','gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax'])

f8vars = np.array(['eventid'])

nf4vars = len(f4vars)
nf8vars = len(f8vars)

nevents = len(datain.events)

variablesin = np.zeros([nf4vars,nevents],np.float)
for ivar in range(0,nf4vars):
	variablesin[ivar,:] = datain[f4vars[ivar]].values

"""
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
"""


if maxspeed == 1:
	zonalspeed = datain['xmaxtravel'].values
else:
	xcenterstart = datain['xcenterstart'].values
	xcenterend = datain['xcenterend'].values
	ycenterstart = datain['ycenterstart'].values
	ycenterend = datain['ycenterend'].values
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
	varforwrite = np.zeros([nf4vars+nf8vars,nevents],np.float)

	if direction == "westerly":
		filetitle = str(int(tbound[ibound+1])) + '-' + str(int(tbound[ibound])) 
	else:
		filetitle = str(int(tbound[ibound])) + '-' + str(int(tbound[ibound+1]))
	if maxspeed == 1:
	        FileO = 'Precip_MaxSpeeds_' + filetitle + 'ms_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	else:
		FileO = 'Precip_Speeds_' + filetitle + 'ms_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

	# Create new file 
	ncfile = Dataset(DirO + FileO, 'w')
	ncfile.createDimension('events', None)

	Ofilevars = []
	for ivar in range(0,nf4vars):
		Ofilevars.append(ncfile.createVariable(f4vars[ivar],'f4',('events'),fill_value=-9999))
	Ofilevars.append(ncfile.createVariable('eventid','f8',('events'),fill_value=-9999))

	writeidx = 0

        #Extract events of startlen-endlen days

	if direction == "westerly":

		for ievent in range(0,nevents):
                	if (zonalspeed[ievent] <= tbound[ibound]) and (zonalspeed[ievent] > tbound[ibound+1]):
				for ivar in range(0,nf4vars):
					varforwrite[ivar,writeidx] = variablesin[ivar][ievent]
				varforwrite[nf4vars,writeidx] = ievent
				writeidx += 1

			if ievent % 10000000 == 0:
				print ievent
	else:
                for ievent in range(0,nevents):
                        if (zonalspeed[ievent] >= tbound[ibound]) and (zonalspeed[ievent] < tbound[ibound+1]):
                                for ivar in range(0,nf4vars):
                                        varforwrite[ivar,writeidx] = variablesin[ivar][ievent]
                                varforwrite[nf4vars,writeidx] = ievent
				writeidx += 1

                        if ievent % 10000000 == 0:
                                print ievent 
			
	print writeidx

	if writeidx > 0:
		for ifvar in range(0,nf4vars+1):
			Ofilevars[ifvar][0:writeidx-1] = varforwrite[ifvar,0:writeidx-1]

	datain.close()
	ncfile.close()








