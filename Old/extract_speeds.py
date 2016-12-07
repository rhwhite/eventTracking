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
speedtspan = 4

# In m/s
tbound = np.array([-30,-6]) #-1000,-30,-10,-6,-3])
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

if maxspeed == 1:
	if speedtspan == 0:
		newvar = "xmaxtravel"
		newvarout = "xmaxspeed"
	else:
		newvar = "xmaxspeed_" + str(speedtspan) + "ts"
		newvarout = newvar

f4vars = np.array([newvarout,'gridboxspanSA','totalprecipSA','uniquegridboxspanSA','gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax'])
f4varsin = np.array([newvar,'gridboxspanSA','totalprecipSA','uniquegridboxspanSA','gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax'])

f8vars = np.array(['eventid'])

nf4vars = len(f4vars)
nf8vars = len(f8vars)

nevents = len(datain.events)

variablesin = np.zeros([nf4vars,nevents],np.float)
for ivar in range(0,nf4vars):
	variablesin[ivar,:] = datain[f4varsin[ivar]].values


if maxspeed == 1:
	if speedtspan == 0:
		zonalspeed = datain['xmaxtravel'].values
	else:
		zonalspeed = datain['xmaxspeed_' + str(speedtspan) +  'ts'].values
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


varforwrite = np.zeros([nf4vars+nf8vars,nevents],np.float)

print nevents
for ibound in range(0,nbounds-1):
        if speedtspan == 0:
		fileadd = ""
	else:
		fileadd = str(speedtspan) + "ts_"

	if tbound[ibound] < 0:
		tboundtitle = str(int(tbound[ibound])) + '-' + str(int(tbound[ibound+1]))
	else:
		tboundtitle = str(int(tbound[ibound+1])) + '-' + str(int(tbound[ibound]))
	if maxspeed == 1:
	        FileO = 'Precip_' + fileadd + 'MaxSpeeds_' + tboundtitle + 'ms_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	else:
		FileO = 'Precip_Speeds_' + tboundtitle + 'ms_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

	# Create new file 
	ncfile = Dataset(DirO + FileO, 'w')
	ncfile.createDimension('events', None)

	Ofilevars = []
	for ivar in range(0,nf4vars):
		Ofilevars.append(ncfile.createVariable(f4vars[ivar],'f4',('events'),fill_value=-9999))
	"""
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

	if maxspeed == 1:	
		Ozonspeed = ncfile.createVariable('maxzonalspeed','f4',('events'),fill_value=-9999)
	else:
                Ozonspeed = ncfile.createVariable('zonalspeed','f4',('events'),fill_value=-9999)

	Oxmin = ncfile.createVariable('xmin','f4',('events'),fill_value=-9999)
	Oxmax = ncfile.createVariable('xmax','f4',('events'),fill_value=-9999)
	Oymin = ncfile.createVariable('ymin','f4',('events'),fill_value=-9999)
	Oymax = ncfile.createVariable('ymax','f4',('events'),fill_value=-9999)
	"""
	Ofilevars.append(ncfile.createVariable('eventid','f8',('events'),fill_value=-9999))

	writeidx = 0


	#Extract events of startlen-endlen days
	if tbound[ibound] < 0:
		for ievent in range(0,nevents):
			if (zonalspeed[ievent] >= tbound[ibound]) and (zonalspeed[ievent] < tbound[ibound+1]):
				for ivar in range(0,nf4vars):
					varforwrite[ivar,writeidx] = variablesin[ivar][ievent]
				varforwrite[nf4vars,writeidx] = ievent
				
				writeidx += 1
			if ievent % 10000000 == 0:
				print ievent 
	else:
                for ievent in range(0,nevents):
                        if (zonalspeed[ievent] <= tbound[ibound]) and (zonalspeed[ievent] > tbound[ibound+1]):
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








