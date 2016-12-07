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
splittype = "day" # day,speed,maxspeed
speedtspan = 0
if splittype == "speed" or splittype == "maxspeed":
        tbound1 = np.array([3,6,10])# [0,1,2,5])
        tbound2 = np.array([6,10,30])# [1,2,5,100]
        unit = "ms"
elif splittype == "day":
        tbound1 = np.array([0,1,2,5])
        tbound2 = np.array([1,2,5,100])
        unit = "day"
nbounds = len(tbound1)


Latsplit1 = np.array([0,5,10,20,30])
Latsplit2 = np.array([5,10,20,30,50])
nlatsplit = len(Latsplit1)
minimalvars = 0
dobasin = 0
# In days

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
        DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
        FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
		DirO = DirI

elif Data == "ERAI":
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	DirO = DirI
        DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
        FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + 'preprocess.nc'
elif Data == "CESM":
	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	DirO = DirI

else:
	print("unexpected data type")
	exit()

if speedtspan == 0:
	fileadd = ""
else:
	fileadd = "_" + str(speedtspan) + "ts"

if splittype == "maxspeed":
	fileadd = fileadd + "_MaxSpeeds_"
	unit = "ms"
        newvar = "xmaxspeed_" + str(speedtspan) + "ts"
        newvarout = newvar
elif splittype == "speed":
	fileadd = fileadd + "_Speeds_"
	unit = "ms"
        newvar = "xmaxtravel"
        newvarout = "xmaxspeed"
elif splittype == "day":
	unit = "day"
	fileadd = fileadd + "_Sizes_" 
        newvar = "xmaxtravel"
        newvarout = "xmaxspeed"

if minimalvars == 1:
	f4vars = np.array([newvar,'gridboxspan','totalprecip','timespan','tstart','tmean'])
else:
	f4vars = np.array([newvar,'gridboxspanSA','totalprecipSA','uniquegridboxspanSA','gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax'])

f8vars = np.array(['eventid'])

nf4vars = len(f4vars) 

# Translate to number of timesteps

tempfile = xray.open_dataset(DirP + FileP)
lons = tempfile['longitude'].values
lats = tempfile['latitude'].values
tempfile.close()

nlons = lons.shape
nlats = lats.shape
print nlons
print nlats

## Define functions

def createfile(filename, f4varsin, f8varsin):
	filevarsALL = []
	filename.createDimension('events', None)
	filevars = []
	for ivar in range(0,len(f4varsin)):
		filevars.append(filename.createVariable(f4varsin[ivar],'f4',('events'),fill_value=-9999))
	for ivar in range(0,len(f8varsin)):
		filevars.append(filename.createVariable(f8varsin[ivar],'f8',('events'),fill_value=-9999))

        return (filevars)


def returnbasin(xicenterin,yicenterin):
        xcenterin = lons[xicenterin]
        ycenterin = lats[yicenterin]
	Tidxout = -1
	for ilat in range(0,nlatsplit):
		if ((ycenterin >= Latsplit1[ilat] and ycenterin < Latsplit2[ilat]) or
			(ycenterin <= -1.0 * Latsplit1[ilat] and ycenterin > -1.0 * Latsplit2[ilat])):
				Tidxout = ilat
	if Tidxout == -1:
		print ycenterin
		exit("problem here")
	return (Tidxout)

## Create new files and put these file handles in arrays

for ibound in range(0,nbounds):
	tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))
	FileI1 = 'Precip' + fileadd + tboundtitle + unit + '_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

	OfilevarsAll = []
	ncfiles = []

	#Read in data from infile
	print DirI + FileI1
	datain = xray.open_dataset(DirI + FileI1)

	timespan = datain['timespan'].values
	xcentermean = datain['xcentermean'].values
	ycentermean = datain['ycentermean'].values

	nevents = timespan.shape[0]
	print nevents
	invars = np.zeros([len(f4vars),nevents],np.float)
	#read in variables
	for ivar in range(0,len(f4vars)):
		invars[ivar,:] = datain[f4vars[ivar]].values


	# Create output files
	for ilat in range(0,nlatsplit):
	        ncfile = Dataset(DirO + 'Region' + str(Latsplit1[ilat]) + "-" + str(Latsplit2[ilat]) + "N+S_" + FileI1, 'w')
	        ncfiles.append(ncfile)

		OfilevarsAll.append(createfile(ncfile,f4vars,f8vars))
	nfiles = len(ncfiles)

	writeidxAll = np.zeros([nlatsplit],np.int)

	filevarsAll = np.zeros([nlatsplit,len(f4vars)+1,nevents],np.float)

	#Extract events of startlen-endlen days
	for ievent in range(0,nevents):
		# Get time index
		# Write to files	
		# first index: latitude section, second index: variable, third index: event
                Tidx = returnbasin(xcentermean[ievent],ycentermean[ievent])

		for ivar in range(0,len(f4vars)):

			filevarsAll[Tidx][ivar][writeidxAll[Tidx]] = invars[ivar,ievent]

		# print last variable, ievent
		
		filevarsAll[Tidx][nf4vars][writeidxAll[Tidx]] = ievent
		writeidxAll[Tidx] +=1

		if ievent % 100000 == 0:
			print ievent 
	for iTidx in range(0,nlatsplit):
		if writeidxAll[iTidx] > 0:
			for ifvar in range(0,nf4vars):
				OfilevarsAll[iTidx][ifvar][0:writeidxAll[iTidx]-1] = filevarsAll[iTidx,ifvar,0:writeidxAll[iTidx]-1]
		
		
			OfilevarsAll[iTidx][nf4vars][0:writeidxAll[iTidx]-1] = filevarsAll[iTidx,nf4vars,0:writeidxAll[iTidx]-1]

			
	datain.close()

	print writeidxAll

	for ifile in range(0,nfiles):
		ncfiles[ifile].close()








