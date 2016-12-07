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

test = 0

Data = "TRMM"
filetimespan = "3hrly"

minimalvars = 0
# In days
tbound1 = np.array([0,1,2,5])
tbound2 = np.array([1,2,5,100])
unit = "day"
nbounds = len(tbound1)

Latsplit1 = np.array([0,5,10,20,30])
Latsplit2 = np.array([5,10,20,30,50])
nlatsplit = len(Latsplit1)

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
        	FileI1 = 'Precip_Sizes_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
		DirO = DirI

elif Data == "ERAI":
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
        FileI1 = 'Precip_Sizes_ERAI' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	DirO = DirI
        DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
        FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + 'preprocess.nc'
elif Data == "CESM":
	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	FileI1 = 'Precip_Sizes_CESM' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	DirO = DirI

else:
	print("unexpected data type")
	exit()

f4vars = np.array(['xmaxspeed_4ts','xmaxtravel','gridboxspan','totalprecip','timespan','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean'])

f8vars = np.array(['eventid'])

nf4vars = len(f4vars) 

# Translate to number of timesteps
if filetimespan == "3hrly":
        tboundper = tbound2 * 8.0
	tboundstart = tbound1[0] * 8.0
print tboundper
tempfile = xray.open_dataset(DirP + FileP)
lons = tempfile['longitude'].values
lats = tempfile['latitude'].values
tempfile.close()

nlons = lons.shape
nlats = lats.shape
print nlons
print nlats


print DirI + FileI1
datain = xray.open_dataset(DirI + FileI1)

print(datain.coords)

timespan = datain['timespan'].values
xcentermean = datain['xcentermean'].values
ycentermean = datain['ycentermean'].values

nevents = timespan.shape[0]

print nevents

invars = np.zeros([len(f4vars),nevents],np.float)
#read in variables
for ivar in range(0,len(f4vars)):
	invars[ivar,:] = datain[f4vars[ivar]].values


## Define functions

def returntime(timespanin):
	
	for ibound in range(0,nbounds):
		if (timespanin < tboundper[ibound]):
			tindex = ibound
			break
	return tindex

def createfile(filenames, f4varsin, f8varsin):
	filevarsAll = []
	for filename in filenames:
        	filename.createDimension('events', None)

        	filevars = []
        	for ivar in range(0,len(f4varsin)):
                	filevars.append(filename.createVariable(f4varsin[ivar],'f4',('events'),fill_value=-9999))

	        for ivar in range(0,len(f8varsin)):
	                filevars.append(filename.createVariable(f8varsin[ivar],'f8',('events'),fill_value=-9999))
		filevarsAll.append(filevars)
        return (filevarsAll)


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


OfilevarsAll = []

ncfilesAll = []

if minimalvars == 1:
	fileOadd = "selvars"
else:
	fileOadd = ""

for ibound in range(0,nbounds):
        tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))

        ncfiles = []
        FileO = 'Precip_Sizes_' + tboundtitle + 'day_' + Data  + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	
	for ilat in range(0,nlatsplit):
		ncfile = Dataset(DirO + 'Region' + str(Latsplit1[ilat]) + "-" + str(Latsplit2[ilat]) + "N+S_" + FileO,'w')
		ncfiles.append(ncfile)
		ncfilesAll.append(ncfile)

	OfilevarsAll.append(createfile(ncfiles,f4vars,f8vars))
	ncfilesAll.append(ncfiles)

nfiles = len(ncfiles)
print nbounds
writeidxAll = np.zeros([nbounds,nlatsplit],np.int)

# Assume that, at a maximum, 75% of all events would end up in one region in the 0-1day timespan chunk!
filevarsAll = np.zeros([nbounds,nlatsplit,len(f4vars)+1,nevents*0.75],np.float)

#Extract events of startlen-endlen days
if test == 1:
	end = 10000
else:
	end = nevents

for ievent in range(0,end):
	# Get time index
	if timespan[ievent] >= tboundstart:
		tindex = returntime(timespan[ievent])
		# get basin index
		Tidx = returnbasin(xcentermean[ievent],ycentermean[ievent])
		
		# Write to files	
		# first index: ndays, second index: midlatitudes[0] vs tropics[1], third index: variable, fourth index: event
		# print most variables
		for ivar in range(0,len(f4vars)):

			filevarsAll[tindex][Tidx][ivar][writeidxAll[tindex][Tidx]] = invars[ivar,ievent]

		# print last variable, ievent
		
		filevarsAll[tindex][Tidx][nf4vars][writeidxAll[tindex][Tidx]] = ievent
		writeidxAll[tindex][Tidx] +=1

	if ievent % 100000 == 0:
		print ievent 

for itindex in range(0,nbounds):
	for iTidx in range(0,nlatsplit):
        
        	if writeidxAll[itindex][iTidx] > 0:
			for ifvar in range(0,nf4vars):
				print itindex
				print iTidx
				print ifvar
				print writeidxAll[itindex][iTidx]
	        	        
				OfilevarsAll[itindex][iTidx][ifvar][0:writeidxAll[itindex][iTidx]-1] = filevarsAll[itindex,iTidx,ifvar,0:writeidxAll[itindex][iTidx]-1]

	        	OfilevarsAll[itindex][iTidx][nf4vars][0:writeidxAll[itindex][iTidx]-1] = filevarsAll[itindex,iTidx,nf4vars,0:writeidxAll[itindex][iTidx]-1]

		
datain.close()

print writeidxAll

for ifile in range(0,nfiles):
        ncfiles[ifile].close()








