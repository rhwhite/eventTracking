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
splittype = "maxspeed" # day,speed,maxspeed
Latsplit1 = np.aray([0,10,30])
Latsplit2 = np.array([5,20,50])
nlatsplit = len(Latsplit1)
minimalvars = 0
dobasin = 0
# In days
tbound1 = np.array([-30]) #np.array([0,1,2,5,100])
tbound2 = np.array([-10])
nbounds = len(tbound1)

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
	fileadd = "_" str(speedtspan) + "ts"

tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))
if splittype == "maxspeed":
	fileadd = fileadd + "MaxSpeeds_" + add
	unit = "ms"
elif splittype == "speed":
	fileadd = fileadd + "Speeds_" + add
	unit = "ms"
elif splittype == "days":
	unit = "day"
	fileadd = fileadd + "Sizes_" + add

FileI1 = 'Precip_' + fileadd + tboundtitle + unit + '_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'


if minimalvars == 1:
	f4vars = np.array([newvar,'gridboxspan','totalprecip','timespan','tstart','tmean'])
else:
	f4vars = np.array([newvar,'gridboxspanSA','totalprecipSA','uniquegridboxspanSA','gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax'])

f8vars = np.array(['eventid'])

nf4vars = len(f4vars) 

# Translate to number of timesteps
if filetimespan == "3hrly":
        tboundper = tbound * 8.0

print tboundper
tempfile = xray.open_dataset(DirP + FileP)
lons = tempfile['longitude'].values
lats = tempfile['latitude'].values
tempfile.close()

nlons = lons.shape
nlats = lats.shape
print nlons
print nlats

print lons

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
	for ibound in range(1,nbounds):
		if (timespanin < tboundper[ibound]):
			tindex = ibound-1
			break
	else:
		tindex = nbounds-1
	return tindex

def createfile(filenames, nlatsregin, f4varsin, f8varsin):
	filevarsALL = []
	for filename in filenames:
		filename.createDimension('events', None)
		filename.createDimension('latRegions',nlatsregin)
		filename.createDimension('events', None)
		filename.createDimension('latRegions',nlatsregin)

		filevars = []
        
		for ivar in range(0,len(f4varsin)):
			filevars.append(filenameM.createVariable(f4varsin[ivar],'f4',('events'),fill_value=-9999))

		for ivar in range(0,len(f8varsin)):
			filevars.append(filenameT.createVariable(f8varsin[ivar],'f8',('events'),fill_value=-9999))
		
		filevarsALL.append(filevars)

        return (filevarsALL)


def returnbasin(xicenterin,yicenterin):
        xcenterin = lons[xicenterin]
        ycenterin = lats[yicenterin]

	for ilat in range(0,nlatsplit):
		if ((ycenterin >= Latsplit1[ilat] and ycenterin < Latsplit2[ilat]) or
			(ycenterin <= -Latsplit1[ilat] and ycenterin > -Latsplit2[ilat])):
			Tidxout = ilat
	
        if ycenterin < -15.0:
                if xcenterin >= Bastart[0][0] or xcenterin < Baend[0][0]:
                        Basidx = 0      #Pacific
                elif xcenterin >= Bastart[1][0] and xcenterin < Baend[1][0]:
                        Basidx = 1      #Atlantic
                elif xcenterin >= Bastart[2][0] and xcenterin < Baend[2][0]:
                        Basidx = 2      #Indian 
                else:
                        print ycenterin,xcenterin
                        os.error("problem with SH basins")

        elif ycenterin >= -15.0 and ycenterin <= 15.0:
                if xcenterin >= Bastart[0][1] or xcenterin < Baend[0][1]:
                        Basidx = 0      #Pacific
                elif xcenterin >= Bastart[1][1] and xcenterin < Baend[1][1]:
                        Basidx = 1      #Atlantic
                elif xcenterin >= Bastart[2][1] and xcenterin < Baend[2][1]:
                        Basidx = 2      #Indian
                else:
                        print ycenterin,xcenterin
                        os.error("problem with Tropical basin")

        elif ycenterin > 15.0:
                if xcenterin >= Bastart[0][2] or xcenterin < Baend[0][2]:
                        Basidx = 0      #Pacific
                elif xcenterin >= Bastart[1][2] and xcenterin < Baend[1][2]:
                        Basidx = 1      #Atlantic
                elif xcenterin >= Bastart[2,2] and xcenterin < Baend[2][2]:
                        Basidx = 2      #Indian
                else:
                        print ycenterin,xcenterin
                        sys.error("problem with NH basins")
        else:
                print lats[ilat]
                sys.error("problem with latitudes")
	return (Basidx,Tidxout)

# Define basins

Bastart = np.array([[135,110,105],[300,300,270],[30,30,30]]) 	# SH, Tropics, NH; Pac, Atl, Ind
Baend = np.array([[300,300,270],[30,30,30],[135,110,105]])

nregions = Bastart.shape[1]
nlatregs = Bastart.shape[0]
if nregions != 3 or nlatregs != 3:
	print nregions, 'basins', nlatregs, 'regions'
	os.err("hard-coded for 3 regions (SH, tropics, NH) and 3 basins (Pacific, Atlantic, Indian")
# Rearrange if lons from -180 to 180
if lons[0] < 0:
	for ireg in range(0,nregions):
		for ilatreg in range(0,nlatregs):
			if Bastart[ireg][ilatreg] > 180:
				Bastart[ireg][ilatreg] -= 360.0
			if Baend[ireg][ilatreg] > 180:
				Baend[ireg][ilatreg] -= 360.0
if np.any(Bastart[0] < Baend[0]):
	os.error("not set up for Pacific NOT over boundaries")
if np.any(Bastart[1] > Baend[1]):
        os.error("not set up for Atlantic over boundaries")
if np.any(Bastart[2] > Baend[2]):
        os.error("not set up for India over boundaries")

## Create new files and put these file handles in arrays

OfilevarsPa = []
OfilevarsAt = []
OfilevarsIn = []
OfilevarsAll = []

ncfiles = []

if minimalvars == 1:
	fileOadd = "selvars"
else:
	fileOadd = ""

for ibound in range(0,nbounds):
	if ibound < nbounds-1:
	        FileO = '_' + fileOadd + 'Precip_Sizes_' + str(int(tbound[ibound])) + '-' + str(int(tbound[ibound+1])) + 'day_' + Data  + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	else:
		FileO = '_' + fileOadd + 'Precip_Sizes_gt' + str(int(tbound[ibound])) + 'day_' + Data  + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

	for ilat in range(0,nlatsplit):
		(ycenterin <= -Latsplit1[ilat] and ycenterin > -Latsplit2[ilat])): 
	        ncfile = Dataset(DirO + 'AllTr' + FileO, 'w')
	ncfileM = Dataset(DirO + 'AllMid' + FileO, 'w')

	OfilevarsAll.append(createfile(ncfiles,nlatregs,f4vars,f8vars))
	ncfiles.append(ncfileT)
	ncfiles.append(ncfileM)

	if dobasin == 1:
		ncfileT = Dataset(DirO + 'PacificTr' + FileO, 'w')
		ncfileM = Dataset(DirO + 'PacificMid' + FileO, 'w')
		OfilevarsPa.append(createfile(ncfileM,ncfileT,nlatregs,f4vars,f8vars))
		ncfiles.append(ncfileT)
		ncfiles.append(ncfileM)

		ncfileT = Dataset(DirO + 'AtlanticTr' + FileO, 'w')
		ncfileM = Dataset(DirO + 'AtlanticMid' + FileO, 'w')
		OfilevarsAt.append(createfile(ncfileM,ncfileT,nlatregs,f4vars,f8vars))
		ncfiles.append(ncfileT)
		ncfiles.append(ncfileM)

		ncfileT = Dataset(DirO + 'IndianTr' + FileO, 'w')
		ncfileM = Dataset(DirO + 'IndianMid' + FileO, 'w')
		OfilevarsIn.append(createfile(ncfileM,ncfileT,nlatregs,f4vars,f8vars))
		ncfiles.append(ncfileT)
		ncfiles.append(ncfileM)

nfiles = len(ncfiles)

writeidxPa = np.zeros([nbounds,2],np.int)
writeidxAt = np.zeros([nbounds,2],np.int)
writeidxIn = np.zeros([nbounds,2],np.int)
writeidxAll = np.zeros([nbounds,2],np.int)

filevarsAll = np.zeros([nbounds-1,2,len(f4vars)+1,nevents],np.float)

#Extract events of startlen-endlen days
for ievent in range(0,nevents):
	# Get time index
	if timespan[ievent] >= tboundper[0]:
		tindex = returntime(timespan[ievent])
		# get basin index
		basindex, Tidx = returnbasin(xcentermean[ievent],ycentermean[ievent])
		
		# Write to files	
		# first index: ndays, second index: midlatitudes[0] vs tropics[1], third index: variable, fourth index: event
		# print most variables
		for ivar in range(0,len(f4vars)):
                        if dobasin == 1:
				if basindex == 0:
					filevarsPa[tindex][Tidx][ivar][writeidxPa[tindex][Tidx]] = invars[ivar,ievent]
				elif basindex == 1:
					filevarsAt[tindex][Tidx][ivar][writeidxAt[tindex][Tidx]] = invars[ivar,ievent]
				elif basindex == 2:
					filevarsIn[tindex][Tidx][ivar][writeidxIn[tindex][Tidx]] = invars[ivar,ievent]
				else:
					print basindex
					os.error("unexpected basindex")

			filevarsAll[tindex][Tidx][ivar][writeidxAll[tindex][Tidx]] = invars[ivar,ievent]

		# print last variable, ievent
		if dobasin == 1:
			if basindex == 0:
				filevarsPa[tindex][Tidx][nf4vars][writeidxPa[tindex][Tidx]] = ievent
				writeidxPa[tindex][Tidx] += 1
			elif basindex == 1:
				filevarsAt[tindex][Tidx][nf4vars][writeidxAt[tindex][Tidx]] = ievent
				writeidxAt[tindex][Tidx] += 1
			elif basindex == 2:
				filevarsIn[tindex][Tidx][nf4vars][writeidxIn[tindex][Tidx]] = ievent
				writeidxIn[tindex][Tidx] += 1
		
		filevarsAll[tindex][Tidx][nf4vars][writeidxAll[tindex][Tidx]] = ievent
		writeidxAll[tindex][Tidx] +=1

	if ievent % 100000 == 0:
		print ievent 

for itindex in range(0,nbounds-1):
	for iTidx in range(0,2):
		print writeidxAll[itindex][iTidx]
        	for ifvar in range(0,nf4vars):
        	        OfilevarsAll[itindex][iTidx][ifvar][0:writeidxAll[itindex][iTidx]-1] = filevarsAll[itindex,iTidx,ifvar,0:writeidxAll[itindex][iTidx]-1]

        	OfilevarsAll[itindex][iTidx][nf4vars][0:writeidxAll[itindex][iTidx]-1] = filevarsAll[itindex,iTidx,nf4vars,0:writeidxAll[itindex][iTidx]-1]

		
datain.close()

print writeidxAll

for ifile in range(0,nfiles):
        ncfiles[ifile].close()








