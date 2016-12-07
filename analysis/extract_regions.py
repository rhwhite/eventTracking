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

f4vars = np.array(['gridboxspanSA','totalprecipSA','uniquegridboxspanSA','gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax'])

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

def createfile(filenameM, nlatsregin, f4varsin, f8varsin):

        filenameM.createDimension('events', None)
        filenameM.createDimension('latRegions',nlatsregin)

        filevarsM = []
        for ivar in range(0,len(f4varsin)):
                filevarsM.append(filenameM.createVariable(f4varsin[ivar],'f4',('events'),fill_value=-9999))

        for ivar in range(0,len(f8varsin)):
                filevarsM.append(filenameM.createVariable(f8varsin[ivar],'f8',('events'),fill_value=-9999))

        return (filevarsM)


def returnbasin(xicenterin,yicenterin):
        xcenterin = lons[xicenterin]
        ycenterin = lats[yicenterin]
        if ycenterin < -15.0:
                Tidxout = 0        # Midlatitudes
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
                Tidxout = 1        # Tropics         
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
                Tidxout = 0
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

filevarsPa = []
filevarsAt = []
filevarsIn = []
filevarsAll = []

ncfiles = []

FileO = '_Precip_Sizes_' + Data  + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

ncfileT = Dataset(DirO + 'AllTr' + FileO, 'w')
ncfileM = Dataset(DirO + 'AllMid' + FileO, 'w')
filevarsAll.append(createfile(ncfileM,nlatregs,f4vars,f8vars))
filevarsAll.append(createfile(ncfileT,nlatregs,f4vars,f8vars))
ncfiles.append(ncfileT)
ncfiles.append(ncfileM)

ncfileT = Dataset(DirO + 'PacificTr' + FileO, 'w')
ncfileM = Dataset(DirO + 'PacificMid' + FileO, 'w')
filevarsPa.append(createfile(ncfileM,nlatregs,f4vars,f8vars))
filevarsPa.append(createfile(ncfileT,nlatregs,f4vars,f8vars))

ncfiles.append(ncfileT)
ncfiles.append(ncfileM)


ncfileT = Dataset(DirO + 'AtlanticTr' + FileO, 'w')
ncfileM = Dataset(DirO + 'AtlanticMid' + FileO, 'w')
filevarsAt.append(createfile(ncfileM,nlatregs,f4vars,f8vars))
filevarsAt.append(createfile(ncfileT,nlatregs,f4vars,f8vars))

ncfiles.append(ncfileT)
ncfiles.append(ncfileM)

ncfileT = Dataset(DirO + 'IndianTr' + FileO, 'w')
ncfileM = Dataset(DirO + 'IndianMid' + FileO, 'w')
filevarsIn.append(createfile(ncfileM,nlatregs,f4vars,f8vars))
filevarsIn.append(createfile(ncfileT,nlatregs,f4vars,f8vars))
ncfiles.append(ncfileT)
ncfiles.append(ncfileM)

nfiles = len(ncfiles)

writeidxPa = np.zeros(2,np.int)
writeidxAt = np.zeros(2,np.int)
writeidxIn = np.zeros(2,np.int)
writeidxAll = np.zeros(2,np.int)

#Extract events of startlen-endlen days
for ievent in range(0,nevents):
	# get basin index
	basindex, Tidx = returnbasin(xcentermean[ievent],ycentermean[ievent])
	
	# Write to files	
	# first index: ndays (not used here), second index: midlatitudes[0] vs tropics[1], third index: variable, fourth index: event
	# print most variables
	for ivar in range(0,len(f4vars)):
		if basindex == 0:
			filevarsPa[Tidx][ivar][writeidxPa[Tidx]] = invars[ivar,ievent]
		elif basindex == 1:
			filevarsAt[Tidx][ivar][writeidxAt[Tidx]] = invars[ivar,ievent]
		elif basindex == 2:
			filevarsIn[Tidx][ivar][writeidxIn[Tidx]] = invars[ivar,ievent]
		else:
			print basindex
			os.error("unexpected basindex")

		filevarsAll[Tidx][ivar][writeidxAll[Tidx]] = invars[ivar,ievent]

	# print last variable, ievent
	if basindex == 0:
		filevarsPa[Tidx][nf4vars][writeidxPa[Tidx]] = ievent
		writeidxPa[Tidx] += 1
	elif basindex == 1:
		filevarsAt[Tidx][nf4vars][writeidxAt[Tidx]] = ievent
		writeidxAt[Tidx] += 1
	elif basindex == 2:
		filevarsIn[Tidx][nf4vars][writeidxIn[Tidx]] = ievent
		writeidxIn[Tidx] += 1
	
	filevarsAll[Tidx][nf4vars][writeidxAll[Tidx]] = ievent
	writeidxAll[Tidx] +=1

	if ievent % 10000 == 0:
		print ievent 
		
datain.close()

print writeidxAll

for ifile in range(0,nfiles):
        ncfiles[ifile].close()








