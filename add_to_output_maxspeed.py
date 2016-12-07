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
import math
import xray

Data = "TRMM"
Version = "6th_from6"
test = 0

filetimespan = "3hrly"

speedtspan = 4

if filetimespan == "3hrly":
	mult = 3.0

startyr = 1998
endyr = 2014 

R = 6371000     # radius of Earth in m

minevent = 100000
if Data == "TRMM":
        Pfilestartyr = 1998
        Pfileendyr = 2014
        DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
        FileP = 'TRMM_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_3B42_3hrly_nonan.nc'
	if Version == "Standard":
		Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
		File1 = 'ts_TRMM1998-2014_final_4Dobjects.nc'
		FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_Standard.nc'
		TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_final_4Dobject_tree.txt'
	elif Version == "7thresh":
		Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip'
		File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_7thresh_4Dobjects.nc'
		TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_7thresh_4Dobject_tree.txt'
		FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_7thresh.nc'
	elif Version == "5thresh_n2z":
		Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
		File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_5th_n2z_4Dobjects.nc'
		TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_5th_n2z_4Dobject_tree.txt'
		FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_5th_n2zero.nc'
	elif Version == "6th_from6" or Version == "5th_from48":
		Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
		File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
		TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
		FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
elif Data == "ERAI":
        Pfilestartyr = 1980
        Pfileendyr = 2015
        startyr = 1998
        endyr = 2014
        DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
        FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_preprocess.nc'

        Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
        File1 = 'ts_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
        TxtFileIn = 'ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
        FileO = 'Precip_Sizes_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
elif Data == "CESM":
        DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
        FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '.nc'
        Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) +'/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
        File1 = 'ts_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
        TxtFileIn = 'CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
        FileO = 'Precip_Sizes_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
else:
	sys.exit(Data + " Version " + Version + " doesn't exist here!")


Filegrid = Data + "_SurfaceArea.nc"


# If output directory doesn't exist, create it!
if not os.path.exists(DirO):
        os.makedirs(DirO)

# Open precip file for latitudes
print DirP + FileP
tempfile = xray.open_dataset(DirP + FileP)
if Data == "CESM":
    lons = tempfile['lon'].values
    lats = tempfile['lat'].values
else:
    lons = tempfile['longitude'].values
    lats = tempfile['latitude'].values

nlats = lats.shape[0]
nlons = lons.shape[0]

textdata = np.loadtxt(Dir + TxtFileIn,skiprows = 1, usecols = (0,))
print textdata.shape
print textdata[0]
print textdata[-1]

event1 = textdata[0]
nevents = textdata[-1] - textdata[0] + 1
print("nevents: ", nevents)


minx = np.zeros(nevents)
maxx = np.zeros(nevents)
miny = np.zeros(nevents)
maxy = np.zeros(nevents)
mint = np.zeros(nevents)
maxt = np.zeros(nevents)
maxzonalspeed = np.zeros(nevents)

#centerxstart = np.zeros(nevents)
#centerystart = np.zeros(nevents)

centerx = np.zeros(nevents)
centery = np.zeros(nevents)
meant = np.zeros(nevents)

nlines = len(textdata)
print("nlines: ", nlines)

preeventnum = -100

with open(Dir + TxtFileIn,"r") as textFile:
	next(textFile)	#Skip header line

	for lines in textFile:	#loop through all lines
		line = lines.split('\t')
	        eventnum = int(float(line[0]))
        	if (eventnum % 1000000 == 0):
                	print "eventnum: " + str(eventnum)
		if eventnum == preeventnum:
			listcenterx.append(int(re.findall(r'\d+',line[6])[0]))
			listcentery.append(int(re.findall(r'\d+',line[6])[1]))
			eventcount += 1
		else:
			if (preeventnum > 0): #then it's not the very first event
				#new event: take mean of center x, center y and t and put into array before updating index
				if eventcount == 1:
					maxzonalspeed[index] = 0.0	# Stationary event by definition of one timestep
				elif eventcount <= speedtspan:
					diffs = listcenterx[eventcount-1] - listcenterx[0]
					maxzonalspeed[index] = np.cos(np.radians(lats[listcentery[int(eventcount/2 -1)]])) * R * 2.0 * math.pi * diffs/(nlons * 3.0 * eventcount * 60.0 * 60.0)   # speed in m/s
				elif eventcount > speedtspan:
					diffs = np.array([x - listcenterx[i-speedtspan] for i,x in enumerate(listcenterx)][speedtspan:])
					coslat = np.cos(np.radians(lats[listcentery[speedtspan:]]))
					maxindex = np.argmax(abs(diffs)*coslat)
					maxzonalspeed[index] = np.cos(np.radians(lats[listcentery[maxindex]])) * R * 2.0 * math.pi * diffs[maxindex]/(nlons * 3.0 * speedtspan * 60.0 * 60.0)	# speed in m/s
				else:
					sys.exit("eventcount is less than 0?")
			#refresh list
			#update index for new event
			index = eventnum-event1
			preeventnum = eventnum
			
                        listcenterx = [int(re.findall(r'\d+',line[6])[0])]
                        listcentery = [int(re.findall(r'\d+',line[6])[1])]
			eventcount = 1

		if test == 1:
			if index > 10000:
				break

textFile.close()

print 'maxzonalspeed'
print maxzonalspeed[0:10]

eventsdata = xray.open_dataset(Dir + File1)

print(eventsdata.coords)

eventsin=eventsdata['value']

eventsin = eventsin.squeeze(dim="z")

ntimes = eventsin.shape[0]
nlats = eventsin.shape[1]
nlons = eventsin.shape[2]
maxevent = int(np.amax(eventsin[ntimes-10:ntimes,:,:]))

nevents2 = maxevent - minevent + 1

if (nevents2 != nevents):
	print nevents2, nevents
	sys.exit("event numbers not equal in netcdf and text file")

if test == 1:
	ncfile = Dataset(DirO + 'test' + FileO, 'w')
else:
	print DirO + FileO
        ncfile = Dataset(DirO + FileO, 'r+')

varname = 'xmaxspeed_' + str(speedtspan) + 'ts' 
try:
	Omaxzonalspeed = ncfile.createVariable(varname,'f4',('events'),fill_value=-9999)
	setattr(Omaxzonalspeed, 'Extra Info', 'maximum zonal speed during lifetime for 3 * ' + str(speedtspan) + 'hours, for non-zero values; minimum value is ~2.5m/s at equator based on spatial and temporal resolution')
except RuntimeError:
	print varname, "already in file"
	Omaxzonalspeed = ncfile[varname]
	print Omaxzonalspeed

setattr(Omaxzonalspeed, 'Extra Info', 'maximum zonal speed during lifetime for 3 * ' + str(speedtspan) + 'hours, for non-zero values; minimum value is ~2.5m/s at equator based on spatial and temporal resolution')
print "starting now"

#vartowrite = np.zeros(nevents,np.float)

if test == 1:
	endloop = minevent + 10
else:
	endloop = maxevent + 1

#for ievent in range(minevent,endloop):

#        index = ievent-minevent
#	if (ievent % 10000 == 0):
#		print "ievent: " + str(ievent)

#	vartowrite[ievent-minevent] = maxzonalspeed[index]

Omaxzonalspeed[:] = maxzonalspeed[:]

#results = results2
ncfile.close()

