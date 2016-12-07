# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os
import errno
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

test = 0
mapprecip = 1
splittype = "maxspeed"	# "speed","maxspeed"
speedtspan = 4
# In m/s
if splittype == "speed" or splittype == "maxspeed":
        tbound = np.array([-1000,-30,-10,-6,3,6,10,30])# [0,1,2,5])
        tbound2 = np.array([-30,-10,-6,-3,6,10,30,1000])# [1,2,5,100]
        unit = "ms"
elif splittype == "days":
        tbound = np.array([0,1,2])#,5])#([0,1,2,5])
        tbound2 = np.array([1,2,5])#,100])#([1,2,5,100]
        unit = "day"
nbounds = len(tbound)

Data = "TRMM"

mapping = 'centre'
Version = 'Standard'
# Version = '5th_nanto25'
# Version = '5th_nantozero'
# Version = '7thresh'
# Version = '6th_from6'
# Version = '5th_from48'

print Version

R = 6371000     # radius of Earth in m

filetimespan = "3hrly"

if filetimespan == "3hrly":
        if Data == "ERAI" or Data == "TRMM" or Data == "TRMM_ERAIgd":
                # convert from mm/hour to mm/3 hours to get total rain over event for 3-hourly data
                mult = 3.0
        elif Data == "CESM":
                # convert from m/s to mm/3 hours to get total rain over event for 3-hourly data
                mult = 1000.0 * 60.0 * 60.0 * 3.0
        else:
                sys.error(Data + " not defined")

if test == 1:
        chunksize = 50
else:
        chunksize = 2000

Seasons = ["MAM", "JJA", "SON", "DJF"]
seaststeps = [(31+30+31)*8, (30+31+31)*8, (30+31+30)*8, (31+31+28)*8]

minevent = 100000

if Data == "TRMM":
    anstartyr = 1998    # year for analysis start
    anendyr = 2014      # year for analysis end

    startyr = 1998
    endyr = 2014
    DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
    FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'
    Filegrid = "SurfaceArea.nc"
    if Version in ["Standard", "7thresh", "6th_from6", "5th_from48"]:
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
        DirO = DirI
        DirEv = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/'
        FileEv = 'ts_TRMM1998-2014_final_4Dobjects.nc'

    else:
        sys.exit('unexpected Version')
elif Data == "ERAI":
    anstartyr = 1980    # year for analysis start
    anendyr = 2014      # year for analysis end

    Pfilestartyr = 1980
    Pfileendyr = 2015
    startyr = 1980
    endyr = 2014
    DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
    FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_preprocess.nc'

    DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
    DirO = DirI
    DirEv = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/'
    FileEv = 'ts_ERAI' + str(startyr) + '-' + str(endyr) + '_Standard_4Dobjects.nc'
    Filegrid = "SurfaceArea.nc"

elif Data == "CESM":
    startyr = 1990   # Don't change - tied to file names!
    endyr = 2014

    Pfilestartyr = 1979
    Pfileendyr = 2012
    DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
    FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '.nc'

    DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
    DirO = DirI
else:
    sys.error("unexpected datatype")

if mapprecip == 0:
    add = "noP"
else:
    add = ""

Seasons = ["MAM", "JJA", "SON", "DJF"]


def initialize(name, numlats, numlons, atype):
    # globals used to create global object
    globals()[name] = np.zeros((numlats, numlons), atype)


def resetvar(name):
    globals()[name][...] = 0



for ibound in range(0,nbounds):
	print  str(int(tbound[ibound])) + '-' + str(int(tbound2[ibound])) + 'ms_'
        if speedtspan == 0:
                fileadd = ""
        else:   
                fileadd = str(speedtspan) + "ts_"

        tboundtitle = str(int(tbound[ibound])) + '-' + str(int(tbound2[ibound]))
	if splittype == "maxspeed":
		fileadd = fileadd + "MaxSpeeds_" + add
	elif splittype == "speed":
		fileadd = fileadd + "Speeds_" + add
	elif splittype == "days":
		fileadd = fileadd + "Sizes_" + add

        FileI1 = 'Precip_' + fileadd + tboundtitle + unit + '_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

	if test == 0:
	    FileO = 'DenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
	else:
	    FileO = 'testDenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

	nyears = anendyr - anstartyr
	nseas = 4
	# Timestep to start at: Assuming file starts at January, this is the first March
	starttsteps = (31+28)*8
	# Make sure we have the same number of years in each season.
	ntsteps = nyears * 365 * 8

	mints = np.zeros([nyears, nseas], np.int)
	maxts = np.zeros([nyears, nseas], np.int)

	filetimespan = "3hrly"
	print DirP + FileP
	tempfile = xray.open_dataset(DirP + FileP)
	if Data == "CESM":
	    lons = tempfile['lon'].values
	    lats = tempfile['lat'].values
	else:
	    lons = tempfile['longitude'].values
	    lats = tempfile['latitude'].values

	nlons = len(lons)
	nlats = len(lats)

	eventsdata = xray.open_dataset(DirEv + FileEv)
	print(eventsdata.coords)
	eventsin = eventsdata['value']

	print 'eventsin shape, ', eventsin.shape
	ntimes = eventsdata.time.size
	print DirP + FileP
	precipdata = xray.open_dataset(DirP + FileP)
	if Data == "TRMM" or Data == "TRMM_ERAIgd":
		precipin = precipdata['pcp']
	elif Data == "ERAI":
		precipin = precipdata['tpnew']
	elif Data == "CESM":
		precipin = precipdata['PRECT']

	ntimespre = len(precipdata['time'])

	# Create new datasets

	print DirO + FileO
	ncfile = Dataset(DirO + FileO, 'r+')

	if splittype == "maxspeed":
		allvars = np.array(['ExmaxSpeed','WxmaxSpeed','TPrecip','EDensity','WDensity','SDensity','ESize','WSize','SSize'])
	elif splittype == "days" or splittype == "speed":
                allvars = np.array(['EPrecip','WPrecip','TPrecip','EDensity','WDensity','SDensity','ESize','WSize','SSize'])

	
	extravars = np.array(['Longitude','Latitude'])

	nvars = len(allvars)
	print nvars
	nXvars = len(extravars)
	Ofilevars = []
	Oextravars = []
	Oextracount = 0
	for ivar in range(0,nXvars):
		if extravars[ivar] == "Longitude":	
			varhere = ncfile.createVariable(extravars[ivar], 'f8', ('lon'), fill_value=-9999)
			varhere[:] = lons[:]
		elif extravars[ivar] == "Latitude":
                        varhere = ncfile.createVariable(extravars[ivar], 'f8', ('lat'), fill_value=-9999)
			varhere[:] = lats[:]
	
		Oextracount += 1
	ncfile.close()



