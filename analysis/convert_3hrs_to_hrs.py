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
import argparse
from rhwhitepackages.readwrite import getunitsdesc
from rhwhitepackages.readwrite import XrayOpen

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')

args = parser.parse_args()
print "here's what I have as arguments: ", args

if args.Data[0] not in ['TRMM','ERAI','CESM']:
        exit("incorrect Data option " + str(args.Data[0]) + " must be TRMM, ERAI, or CESM")

Data = args.Data[0]
Version = args.Version[0]
startyr = args.startyr[0]
endyr = args.endyr[0]

nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

minevent = 100000

if Data == "TRMM":
        Fstartyr = 1998
        Fendyr = 2014
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
        	FileI1 = 'Precip_Sizes_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
		DirO = DirI
	elif Version == 'ERAIgd':
                DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/TRMM_output/' + Version + '/Precip/'
                FileI1 = 'Precip_Sizes_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
                DirO = DirI
        FileInLats = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc'

elif Data == "ERAI":
        Fstartyr = 1980
        Fendyr = 2015
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
        FileI1 = 'Precip_Sizes_ERAI' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	DirO = DirI
	FileInLats = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/SeasAnn_ERAI_Totalprecip_1980-2015_preprocess.nc' 
elif Data == "CESM":
        Fstartyr = 1990
        Fendyr = 2014
	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	FileI1 = 'Precip_Sizes_CESM' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
	DirO = DirI
	FileInLats = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'
else:
	print("unexpected data type")
	exit()

#Get lons and lats
iday = 0
print FileInLats
FileIn = XrayOpen(FileInLats)

if Data == "CESM":
	lats = FileIn['lat'].values
	lons = FileIn['lon'].values
else:
	lats = FileIn['Latitude'].values
	lons = FileIn['Longitude'].values

nlats = len(lats)
nlons = len(lons)


print DirI + FileI1
datain = Dataset(DirI + FileI1,'r+')

varin = datain['timespan'][...]
print varin[0:10]
if np.any(varin == 1):
	print "can't have one hour in 3hrly data, need to convert to hours"
	varout = varin * 3
	datain['timespan'][...] = varout[...]
datain.close()








