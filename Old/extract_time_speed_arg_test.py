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
from rhwhitepackages.readwhite import getunitdesc

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs=1,help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound',type=int,nargs='+',help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=int,nargs="+",help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--filetspan',type=str,nargs='?',default=['3hrly'],help='string for file time resolution, 3hrly etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--SA',type=int,nargs='?',default=[1],help='do surface area variables (1) or not (0)?')


args = parser.parse_args()

print "here's what I have as arguments: ", args

if args.splittype[0] not in ["day","speed","maxspeed"]:
        exit("incorrect splittype " + str(args.splittype[0]) + " must be day, speed, or maxspeed")
if args.speedtspan[0] not in [0,1,4]:
        exit("incorrect speedtspan " + str(args.speedtspan[0]) + " must be 0 or 1")
if args.Data[0] not in ['TRMM','ERAI','CESM']:
        exit("incorrect Data option " + str(args.Data[0]) + " must be TRMM, ERAI, or CESM")

splittype = args.splittype[0]
speedtspan = args.speedtspan[0]
tbound1 = args.tbound1
tbound2 = args.tbound2
unit = args.unit[0]
Data = args.Data[0]
Version = args.Version[0]
filetimespan = args.filetspan[0]
startyr = args.startyr[0]
endyr = args.endyr[0]

print tbound1
print filetimespan

SA = args.SA


def runchunk(nevents,chunksize,writeidxchunk):
	for chunkevent in range(0,nevents,chunksize):
		chunkmax = np.amin([chunkevent+chunksize,nevents])
		variablesin = np.zeros([nf4vars,chunksize],np.float)

		for ivar in range(0,nf4vars):
			variablesin[ivar,:] = datain[f4varsin[ivar]][chunkevent:chunkmax].values

		print "chunksize: ", chunksize, ", chunkevent: ",chunkevent
		# Loop attempts 10 times to do this, multplying chunksize by 0.25 it fails on memory, by last attempt chunksize will be ~0.000001 of original
		varforwrite = np.zeros([nf4vars+nf8vars,chunksize],np.float)
		writeidx = 0
		for ievent in range(chunkevent,chunkmax):
			if tbound1[ibound] < 0:
				if (splitvar[ievent] >= tbound1[ibound]) and (splitvar[ievent] < tbound2[ibound]):
					for ivar in range(0,nf4vars):
						varforwrite[ivar,writeidx] = variablesin[ivar][ievent-chunkevent]
					varforwrite[nf4vars,writeidx] = ievent
					writeidx += 1
			else:   
				if (splitvar[ievent] > tbound1[ibound]) and (splitvar[ievent] <= tbound2[ibound]):
					for ivar in range(0,nf4vars):
						varforwrite[ivar,writeidx] = variablesin[ivar][ievent]
					varforwrite[nf4vars,writeidx] = ievent
					writeidx += 1

		if writeidx > 0:
			if writeidx == 1:
				for ifvar in range(0,nf4vars+1):
					Ofilevars[ifvar][writeidxchunk] = varforwrite[ifvar,0]
			else:
				for ifvar in range(0,nf4vars+1):
					Ofilevars[ifvar][writeidxchunk:writeidxchunk+writeidx] = varforwrite[ifvar,0:writeidx]
		writeidxchunk += writeidx

nbounds = len(tbound1)

R = 6371000     # radius of Earth in m

nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

minevent = 100000

if Data == "TRMM":
        Fstartyr = 1998
        Fendyr = 2014
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
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
FileIn = xray.open_dataset(FileInLats)

if Data == "CESM":
	lats = FileIn['lat'].values
	lons = FileIn['lon'].values
else:
	lats = FileIn['Latitude'].values
	lons = FileIn['Longitude'].values

nlats = len(lats)
nlons = len(lons)


print DirI + FileI1
datain = xray.open_dataset(DirI + FileI1)

varlistin = []

if SA == 1:
        varlistin.append(['gridboxspanSA','totalprecipSA','uniquegridboxspanSA'])

newvar = "xmaxspeed_" + str(speedtspan) + "ts"
varlistin.append(newvar)

varlistin.extend(('gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax')) 

#if speedtspan == 0:	# have updated process script to not use this variable name now
#	newvar = "xmaxtravel"
#	newvarout = "xmaxspeed"
	
f4vars = np.array(varlistin)
f4varsin = np.array(varlistin)

f8vars = np.array(['eventid'])

nf4vars = len(f4vars)
nf8vars = len(f8vars)

nevents = len(datain.events)
chunksize = nevents

if splittype == "day":
	if filetimespan == '3hrly':
		splitvar = datain['timespan'].values /8.0	# Divide by 8 to get values in days instead of 3hours 
	else:
		exit("not set up for anything other than 3 hourly right now")
elif splittype == "speed":
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
        splitvar = direction * distancez/(timespan*3.0*60.0*60.0)

elif splittype == "maxspeed":	
	#if speedtspan == 0:
	#	splitvar = datain['xmaxtravel'].values
	#else:
	splitvar = datain['xmaxspeed_' + str(speedtspan) +  'ts'].values
else:
	exit("unacceptable splittype: " + str(splittype))

varforwrite = np.zeros([nf4vars+nf8vars,nevents],np.float)

for ibound in range(0,nbounds):
        if splittype == "maxspeed":
                fileadd = fileadd + "MaxSpeeds_" + str(speedtspan) + "ts_"
        elif splittype == "speed":
                fileadd = fileadd + "Speeds_"
        elif splittype == "day":
                fileadd = fileadd + "Sizes_"

	if tbound1[ibound] < 0:
		tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))
	else:
		tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))

        FileO = 'Test_Precip_' + fileadd + tboundtitle + unit + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

	print FileO
	# Create new file 
	ncfile = Dataset(DirO + FileO, 'w')
	ncfile.createDimension('events', None)

	Ofilevars = []
	for ivar in range(0,nf4vars):
		Ofilevars.append(ncfile.createVariable(f4vars[ivar],'f4',('events'),fill_value=-9999))

		units,desc = getunitsdesc(f4vars[ivar])
		setattr(Ofilevars[-1],'units',units)
                setattr(Ofilevars[-1],'description',desc)
	
	Ofilevars.append(ncfile.createVariable('eventid','f8',('events'),fill_value=-9999))

	#Extract events of startlen-endlen days
	# Do in chunks
        # Do in chunks
	readidxchunk = 0
        writeidxchunk = 0
        first = True
	for attempt in range(10):
		try:
			runchunk(nevents,chunksize,writeidxchunk)

		except MemoryError:
			chunksize = int(chunksize/4)
			print "failed on memory, quartering chunksize"
			continue
		break
				
	datain.close()
	ncfile.close()








