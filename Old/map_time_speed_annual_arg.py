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
#import datetime as dt
import re
import sys
import Ngl
import xray
import math
from rhwhitepackages.readwrite import XrayOpen
import argparse
import resource

rsrcV = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit starts as  :', soft
print 'Hard limit starts as :', hard
resource.setrlimit(rsrcV, (50000000000, hard)) #limit memory usage
#                          137438953472
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit changed to :', soft

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs=1,help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound1',type=int,nargs='+',help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=int,nargs="+",help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')

args = parser.parse_args()

print args.Data

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
startyr = args.startyr[0]
endyr = args.endyr[0]

nbounds = len(tbound1)

filetimespan = "3hrly"
test = 0

def getindices(fromevent,toevent):
	indices = eventid[fromevent:toevent] + minevent
	indexend = eventid[toevent-1] + minevent
	tmin = int(max(0, tstart[fromevent]-1))
	tmax = int(min(ntimes, tend[toevent-1]+1))

	return(indices, indexend, tmin, tmax)

def geteventprecip(tmin,tmax,TPrecip):
	eventschunk = eventsin.isel(time=slice(tmin, tmax)).values
	# Here we include a multiplier to get from mm/hr (or mm/day) to mm/timestep (usually 3 hrs)
	precipchunk = mult * precipin.isel(time=slice(tmin, tmax)).sel(latitude=slice(minlat,maxlat)).values

	data_mask = np.ma.array(eventschunk[:,0,:,:], mask=(np.in1d(eventschunk[:,0,:,:],indices).reshape(eventschunk[:,0,:,:].shape)))
	
	TPrecip[:, :] += np.sum(data_mask.mask * precipchunk,axis=0)    # Can't work out east vs west now! But this is kinda ok for speeds!
	return(TPrecip)

def runchunk(nchunks,fromevent,toevent):
	print 'nchunks', nchunks 
	newfromevent = fromevent
	chunksize = (toevent-fromevent)/nchunks

	for ichunk in range(nchunks-1):
        	newtoevent = newfromevent + chunksize
		indices, indexend, tmin, tmax = getindices(newfromevent,newtoevent)
		print ichunk, newfromevent, newtoevent, tmin,tmax
                geteventprecip(tmin,tmax,TPrecip)
            	newfromevent = newtoevent
		newtoevent = newfromevent + chunksize  
	
        # Last chunk
        indices, indexend, tmin, tmax = getindices(newfromevent,toevent)
	print newfromevent, toevent, tmin, tmax
        geteventprecip(tmin,tmax,TPrecip)		
	print "finished"
	return nchunks

def runmemory(nchunks):
	for attempt in range(10):
		try:
			try:
				del(eventschunk)
				del(precipchunk)
				del(data_mask)
			except NameError:
				pass
			nchunks = runchunk(nchunks,fromevent,toevent)
			print "done"
		except MemoryError:
			nchunks = nchunks * 2
			continue
		break
	else:
		exit("tried re-chunking 10 times but we just keep running out of memory!")


nchunks = 4	# number of chunks to start with

mapping = 'centre'

R = 6371000     # radius of Earth in m

if filetimespan == "3hrly":
        if Data == "ERAI" or Data == "TRMM" or Data == "TRMM_ERAIgd":
                # convert from mm/hour to mm/3 hours to get total rain over event for 3-hourly data
                mult = 3.0
		anntsteps = 365 * 8.0
        elif Data == "CESM":
                # convert from m/s to mm/3 hours to get total rain over event for 3-hourly data
                mult = 1000.0 * 60.0 * 60.0 * 3.0
		anntsteps = 365 * 8.0
        else:
                sys.error(Data + " not defined")

if test == 1:
        chunksize = 50
else:
        chunksize = 2000

Seasons = ["MAM", "JJA", "SON", "DJF"]
seaststeps = [(31+30+31)*8, (30+31+31)*8, (30+31+30)*8, (31+31+28)*8]

minevent = 100000

runinchunks = False

if Data == "TRMM":
    anstartyr = 1998    # year for analysis start
    anendyr = 2015      # year for analysis end

    startyr = 1998
    endyr = 2014
    if Version in ["Standard", "7thresh", "6th_from6", "5th_from48"]:
        DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
        FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'
        Filegrid = "SurfaceArea.nc"
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
        DirO = DirI
        DirEv = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/'
        FileEv = 'ts_TRMM1998-2014_final_4Dobjects.nc'
    elif Version in ['ERAIgd']:
        DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
        FileP = 'regrid2ERAI_TRMM_3B42_1998-2014.nc'
        Filegrid = "SurfaceArea.nc"
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/TRMM_output/' + Version + '/Precip/'
	FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
	DirO = DirI
        DirEv = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/TRMM_output/' + Version + '/'
        FileEv = 'ts_TRMM_ERAIgd_' + str(startyr) + '-' + str(endyr) + '_4Dobjects.nc'

    else:
        sys.exit('unexpected Version')
elif Data == "ERAI":
    anstartyr = 1980    # normally 1980, test 2012 year for analysis start
    anendyr = 2015      # year for analysis end

    Pfilestartyr = 1980
    Pfileendyr = 2015
    startyr = 1980
    endyr = 2015
    DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
    FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_preprocess.nc'

    DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
    DirO = DirI
    DirEv = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/'
    FileEv = 'ts_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
    Filegrid = "SurfaceArea.nc"

elif Data == "CESM":
    startyr = 1990   # Don't change - tied to file names!
    endyr = 2014

    anstartyr = 1990    # year for analysis start
    anendyr = 2014      # year for analysis end

    Pfilestartyr = 1979
    Pfileendyr = 2012
    DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
    FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '.nc'

    DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
    DirO = DirI
    DirEv = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/'
    FileEv = 'ts_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'

    Filegrid = "CESM_SurfaceArea.nc"
else:
    sys.error("unexpected datatype")


def initialize(name, numlats, numlons, atype):
    # globals used to create global object
    globals()[name] = np.zeros((numlats, numlons), atype)


def resetvar(name):
    globals()[name][...] = 0

def delvar(name):
    del(name)

writeFiles = []

for ibound in range(0,nbounds):
        if splittype == "maxspeed":
                fileadd = "MaxSpeeds_" + str(speedtspan) + "ts_"
        elif splittype == "speed":
                fileadd = "Speeds_"
        elif splittype == "day":
                fileadd = "Sizes_"

        if tbound1[ibound] < 0:
                tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))
        else:
                tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))

        FileI1 = 'Precip_' + fileadd + tboundtitle + unit + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

	if test == 0:
	    FileO = 'DenDirSpd_Map_Ann_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(anstartyr) + '-' + str(anendyr) + '_' + Version + '.nc'
	else:
	    FileO = 'testDenDirSpd_Map_Ann_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(anstartyr) + '-' + str(anendyr) + '_' + Version + '.nc'

	nyears = anendyr - anstartyr
	print nyears
	
	skipsteps = (anstartyr-startyr) * anntsteps
	# Timestep to start at: Assuming file starts at January

	# Make sure we have the same number of years in each season.
	ntsteps = nyears 

	mints = np.zeros([nyears], np.int)
	maxts = np.zeros([nyears], np.int)

	filetimespan = "3hrly"
	print DirP + FileP
	tempfile = XrayOpen(DirP + FileP)
	if Data == "CESM":
	    lons = tempfile['lon'].values
	    lats = tempfile['lat'].values
	else:
	    lons = tempfile['longitude'].values
	    lats = tempfile['latitude'].values

	nlons = len(lons)
	nlats = len(lats)
	minlat = np.amin(lats)
	maxlat = np.amax(lats)

	griddata = XrayOpen(DirP + Filegrid)
	SurfA = griddata['SurfaceArea']

	eventsdata = XrayOpen(DirEv + FileEv)
	eventsin = eventsdata['value']

	ntimes = eventsdata.time.size
	precipdata = XrayOpen(DirP + FileP)
	if Data == "TRMM" or Data == "TRMM_ERAIgd":
		precipin = precipdata['pcp']

	elif Data == "ERAI":
		precipin = precipdata['tpnew']
	elif Data == "CESM":
		precipin = precipdata['PRECT']

	ntimespre = len(precipdata['time'])

	# Create new datasets

	ncfile = Dataset(DirO + FileO, 'w')
	ncfile.createDimension('years', nyears)
	ncfile.createDimension('lon', nlons)
	ncfile.createDimension('lat', nlats)

	if splittype == "maxspeed":
		allvars = np.array(['ExmaxSpeed','WxmaxSpeed','TPrecip','EDensity','WDensity','SDensity','TDensity','ESize','WSize','SSize'])
	elif splittype == "speed":
                allvars = np.array(['TPrecip','EDensity','WDensity','SDensity','TDensity','ESize','WSize','SSize'])
	elif splittype == "day" or splittype == "speed":
                allvars = np.array(['TPrecip','EDensity','WDensity','SDensity','TDensity','ESize','WSize','SSize'])

	
	extravars = np.array(['lon','lat','years'])
	initvars = ['TPrecip', 'EPrecip', 'WPrecip', 'SPrecip', 'EDensity', 'WDensity','SDensity', 'TDensity', 'ExmaxSpeed', 'WxmaxSpeed', 'EDistance', 'WDistance', 'SDistance', 'ESize', 'WSize', 'SSize']

	nvars = len(allvars)
	nXvars = len(extravars)
	Ofilevars = []
	Oextravars = []

	Ocount = 0
	for ivar in range(0,nvars):
		Ofilevars.append(ncfile.createVariable(allvars[ivar], 'f8', ('years', 'lat', 'lon'), fill_value=-9999))

                if allvars[ivar] == "EPerc":
                        setattr(Ofilevars[Ocount], 'Extra Info', 'E Perc in percentage of Easterly; excluding stationary in whole calculation')

		elif allvars[ivar][0] == 'E': 
			setattr(Ofilevars[Ocount], 'Extra Info', 'Moving from the east, i.e. easterly; excluding stationary')
                elif allvars[ivar][0] == 'W':
                        setattr(Ofilevars[Ocount], 'Extra Info', 'Moving from the west, i.e. westerly; excluding stationary')

		Ocount += 1

	Oextracount = 0
	for ivar in range(0,nXvars):
		if extravars[ivar] == "years":
			Oextravars.append(ncfile.createVariable(extravars[ivar], 'f4', ('years'), fill_value=-9999))
		        for iyear in range(0,nyears):
            			Oextravars[Oextracount][iyear] = anstartyr + iyear
		elif extravars[ivar] in ["Longitude","lon"]:	
			Oextravars.append(ncfile.createVariable(extravars[ivar], 'f8', ('lon'), fill_value=-9999))
			Oextravars[Oextracount][:] = lons[:]
		elif extravars[ivar] in ["Latitude","lat"]:
                        Oextravars.append(ncfile.createVariable(extravars[ivar], 'f8', ('lat'), fill_value=-9999))
			Oextravars[Oextracount][:] = lats[:]
		else:
			print extravars[ivar]
			sys.exit("not set up for this extra variable" + extravars[ivar])


		Oextracount += 1

	datain = XrayOpen(DirI + FileI1,decodetimes=False)

	eventid = datain['eventid'].values

	#tspan = datain['timespan'].values.astype('timedelta64[h]').astype('int')
	tspan = datain['timespan'].values
	tstart = datain['tstart'].values
	tend = tstart + datain['timespan'].values

	xstart = datain['xcenterstart']
	xend = datain['xcenterend']
	ystart = datain['ycenterstart']
	yend = datain['ycenterend']

	xcenter = datain['xcentermean'].values
	ycenter = datain['ycentermean'].values

	totprecip = datain['totalprecip'].values
	totsize = datain['gridboxspan'].values
	totunisize = datain['uniquegridboxspan'].values

#	zonalspeed = datain['zonalspeed'].values
	if splittype == "maxspeed":
		if speedtspan == 0:
			xmaxspeed = datain['maxzonalspeed'].values
		elif speedtspan == 4:
			xmaxspeed = datain['xmaxspeed_4ts'].values
	

	nevents = tspan.shape[0]

	ievent = 1

	print "assuming we start at the beginning of a year"
	starttsteps = 0

	n = 0
	iyear = 0
	# find first timestep as defined above. 0 if beginning of file
	for n in range(0, nevents):
	    curtime = tstart[n]
	    if curtime >= starttsteps:
		mints[iyear] = n
		print('mints found', n)
		break

	# Loop through years and months to find the first and last timestep for each year
	totaldays = starttsteps     # keep track of days of simulation

	for iyear in range(0, nyears):
		totaldays = totaldays + (365.25 * 8)
		if mints[iyear] == nevents:
			maxts[iyear] = nevents
			if iyear < nyears-1:
				mints[iyear+1] = nevents
		else:
			for n in range(int(mints[iyear]), nevents):  # Start loop from beginning of last

			    curtime = tstart[n]
			    if curtime >= totaldays:        # End of this month
				if iyear < nyears-1:
					mints[iyear+1] = n

				maxts[iyear] = n 
				break
			else:
			    maxts[iyear] = nevents
			    if iyear < nyears-1:
			        mints[iyear+1] = nevents
	if maxts[nyears-1] == 0:
		maxts[nyears-1] = nevents 

	print mints
	print maxts
	print maxts-mints
	
	if test == 1:
	    nloopyears = 1
	else:
	    nloopyears = nyears
	tminchunk = 0
	tmaxchunk = 0
	precipchunk = []
	eventschunk = []

	for var in initvars:
	    initialize(var, nlats, nlons, np.float)

	#print maxts-mints
	print datain['xcenterstart'].values.astype(int)
	
	#lats1 = np.radians(lats[datain['ycenterstart']])
	#lats2 = np.radians(lats[datain['ycenterend']])
	lons1 = np.radians(lons[datain['xcenterstart'].values.astype(int)])
	lons2 = np.radians(lons[datain['xcenterend'].values.astype(int)])

	# Loop through years and months to count the number
	for iyear in range(0, nloopyears):
	    	runningtotal = 0
	    	print 'year ',  iyear
		# Re-initialize all arrays
		# Create new datasets
		
		for var in allvars:
		    initialize(var, nlats, nlons, np.float)
		for var in initvars:
		    initialize(var, nlats, nlons, np.float)
		
		resetvar("WDensity")
		resetvar("SDensity")
		resetvar("EDensity")
		resetvar("TPrecip")

		fromevent = mints[iyear]
		if test == 1:
		    toevent = min(mints[iyear] + 50,maxts[iyear])
		else:
		    toevent = maxts[iyear]

		counting = 0

		print fromevent, toevent
		if toevent > fromevent:
			for ievent in range(fromevent, toevent):

				if (ievent % 5000000 == 0):
					print "ievent: " + str(ievent)  # keep track of program progress
				if mapping == 'genesis':
					yidx = int(round(ystart[ievent]))
					xidx = int(round(xstart[ievent]))
				elif mapping == 'centre':
					yidx = int(round(ycenter[ievent]))
					xidx = int(round(xcenter[ievent]))

				if lons2[ievent] > lons1[ievent]:  # If end lon > start lon it's moving to the east,  i.e. it's Westerly
					WDensity[yidx, xidx] += 1
					#WTotalPrecip[yidx, xidx] += totprecip[ievent]
					WSize[yidx, xidx] += totsize[ievent]
					#WTSpan[yidx, xidx] += tspan[ievent]
					if splittype == "maxspeed":
						WxmaxSpeed[yidx,xidx] += xmaxspeed[ievent]
					#WzonalSpeed[yidx,xidx] += zonalspeed[ievent]
				elif lons2[ievent] < lons1[ievent]:
					EDensity[yidx, xidx] += 1
					#ETotalPrecip[yidx, xidx] += totprecip[ievent]
					ESize[yidx, xidx] += totsize[ievent]
					#ETSpan[yidx, xidx] += tspan[ievent]
					if splittype == "maxspeed":
						ExmaxSpeed[yidx,xidx] += xmaxspeed[ievent]
					#EzonalSpeed[yidx,xidx] += zonalspeed[ievent]
				else:
					SDensity[yidx, xidx] += 1
					#ETotalPrecip[yidx, xidx] += totprecip[ievent]
					SSize[yidx, xidx] += totsize[ievent]
					#STSpan[yidx, xidx] += tspan[ievent]

			for ivar in range(0,nvars):
				if allvars[ivar] == 'ExmaxSpeed':
					Ofilevars[ivar][iyear,:,:] = ExmaxSpeed[:,:]
				elif allvars[ivar] == 'WxmaxSpeed':
					Ofilevars[ivar][iyear,:,:] = WxmaxSpeed[:,:]
				elif allvars[ivar] == 'ETotalPrecip':
					Ofilevars[ivar][iyear,:,:] = ETotalPrecip[:,:]
				elif allvars[ivar] == 'WTotalPrecip':
					Ofilevars[ivar][iyear,:,:] = WTotalPrecip[:,:]
				elif allvars[ivar] == 'STotalPrecip':
					Ofilevars[ivar][iyear,:,:] = STotalPrecip[:,:]
				elif allvars[ivar] == 'EDensity':
					Ofilevars[ivar][iyear,:,:] = EDensity[:,:]
				elif allvars[ivar] == 'WDensity':
					Ofilevars[ivar][iyear,:,:] = WDensity[:,:]
				elif allvars[ivar] == 'SDensity':
					Ofilevars[ivar][iyear,:,:] = SDensity[:,:]
				elif allvars[ivar] == 'TDensity':
					Ofilevars[ivar][iyear,:,:] = np.nansum([EDensity,WDensity,SDensity],axis=0)
				elif allvars[ivar] == 'ESize':
					Ofilevars[ivar][iyear,:,:] = ESize[:,:]
				elif allvars[ivar] == 'WSize':
					Ofilevars[ivar][iyear,:,:] = WSize[:,:]
				elif allvars[ivar] == 'ETSpan':
					Ofilevars[ivar][iyear,:,:] = ETSpan[:,:]
				elif allvars[ivar] == 'WTSpan':
					Ofilevars[ivar][iyear,:,:] = WTSpan[:,:]


			for var in initvars:
			    delvar(var)

			indices, indexend, tmin, tmax = getindices(fromevent,toevent)
	#		indices = eventid[fromevent:toevent] + minevent
	#		indexend = eventid[toevent-1] + minevent
	#		tmin = int(max(0, tstart[fromevent]-1))
	#		tmax = int(min(ntimes, tend[toevent-1]+1))

			print "doing precip bit"

			if not runinchunks:
				try:
					geteventprecip(tmin,tmax,TPrecip)
					chunksize = toevent-fromevent
				except MemoryError:
					print "had memory issues, needed to chunk up!"
					# if necessary, do this in chunks
					runinchunks = True
					runmemory(nchunks)
			else:
				runmemory(nchunks)


			for ivar in range(0,nvars):
				if allvars[ivar] == 'TPrecip':
					Ofilevars[ivar][iyear,:,:] = TPrecip[:,:]
				elif allvars[ivar] == 'EPrecip':
					Ofilevars[ivar][iyear,:,:] = EPrecip[:,:]
				elif allvars[ivar] == 'WPrecip':
					Ofilevars[ivar][iyear,:,:] = WPrecip[:,:]
				elif allvars[ivar] == 'SPrecip':
					Ofilevars[ivar][iyear,:,:] = SPrecip[:,:]

	datain.close()
	ncfile.close()



