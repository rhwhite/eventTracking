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
test =0

def getindices(ibound,fromevent,toevent):
	print 'ibound, fromevent, toevent, inclusive: ', ibound, fromevent, toevent
	tstart = datain[ibound]['tstart'][fromevent].values
	tend = datain[ibound]['tstart'][toevent].values + datain[ibound]['timespan'][toevent].values
	print datain[ibound]['tstart'][toevent].values
	print datain[ibound]['timespan'][toevent].values

        eventid = datain[ibound]['eventid'].values
	giindices = eventid[fromevent:toevent+1] + minevent	# +1 as toevent is inclusive
	gitmin = int(max(0, tstart))
	gitmax = int(min(ntimes, tend+1))	# plus one as want to include this last time		

	return(giindices, gitmin, gitmax)

def geteventprecip(indicesin,tmin,tmax,TPrecip):
        # Here we include a multiplier to get from mm/hr (or mm/day) to mm/timestep (usually 3 hrs)
        precipchunk = mult * precipin.isel(time=slice(tmin, tmax)).sel(latitude=slice(minlat,maxlat)).values
        eventschunk = eventsin.isel(time=slice(tmin, tmax)).values
	for ibound in range(0,nbounds):
		data_mask = np.ma.array(eventschunk[:,0,:,:], mask=(np.in1d(eventschunk[:,0,:,:],indicesin[ibound]).reshape(eventschunk[:,0,:,:].shape)))
		TPrecip[ibound,:, :] += np.nansum(data_mask.mask * precipchunk,axis=0)    # Can't work out east vs west now! But this is kinda ok for speeds!

	return(TPrecip)

def runchunk(nchunks,fromevent,toevent,TPrecip):
	newtmin = tmin
	chunksizes = np.int((tmax-tmin)/nchunks)	# keep chunksizes as integer
        newtmax = newtmin + chunksizes

	for ichunk in range(nchunks-1):
		print ichunk, newtmin,newtmax
                TPrecip = geteventprecip(indicesL,newtmin,newtmax,TPrecip)

            	newtmin = newtmax
		newtmax = newtmin + chunksizes  
	
        # Last chunk
	print ichunk, newtmin, tmax
        Tprecip = geteventprecip(indicesL,newtmin,tmax,TPrecip)		
	return nchunks

def runmemory(nchunks,TPrecip):
	for attempt in range(10):
		try:
			try:
				del(eventschunk)
				del(precipchunk)
				del(data_mask)
			except NameError:
				pass
			nchunks = runchunk(nchunks,fromevent,toevent,TPrecip)
		except MemoryError:
			nchunks = nchunks * 2
			continue
		break
	else:
		exit("tried re-chunking 10 times but we just keep running out of memory!")

def writeall(ibound,allvars):
	for ivar in allvars:
		if ivar == 'ExmaxSpeed':
			writeFiles[ibound][ivar][iyear,:,:] = ExmaxSpeed[:,:]
		elif ivar == 'WxmaxSpeed':
			writeFiles[ibound][ivar][iyear,:,:] = WxmaxSpeed[:,:]
		elif ivar == 'ETotalPrecip':
			writeFiles[ibound][ivar][iyear,:,:] = ETotalPrecip[:,:]
		elif ivar == 'WTotalPrecip':
			writeFiles[ibound][ivar][iyear,:,:] = WTotalPrecip[:,:]
		elif ivar == 'STotalPrecip':
			writeFiles[ibound][ivar][iyear,:,:] = STotalPrecip[:,:]
		elif ivar == 'EDensity':
			writeFiles[ibound][ivar][iyear,:,:] = EDensity[:,:]
		elif ivar == 'WDensity':
			writeFiles[ibound][ivar][iyear,:,:] = WDensity[:,:]
		elif ivar == 'SDensity':
			writeFiles[ibound][ivar][iyear,:,:] = SDensity[:,:]
		elif ivar == 'TDensity':
			writeFiles[ibound][ivar][iyear,:,:] = np.nansum([EDensity[:,:],WDensity[:,:],SDensity[:,:]],axis=0)
		elif ivar == 'ESize': 
			writeFiles[ibound][ivar][iyear,:,:] = ESize[:,:]
		elif ivar == 'WSize': 
			writeFiles[ibound][ivar][iyear,:,:] = WSize[:,:]
		elif ivar == 'ETSpan':
			writeFiles[ibound][ivar][iyear,:,:] = ETSpan[:,:]
		elif ivar == 'WTSpan':
			writeFiles[ibound][ivar][iyear,:,:] = WTSpan[:,:]
                elif ivar == 'TPrecip':
                        writeFiles[ibound][ivar][iyear,:,:] = TPrecip[ibound,:,:]
                elif ivar == 'EPrecip':
                        writeFiles[ibound][ivar][iyear,:,:] = EPrecip[:,:]
                elif ivar == 'WPrecip':
                        writeFiles[ibound][ivar][iyear,:,:] = WPrecip[:,:]
                elif ivar == 'SPrecip':
                        writeFiles[ibound][ivar][iyear,:,:] = SPrecip[:,:]



nchunks = 2	# number of chunks to start with

mapping = 'center'	# options: genesis, center, end

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


def initialize(name, dims, atype):
    # globals used to create global object
    globals()[name] = np.zeros(dims, atype)


def resetvar(name):
    globals()[name][...] = 0

def delvar(name):
    del(name)


# Common to all tbounds
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

# Define which variables to print out, and which to initialize (all possible)
nonprecipvars = np.array(['ExmaxSpeed','WxmaxSpeed','EDensity','WDensity','SDensity','TDensity','ESize','WSize','SSize'])

precipvars = np.array(['TPrecip'])

allvars = np.concatenate((nonprecipvars,precipvars),axis=0)
extravars = np.array(['lon','lat','years'])
initvars = ['TPrecip', 'EPrecip', 'WPrecip', 'SPrecip', 'EDensity', 'WDensity','SDensity', 'TDensity', 'ExmaxSpeed', 'WxmaxSpeed', 'EDistance', 'WDistance', 'SDistance', 'ESize', 'WSize', 'SSize']

nvars = len(allvars)
nXvars = len(extravars)

nyears = anendyr - anstartyr

minev = np.zeros([nbounds,nyears], np.int)
maxev = np.zeros([nbounds,nyears], np.int)



# To run all ibounds at once, need to keep track of a. Input data, and b. Output data
writeFiles = []
datain = []
lons1 = []
lons2 = []


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

	# Create new datasets

	writeFiles.append(Dataset(DirO + FileO, 'w'))

	writeFiles[ibound].createDimension('years', nyears)
	writeFiles[ibound].createDimension('lon', nlons)
	writeFiles[ibound].createDimension('lat', nlats)

	for ivar in allvars:
		writeFiles[ibound].createVariable(ivar, 'f8', ('years', 'lat', 'lon'), fill_value=-9999)

                if ivar == "EPerc":
                        setattr(writeFiles[ibound][ivar], 'Extra Info', 'E Perc in percentage of Easterly; excluding stationary in whole calculation')
		elif ivar[0] == 'E': 
			setattr(writeFiles[ibound][ivar], 'Extra Info', 'Moving from the east, i.e. easterly; excluding stationary')
                elif ivar[0] == 'W':
                        setattr(writeFiles[ibound][ivar], 'Extra Info', 'Moving from the west, i.e. westerly; excluding stationary')

	for ivar in extravars:
		if ivar == "years":
			writeFiles[ibound].createVariable(ivar, 'f4', ('years'), fill_value=-9999)
		        for iyear in range(0,nyears):
            			writeFiles[ibound][ivar][iyear] = anstartyr + iyear
		elif ivar in ["Longitude","lon"]:	
			writeFiles[ibound].createVariable(ivar, 'f8', ('lon'), fill_value=-9999)
			writeFiles[ibound][ivar][:] = lons[:]
		elif ivar in ["Latitude","lat"]:
                        writeFiles[ibound].createVariable(ivar, 'f8', ('lat'), fill_value=-9999)
			writeFiles[ibound][ivar][:] = lats[:]
		else:
			print ivar
			sys.exit("not set up for this extra variable" + ivar)


	datain.append(XrayOpen(DirI + FileI1,decodetimes=False))

	#tspan = datain['timespan'].values
	tstart = datain[ibound]['tstart'].values

	nevents = tstart.shape[0]
	ievent = 1

	starttsteps = 0		# Can alter this to not start on first year
	if starttsteps == 0: print "******** starting at beginning of file, you should confirm this if you want accurate years ********"
	n = 0
	# find first timestep as defined above. 0 if beginning of file
	for n in range(0, nevents):
		curtime = tstart[n]
		if curtime >= starttsteps:
			minev[ibound,0] = n
			break

	# Loop through years and months to find the first and last timestep for each year
	totaldays = starttsteps     # keep track of days of simulation

	for iyear in range(0, nyears):
		totaldays = totaldays + (365.25 * 8)
		if minev[ibound,iyear] == nevents-1:	# if last minimum was nevents, then all subsequent mins and maxes will be this
			maxev[bound,iyear] = nevents-1
			if iyear < nyears-1:
				minev[ibound,iyear+1] = nevents-1
		else:
			for n in range(int(minev[ibound,iyear]), nevents):  # Start loop from beginning of last
			    curtime = tstart[n]
			    if curtime >= totaldays:        # End of this year
				if iyear < nyears-1:
					minev[ibound,iyear+1] = n

				if minev[ibound,iyear] == n:	# then no events in this last year
					maxev[ibound,iyear] = n
				else:
					maxev[ibound,iyear] = n-1	
				break
			else:
			    maxev[ibound,iyear] = nevents-1	# if didn't find it, then last event is max event
			    if iyear < nyears-1:
			        minev[ibound,iyear+1] = nevents-1

	# Get max and min lons for each event
	lons1.append(np.radians(lons[datain[ibound]['xcenterstart'].values.astype(int)]))
	lons2.append(np.radians(lons[datain[ibound]['xcenterend'].values.astype(int)]))

if np.any(maxev-minev < 0):
	print minev, maxev
	print maxev-minev
	exit("maxev-minev less than zero")

if test == 1:
    nloopyears = 1
else:
    nloopyears = nyears
tminchunk = 0
tmaxchunk = 0
precipchunk = []
eventschunk = []

print minev
print maxev

# Loop through years and months to count the number
for iyear in range(0, nloopyears):
	runningtotal = 0
	print 'year ',  iyear
	# Re-initialize all arrays
	# Create new datasets

	for ibound in range(0,nbounds):
		for var in allvars:
		    initialize(var, (nlats, nlons), np.float)
		for var in initvars:
		    initialize(var, (nlats, nlons), np.float)
		fromevent = minev[ibound][iyear]
		if test == 1:
		    toevent = min(minev[ibound][iyear] + 50,maxev[ibound][iyear])
		else:
		    toevent = maxev[ibound][iyear]

		counting = 0
		if mapping == 'genesis':
			ys = datain[ibound]['ycenterstart'].values
			xs = datain[ibound]['xcenterstart'].values
		elif mapping == 'center':
                        ys = datain[ibound]['ycentermean'].values
                        xs = datain[ibound]['xcentermean'].values
                elif mapping == 'end':
                        ys = datain[ibound]['ycenterend'].values
                        xs = datain[ibound]['xcenterend'].values
		
		gridboxspan = datain[ibound]['gridboxspan'].values
                xmaxspeed_4ts = datain[ibound]['xmaxspeed_4ts'].values


		if toevent > fromevent:
			for ievent in range(fromevent, toevent + 1):

				if (ievent % 5000000 == 0):
					print "ievent: " + str(ievent)  # keep track of program progress
			
				yidx = int(round(ys[ievent]))
				xidx = int(round(xs[ievent]))

				if lons2[ibound][ievent] > lons1[ibound][ievent]:  # If end lon > start lon it's moving to the east,  i.e. it's Westerly
					WDensity[yidx, xidx] += 1
					#WTotalPrecip[yidx, xidx] += totprecip[ievent]
					WSize[yidx, xidx] += gridboxspan[ievent]
					#WTSpan[yidx, xidx] += tspan[ievent]
					WxmaxSpeed[yidx,xidx] += xmaxspeed_4ts[ievent]
					#WzonalSpeed[yidx,xidx] += zonalspeed[ievent]
				elif lons2[ibound][ievent] < lons1[ibound][ievent]:
					EDensity[yidx, xidx] += 1
					#ETotalPrecip[yidx, xidx] += totprecip[ievent]
					ESize[yidx, xidx] += gridboxspan[ievent]
					#ETSpan[yidx, xidx] += tspan[ievent]
					ExmaxSpeed[yidx,xidx] += xmaxspeed_4ts[ievent]
					#EzonalSpeed[yidx,xidx] += zonalspeed[ievent]
				else:
					SDensity[yidx, xidx] += 1
					#ETotalPrecip[yidx, xidx] += totprecip[ievent]
					SSize[yidx, xidx] += gridboxspan[ievent]
					#STSpan[yidx, xidx] += tspan[ievent]

			writeall(ibound,nonprecipvars)

		for var in initvars:
		    delvar(var)

#		indices = eventid[fromevent:toevent] + minevent
#		indexend = eventid[toevent-1] + minevent
#		tmin = int(max(0, tstart[fromevent]-1))
#		tmax = int(min(ntimes, tend[toevent-1]+1))

	# Find minimum t and maximum t for all ibound
	for var in precipvars:
	    initialize(var, (nbounds,nlats, nlons), np.float)

        tminsA = np.zeros((nbounds),int)
        tmaxsA = np.zeros((nbounds),int)
	fromevents = np.zeros((nbounds),'int64')
	toevents = np.zeros((nbounds),'int64')
        indicesL = []
	for ibound in range(0,nbounds):

		fromevents[ibound] = minev[ibound][iyear]
		if test == 1:
		    toevents[ibound] = min(minev[ibound][iyear] + 50,maxev[ibound][iyear])
		else:
		    toevents[ibound] = maxev[ibound][iyear]
		print fromevents[ibound], toevents[ibound]
		indices, tminsA[ibound],tmaxsA[ibound] = getindices(ibound,fromevents[ibound],toevents[ibound])
		indicesL.append(np.array(indices))

	print tminsA
	print tmaxsA

	tmax = np.amax(tmaxsA)
	tmin = np.amin(tminsA)
	print tmin, tmax

	if not runinchunks:
		try:
			TPrecip = geteventprecip(indicesL,tmin,tmax,TPrecip)
			chunksize = toevent-fromevent
		except MemoryError:
			print "had memory issues, needed to chunk up!"
			# if necessary, do this in chunks
			runinchunks = True
			runmemory(nchunks,TPrecip)
	else:
		runmemory(nchunks,TPrecip)

	for ibound in range(0,nbounds):
		writeall(ibound,precipvars)

for ibound in range(0,nbounds):
        writeFiles[ibound].close()
	datain[ibound].close()



