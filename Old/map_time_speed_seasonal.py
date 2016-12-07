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
from rhwhitepackages.readwrite import XrayOpen


test = 0
mapprecip = 1
splittype = "day"	# "speed","maxspeed"
speedtspan = 0
# In m/s
if splittype == "speed" or splittype == "maxspeed":
	tbound = np.array([-30]) #[-1000,-30,-10,-6,-3,3,6,10,30])
	tbound2 = np.array([-6]) #[-30,-10,-6,-3,3,6,10,30,1000])
	unit = "ms"
elif splittype == "day":
        tbound = np.array([0,1,2,5,1])
        tbound2 = np.array([1,2,5,100,5])
	unit = "day"
nbounds = len(tbound)

Data = "ERAI"

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

if Data == "TRMM":
    anstartyr = 1998    # year for analysis start
    anendyr = 2015      # year for analysis end

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
    anstartyr = 1980    # normally 1980, test 2012 year for analysis start
    anendyr = 2015      # year for analysis end

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
        if speedtspan == 0:
                fileadd = ""
        else:   
                fileadd = str(speedtspan) + "ts_"

        tboundtitle = str(int(tbound[ibound])) + '-' + str(int(tbound2[ibound]))
	if splittype == "maxspeed":
		fileadd = fileadd + "MaxSpeeds_" + add
	elif splittype == "speed":
		fileadd = fileadd + "Speeds_" + add
	elif splittype == "day":
		fileadd = fileadd + "Sizes_" + add

        FileI1 = 'Precip_' + fileadd + tboundtitle + unit + '_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
	print FileI1
	if test == 0:
	    FileO = 'DenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(anstartyr) + '-' + str(anendyr) + '_' + Version + '.nc'
	else:
	    FileO = 'testDenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(anstartyr) + '-' + str(anendyr) + '_' + Version + '.nc'

	nyears = anendyr - anstartyr
	print nyears
	
	skipsteps = (anstartyr-startyr) * anntsteps
	nseas = 4
	# Timestep to start at: Assuming file starts at January, this is the first March
	if filetimespan == "3hrly":
		starttsteps = (31+28)*8 + skipsteps

	# Make sure we have the same number of years in each season.
	ntsteps = nyears * anntsteps

	mints = np.zeros([nyears, nseas], np.int)
	maxts = np.zeros([nyears, nseas], np.int)

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

	griddata = XrayOpen(DirP + Filegrid)
	SurfA = griddata['SurfaceArea']

	print DirEv + FileEv

	eventsdata = XrayOpen(DirEv + FileEv)
	eventsin = eventsdata['value']

	ntimes = eventsdata.time.size
	print DirP + FileP
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
	ncfile.createDimension('season', nseas)
	ncfile.createDimension('lon', nlons)
	ncfile.createDimension('lat', nlats)

	if splittype == "maxspeed":
		allvars = np.array(['ExmaxSpeed','WxmaxSpeed','TPrecip','EDensity','WDensity','SDensity','TDensity','ESize','WSize','SSize'])
	elif splittype == "speed":
                allvars = np.array(['TPrecip','EDensity','WDensity','SDensity','TDensity','ESize','WSize','SSize'])
	elif splittype == "day" or splittype == "speed":
                allvars = np.array(['TPrecip','EDensity','WDensity','SDensity','TDensity','ESize','WSize','SSize'])

	
	extravars = np.array(['Longitude','Latitude','years','Seasons'])
	initvars = ['ETotalPrecip', 'WTotalPrecip', 'STotalPrecip','TPrecip', 'EPrecip', 'WPrecip', 'SPrecip', 'EDensity', 'WDensity','SDensity', 'TDensity', 'ExmaxSpeed', 'WxmaxSpeed', 'EDistance', 'WDistance', 'SDistance', 'ESize', 'WSize', 'SSize', 'EUSize', 'WUSize', 'SUSize', 'ESizeSA', 'WSizeSA', 'SSizeSA', 'EUSizeSA', 'WUSizeSA', 'SUSizeSA', 'ETSpan', 'WTSpan', 'STSpan']

	nvars = len(allvars)
	nXvars = len(extravars)
	Ofilevars = []
	Oextravars = []

	Ocount = 0
	for ivar in range(0,nvars):
		Ofilevars.append(ncfile.createVariable(allvars[ivar], 'f8', ('years', 'season', 'lat', 'lon'), fill_value=-9999))

                if allvars[ivar] == "EPerc":
                        setattr(Ofilevars[Ocount], 'Extra Info', 'E Perc in percentage of Easterly; excluding stationary in whole calculation')

		elif allvars[ivar][0] == 'E': 
			setattr(Ofilevars[Ocount], 'Extra Info', 'Moving from the east, i.e. easterly; excluding stationary')
                elif allvars[ivar][0] == 'W':
                        setattr(Ofilevars[Ocount], 'Extra Info', 'Moving from the west, i.e. westerly; excluding stationary')

		Ocount += 1

	Oextracount = 0
	for ivar in range(0,nXvars):
		if extravars[ivar] == "Seasons":
			Oextravars.append(ncfile.createVariable(extravars[ivar], 'a3', ('season')))
                        for iseas in range(0,nseas):
                                Oextravars[Oextracount][iseas] = Seasons[iseas]
		elif extravars[ivar] == "years":
			Oextravars.append(ncfile.createVariable(extravars[ivar], 'f4', ('years'), fill_value=-9999))
		        for iyear in range(0,nyears):
            			Oextravars[Oextracount][iyear] = anstartyr + iyear
		elif extravars[ivar] == "Longitude":	
			Oextravars.append(ncfile.createVariable(extravars[ivar], 'f8', ('lon'), fill_value=-9999))
			Oextravars[Oextracount][:] = lons[:]
		elif extravars[ivar] == "Latitude":
                        Oextravars.append(ncfile.createVariable(extravars[ivar], 'f8', ('lat'), fill_value=-9999))
			Oextravars[Oextracount][:] = lats[:]
		else:
			sys.error("not set up for this extra variable" + extravars[ivar])


		Oextracount += 1


	print DirI + FileI1

	datain = XrayOpen(DirI + FileI1)

	eventid = datain['eventid'].values

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
	print nevents

	ievent = 1

	n = 0
	iyear = 0
	iseas = 0
	# find first timestep as defined above. 0 if beginning of file
	for n in range(0, nevents):
	    curtime = tstart[n]
	    if curtime >= starttsteps:
		mints[iyear, iseas] = n
		print('mints found', n)
		break

	# Loop through years and months to find the first and last timestep for each month
	totaldays = starttsteps     # keep track of days of simulation

	for iyear in range(0, nyears):
	    for iseas in range(0, nseas):
		totaldays = totaldays + seaststeps[iseas]

		for n in range(int(mints[iyear, iseas]), nevents):  # Start loop from beginning of last
		    curtime = tstart[n]
		    if curtime >= totaldays:        # End of this month
			if iseas < 3:
			    mints[iyear, iseas+1] = n
			else:
			    if iyear < nyears-1:
				mints[iyear+1, 0] = n

			maxts[iyear, iseas] = n  # switched to n instead of n-1 because python doesn't use the last index in a range
			break

	if maxts[nyears-1,3] == 0:
		maxts[nyears-1,3] = nevents 

	print maxts-mints
	print maxts

	if test == 1:
	    nloopyears = 1
	    nloopseas = 4
	else:
	    nloopyears = nyears
	    nloopseas = nseas
	tminchunk = 0
	tmaxchunk = 0
	precipchunk = []
	eventschunk = []

	for var in initvars:
	    initialize(var, nlats, nlons, np.float)

	#print maxts-mints
	lats1 = np.radians(lats[datain['ycenterstart']])
	lats2 = np.radians(lats[datain['ycenterend']])
	lons1 = np.radians(lons[datain['xcenterstart']])
	lons2 = np.radians(lons[datain['xcenterend']])

	# Loop through years and months to count the number
	for iyear in range(0, nloopyears):
	    runningtotal = 0
	    print 'year ',  iyear
	    for iseas in range(0, nloopseas):
		# Re-initialize all arrays
		# Create new datasets
		for var in allvars:
		    resetvar(var)
		for var in initvars:
		    resetvar(var)
		
		resetvar("WDensity")
		resetvar("SDensity")
		resetvar("EDensity")
		resetvar("TPrecip")

		print 'season ',  iseas
	#       print mints[iyear, iseas]
	#       print maxts[iyear, iseas]
		fromevent = mints[iyear,iseas]
		if test == 1:
		    toevent = min(mints[iyear, iseas] + 50,maxts[iyear,iseas])
		else:
		    toevent = maxts[iyear, iseas]

		counting = 0

		print fromevent, toevent
		for ievent in range(fromevent, toevent):

			if (ievent % 10000 == 0):
				print "ievent: " + str(ievent)  # keep track of program progress
			if mapping == 'genesis':
				yidx = int(round(ystart[ievent]))
				xidx = int(round(xstart[ievent]))
			elif mapping == 'centre':
				yidx = int(round(ycenter[ievent]))
				xidx = int(round(xcenter[ievent]))

			if lons2[ievent] > lons1[ievent]:  # If end lon > start lon it's moving to the east,  i.e. it's Westerly
				WDensity[yidx, xidx] += 1
				WTotalPrecip[yidx, xidx] += totprecip[ievent]
				WSize[yidx, xidx] += totsize[ievent]
				WTSpan[yidx, xidx] += tspan[ievent]
				if splittype == "maxspeed":
					WxmaxSpeed[yidx,xidx] += xmaxspeed[ievent]
				#WzonalSpeed[yidx,xidx] += zonalspeed[ievent]
			elif lons2[ievent] < lons1[ievent]:
				EDensity[yidx, xidx] += 1
				ETotalPrecip[yidx, xidx] += totprecip[ievent]
				ESize[yidx, xidx] += totsize[ievent]
				ETSpan[yidx, xidx] += tspan[ievent]
				if splittype == "maxspeed":
					ExmaxSpeed[yidx,xidx] += xmaxspeed[ievent]
				#EzonalSpeed[yidx,xidx] += zonalspeed[ievent]
			else:
				SDensity[yidx, xidx] += 1
				ETotalPrecip[yidx, xidx] += totprecip[ievent]
				SSize[yidx, xidx] += totsize[ievent]
				STSpan[yidx, xidx] += tspan[ievent]

		if mapprecip == 1:
			indices = eventid[fromevent:toevent] + minevent
			indexend = eventid[toevent] + minevent
			tmin = int(max(0, tstart[fromevent]-1))
			tmax = int(min(ntimes, tend[toevent]+1))

			eventschunk = eventsin.isel(time=slice(tmin, tmax)).values
			# Here we include a multiplier to get from mm/hr (or mm/day) to mm/timestep (usually 3 hrs)
			precipchunk = mult * precipin.isel(time=slice(tmin, tmax)).values

			ntimeshere = eventschunk.shape[0]

			runningtot = 0

			data_mask = np.ma.array(eventschunk[:,0,:,:], mask=(np.in1d(eventschunk[:,0,:,:],indices).reshape(eventschunk[:,0,:,:].shape)))
			#runningtotal = runningtotal + np.nansum(np.squeeze(data_mask_small.mask * precipin_small),axis=None)	# sum over all

			TPrecip[:, :] += np.sum(data_mask.mask * precipchunk,axis=0)	# Can't work out east vs west now! But this is kinda ok for speeds!
			#runningtot += np.nansum(np.squeeze(data_mask.mask * precipin_chunk))
			#if lons2[ievent] > lons1[ievent]:
			#	WPrecip[:,:] += np.sum(data_mask.mask * precipchunk,axis=0)
			#elif lons2[ievent] < lons1[ievent]:
                        #        EPrecip[:,:] += np.sum(data_mask.mask * precipchunk,axis=0)

		for ivar in range(0,nvars):
			if allvars[ivar] == 'ExmaxSpeed':
				Ofilevars[ivar][iyear,iseas,:,:] = ExmaxSpeed[:,:]
			elif allvars[ivar] == 'WxmaxSpeed':
                                Ofilevars[ivar][iyear,iseas,:,:] = WxmaxSpeed[:,:]
			elif allvars[ivar] == 'ETotalPrecip':
				Ofilevars[ivar][iyear,iseas,:,:] = ETotalPrecip[:,:]
			elif allvars[ivar] == 'WTotalPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = WTotalPrecip[:,:]
                        elif allvars[ivar] == 'STotalPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = STotalPrecip[:,:]
                        elif allvars[ivar] == 'EDensity':
                                Ofilevars[ivar][iyear,iseas,:,:] = EDensity[:,:]
                        elif allvars[ivar] == 'WDensity':
                                Ofilevars[ivar][iyear,iseas,:,:] = WDensity[:,:]
                        elif allvars[ivar] == 'SDensity':
                                Ofilevars[ivar][iyear,iseas,:,:] = SDensity[:,:]
                        elif allvars[ivar] == 'TDensity':
				Ofilevars[ivar][iyear,iseas,:,:] = np.nansum([EDensity,WDensity,SDensity],axis=0)
			elif allvars[ivar] == 'TPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = TPrecip[:,:]
                        elif allvars[ivar] == 'EPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = EPrecip[:,:]
                        elif allvars[ivar] == 'WPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = WPrecip[:,:]
                        elif allvars[ivar] == 'SPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = SPrecip[:,:]
                        elif allvars[ivar] == 'ESize':
                                Ofilevars[ivar][iyear,iseas,:,:] = ESize[:,:]
                        elif allvars[ivar] == 'WSize':
                                Ofilevars[ivar][iyear,iseas,:,:] = WSize[:,:]
                        elif allvars[ivar] == 'ETSpan':
                                Ofilevars[ivar][iyear,iseas,:,:] = ETSpan[:,:]
                        elif allvars[ivar] == 'WTSpan':
                                Ofilevars[ivar][iyear,iseas,:,:] = WTSpan[:,:]

	datain.close()
	ncfile.close()



