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
maxspeed = 1

# In m/s
tbound = np.array([1000,30,10,6,3])
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
        chunksize = 4000

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


for ibound in range(0,nbounds-1):
	print  str(int(tbound[ibound])) + '-' + str(int(tbound[ibound+1])) + 'ms_'
	if maxspeed ==1:
		fileadd = "MaxSpeeds_" + add
	else:
		fileadd = "Speed_" + add
        FileI1 = 'Precip_' + fileadd + str(int(tbound[ibound])) + '-' + str(int(tbound[ibound+1])) + 'ms_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

	if test == 0:
	    FileO = 'DenDirSpd_Map_' + fileadd + str(int(tbound[ibound])) + '-' + str(int(tbound[ibound+1])) + 'ms_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
	else:
	    FileO = 'testDenDirSpd_Map_' + fileadd + str(int(tbound[ibound])) + '-' + str(int(tbound[ibound+1])) + 'ms_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

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

	griddata = xray.open_dataset(DirP + Filegrid)
	SurfA = griddata['SurfaceArea']


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

	ncfile = Dataset(DirO + FileO, 'w')
	ncfile.createDimension('years', nyears)
	ncfile.createDimension('season', nseas)
	ncfile.createDimension('size', 4)
	ncfile.createDimension('lon', nlons)
	ncfile.createDimension('lat', nlats)

	allvars = np.array(['ExmaxSpeed','WxmaxSpeed','EPrecip','WPrecip','SPrecip','EDensity','WDensity','SDensity','ESize','WSize','SSize'])
	extravars = np.array(['Longitude','Latitude','years','Seasons'])
#	allvars = ['ETotalPrecip', 'WTotalPrecip', 'STotalPrecip', 'EPrecip', 'WPrecip', 'SPrecip', 'EDensity', 'WDensity', 'TDensity', 'EZSpeed', 'WZSpeed', 'SZSpeed', 'EDistance', 'WDistance', 'SDistance', 'ESize', 'WSize', 'SSize', 'EUSize', 'WUSize', 'SUSize', 'ESizeSA', 'WSizeSA', 'SSizeSA', 'EUSizeSA', 'WUSizeSA', 'SUSizeSA', 'ETSpan', 'WTSpan', 'STSpan']

	nvars = len(allvars)
	print nvars
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
			Oextravars.append(ncfile.createVariable(extravars[ivar], 'f8', ('years', 'season', 'lat', 'lon'), fill_value=-9999))
			Oextravars[Oextracount] = lons
		elif extravars[ivar] == "Latitude":
                        Oextravars.append(ncfile.createVariable(extravars[ivar], 'f8', ('years', 'season', 'lat', 'lon'), fill_value=-9999))
			Oextravars[Oextracount] = lats
		else:
			sys.error("not set up for this extra variable" + extravars[ivar])


		Oextracount += 1


	print DirI + FileI1

	datain = xray.open_dataset(DirI + FileI1)

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
	xmaxspeed = datain['maxzonalspeed'].values

	nevents = tspan.shape[0]
	print nevents

	ievent = 1

	curtime = tstart[starttsteps]
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

	if test == 1:
	    nloopyears = 1
	    nloopseas = 1
	else:
	    nloopyears = nyears
	    nloopseas = nseas
	tminchunk = 0
	tmaxchunk = 0
	precipchunk = []
	eventschunk = []


	def initialize(name, numlats, numlons, atype):
	    # globals used to create global object
	    globals()[name] = np.zeros((numlats, numlons), atype)


	def resetvar(name):
	    globals()[name][...] = 0


	for var in allvars:
	    initialize(var, nlats, nlons, np.float)

	print maxts-mints

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

		print 'season ',  iseas
	#       print mints[iyear, iseas]
	#       print maxts[iyear, iseas]
		if test == 1:
		    toevent = mints[iyear, iseas] + 10
		else:
		    toevent = maxts[iyear, iseas]

		counting = 0

		for ievent in range(mints[iyear, iseas], toevent):

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
			#WTotalPrecip[yidx, xidx] += totprecip[ievent]
			WSize[yidx, xidx] += totsize[ievent]
			#WTSpan[yidx, xidx] += tspan[ievent]
			WxmaxSpeed[yidx,xidx] += xmaxspeed[ievent]
			#WzonalSpeed[yidx,xidx] += zonalspeed[ievent]
		    elif lons2[ievent] < lons1[ievent]:
			EDensity[yidx, xidx] += 1
			#ETotalPrecip[yidx, xidx] += totprecip[ievent]
			ESize[yidx, xidx] += totsize[ievent]
			#ETSpan[yidx, xidx] += tspan[ievent]
                        ExmaxSpeed[yidx,xidx] += xmaxspeed[ievent]
                        #EzonalSpeed[yidx,xidx] += zonalspeed[ievent]
		    else:
                        SDensity[yidx, xidx] += 1
                        #ETotalPrecip[yidx, xidx] += totprecip[ievent]
                        #ESize[yidx, xidx] += totsize[ievent]
                        #ETSpan[yidx, xidx] += tspan[ievent]
                        #ExmaxSpeed[yidx,xidx] += xmaxspeed[ievent]

		    if mapprecip == 1:

			index = eventid[ievent] + minevent
			tmin = int(max(0, tstart[ievent]-1))
			tmax = int(min(ntimes, tend[ievent]+1))

			if (tmax > tmaxchunk):
				print 'rereading'
				del(precipchunk)
				del(eventschunk)
				# print "tmax is " + str(tmax) + " whilst tmaxchunk is " + str(tmaxchunk)
				tminchunk = max(tmin - 10, 0)
				tmaxchunk = min(tmax + chunksize + 1, ntimes)
				# print "tminchunk is now " + str(tminchunk) + "whilst tmin is " + str(tmin)
				# print "tmaxchunk is now " + str(tmaxchunk) + "whilst tmax is " + str(tmax)
				eventschunk = eventsin.isel(time=slice(tminchunk, tmaxchunk)).values
				# Here we include a multiplier to get from mm/hr (or mm/day) to mm/timestep (usually 3 hrs)
				precipchunk = mult * precipin.isel(time=slice(tminchunk, tmaxchunk)).values

				ntimeshere = eventschunk.shape[0]

				# print eventschunk.shape
				# print precipchunk.shape
				# print np.nansum(precipchunk)

			tminsel = max((tmin-tminchunk)-1, 0)
			tmaxsel = min((tmax-tminchunk)+2, ntimeshere)
			runningtot = 0
	
			eventsin_small = eventschunk[tminsel:tmaxsel, 0, :, :]
			
			data_mask_small = np.ma.array(eventsin_small, mask=(eventsin_small == index))
			
			precipin_small = precipchunk[tminsel:tmaxsel, :, :]
			runningtotal = runningtotal + np.nansum(np.squeeze(data_mask_small.mask * precipin_small),axis=None)	# sum over all

			if lons2[ievent] > lons1[ievent]:
			    WPrecip[:, :] += np.sum(data_mask_small.mask * precipin_small,axis=0)
			elif lons2[ievent] < lons1[ievent]:
			    EPrecip[:, :] += np.sum(data_mask_small.mask * precipin_small,axis=0)
		        else:
			    SPrecip[:,:] += np.sum(data_mask_small.mask * precipin_small)
			runningtot += np.nansum(np.squeeze(data_mask_small.mask * precipin_small))

			if abs(runningtot - totprecip[ievent]) > runningtot * 0.01:
				print 'ERROR!!'
				print runningtot
				print totprecip[ievent]
				tmid = (tminsel + tmaxsel) * 0.5
				eventtemp = eventschunk[tminsel:tmaxsel,0,:,:]
				print np.amin(eventtemp[np.nonzero(eventtemp)]), np.amax(eventtemp)
				print np.where(eventtemp == index)	
				np.set_printoptions(threshold=np.nan)
				
				print 'index is ', index
				print 'ievent is ', ievent
				print tmid
				print eventschunk[:,0,:,:].shape
				print lons2[ievent], lons1[ievent], lats2[ievent], lats1[ievent]
				print lats[ystart[ievent]], lats[yend[ievent]], lons[xstart[ievent]],lons[xend[ievent]]
				print tminsel, tmaxsel
				print np.nansum(precipin[tminsel:tmaxsel, :, :])
				print tspan[ievent], tstart[ievent], tend[ievent]
				exit()
	#       print np.nansum(ETotalPrecip),  np.nansum(EPrecip*SurfA),  np.nansum(EPrecip)
	#       print np.nansum(WTotalPrecip),  np.nansum(WPrecip*SurfA),  np.nansum(WPrecip)
	#       print np.nansum(STotalPrecip),  np.nansum(SPrecip*SurfA),  np.nansum(SPrecip)

		for ivar in range(0,nvars):
			if allvars[ivar] == 'ExmaxSpeed':
				Ofilevars[ivar][iyear,iseas,:,:] = ExmaxSpeed
			elif allvars[ivar] == 'WxmaxSpeed':
                                Ofilevars[ivar][iyear,iseas,:,:] = WxmaxSpeed
			elif allvars[ivar] == 'ETotalPrecip':
				Ofilevars[ivar][iyear,iseas,:,:] = ETotalPrecip
			elif allvars[ivar] == 'WTotalPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = WTotalPrecip
                        elif allvars[ivar] == 'STotalPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = STotalPrecip
                        elif allvars[ivar] == 'EDensity':
                                Ofilevars[ivar][iyear,iseas,:,:] = EDensity
                        elif allvars[ivar] == 'WDensity':
                                Ofilevars[ivar][iyear,iseas,:,:] = WDensity
                        elif allvars[ivar] == 'SDensity':
                                Ofilevars[ivar][iyear,iseas,:,:] = SDensity
                        elif allvars[ivar] == 'EPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = EPrecip
                        elif allvars[ivar] == 'WPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = WPrecip
                        elif allvars[ivar] == 'SPrecip':
                                Ofilevars[ivar][iyear,iseas,:,:] = SPrecip
                        elif allvars[ivar] == 'ESize':
                                Ofilevars[ivar][iyear,iseas,:,:] = ESize
                        elif allvars[ivar] == 'WSize':
                                Ofilevars[ivar][iyear,iseas,:,:] = WSize
                        elif allvars[ivar] == 'ETSpan':
                                Ofilevars[ivar][iyear,iseas,:,:] = ETSpan
                        elif allvars[ivar] == 'WTSpan':
                                Ofilevars[ivar][iyear,iseas,:,:] = WTSpan

	datain.close()
	ncfile.close()



