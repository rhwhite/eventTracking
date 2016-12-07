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

test =0

day1=2
day2=5

mapping = 'centre'
Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

anstartyr = 1998 #year for analysis start
anendyr = 2014 #year for analysis end
nyears = anendyr - anstartyr
nseas = 4
# Timestep to start at: Assuming file starts at January, this is the first March
starttsteps = (31+28)*8
# Make sure we have the same number of years in each season.
ntsteps = nyears * 365 * 8

print Version

R = 6371000	# radius of Earth in m

Data = "TRMM"

filetimespan = "3hrly"

if filetimespan == "3hrly":
        if Data == "ERAI" or Data == "TRMM" or Data == "TRMM_ERAIgd":
                # convert from mm/hour to mm/3 hours to get total rain over event for 3-hourly data
                mult = 3.0
        elif Data == "CESM":
                # convert from m/s to mm/3 hours to get total rain over event for 3-hourly data
                mult =  1000.0 * 60.0 * 60.0 * 3.0
        else:
                sys.error(Data + " not defined")

if test == 1:
        chunksize = 50
else:
        chunksize = 1000

Seasons = ["MAM","JJA","SON","DJF"]
seaststeps = [(31+30+31)*8,(30+31+31)*8,(30+31+30)*8,(31+31+28)*8]

mints = np.zeros([nyears,nseas],np.int)
maxts = np.zeros([nyears,nseas],np.int)

minevent = 100000

if Data == "TRMM":
	startyr = 1998
	endyr = 2014
	DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
	FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'

	if Version in ["Standard","7thresh","6th_from6","5th_from48"]:
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
		DirO = DirI
                DirEv = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/'
                FileEv = 'ts_TRMM1998-2014_final_4Dobjects.nc'
		FileI1 = 'Precip_Sizes_2-5day_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

	elif Version == '5th_nanto25':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'

		FileI1 = 'Precip_Sizes_2-5day_' + str(startyr) + '-' + str(endyr) + '_5th_nanto25.nc'

	elif Version == '5th_nantozero':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'

		FileI1 = 'Precip_Sizes_2-5day_' + str(startyr) + '-' + str(endyr) + '_5th_n2zero.nc'
	else:
		sys.exit('unexpected Version')
elif Data == "ERAI":
        DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
        FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_preprocess.nc'

	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	DirO = DirI
	FileI1 = 'Precip_Sizes_ERAI_2-5day_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

elif Data == "CESM":
	startyr = 1990 # Don't change - tied to file names!
	endyr = 2014

	Pfilestartyr = 1979
	Pfileendyr = 2012
        DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
        FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '.nc'

        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
	DirO = DirI
        FileI1 = 'Precip_Sizes_CESM_2-5day_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

if test == 0:
	FileO = 'DenDirSpd_Map_' + str(day1) + '-' + str(day2) +'Day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
else:
        FileO = 'testDenDirSpd_Map_' + str(day1) + '-' + str(day2) + 'Day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'


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
eventsin=eventsdata['value']

ntimes = eventsdata.time.size
print DirP + FileP
precipdata = xray.open_dataset(DirP + FileP)
if Data == "TRMM" or Data == "TRMM_ERAIgd":
        precipin = precipdata['pcp']
elif Data == "ERAI":
        precipin = precipdata['tpnew']
elif Data == "CESM":
        # Convert from m/s to mm/day
        precipin = precipdata['PRECT']

ntimespre = len(precipdata['time'])

#Create new datasets
#Total precip in event
ETotalPrecip = np.zeros((nyears,nseas,nlats,nlons),np.float)
WTotalPrecip = np.zeros((nyears,nseas,nlats,nlons),np.float)
STotalPrecip = np.zeros((nyears,nseas,nlats,nlons),np.float)

# Precip (per gridbox)
EPrecip = np.zeros((nyears,nseas,nlats,nlons),np.float)
WPrecip = np.zeros((nyears,nseas,nlats,nlons),np.float)
SPrecip = np.zeros((nyears,nseas,nlats,nlons),np.float)

#Number of events
EDensity = np.zeros((nyears,nseas,nlats,nlons),np.float)
WDensity = np.zeros((nyears,nseas,nlats,nlons),np.float)
SDensity = np.zeros((nyears,nseas,nlats,nlons),np.float)
TDensity = np.zeros((nyears,nseas,nlats,nlons),np.float)

#Average speed
ESpeed =  np.zeros((nyears,nseas,nlats,nlons),np.float)
WSpeed =  np.zeros((nyears,nseas,nlats,nlons),np.float)
SSpeed =  np.zeros((nyears,nseas,nlats,nlons),np.float)

# Average distance travelled
EDistance =  np.zeros((nyears,nseas,nlats,nlons),np.float)
WDistance =  np.zeros((nyears,nseas,nlats,nlons),np.float)
SDistance =  np.zeros((nyears,nseas,nlats,nlons),np.float)

# Average size (number of gridboxes (not unique))
ESize =  np.zeros((nyears,nseas,nlats,nlons),np.float)
WSize =  np.zeros((nyears,nseas,nlats,nlons),np.float)
SSize =  np.zeros((nyears,nseas,nlats,nlons),np.float)

# Average size (number of gridboxes (unique))
EUSize =  np.zeros((nyears,nseas,nlats,nlons),np.float)
WUSize =  np.zeros((nyears,nseas,nlats,nlons),np.float)
SUSize =  np.zeros((nyears,nseas,nlats,nlons),np.float)

# Average timespan
ETSpan =  np.zeros((nyears,nseas,nlats,nlons),np.float)
WTSpan =  np.zeros((nyears,nseas,nlats,nlons),np.float)
STSpan =  np.zeros((nyears,nseas,nlats,nlons),np.float)

# Ratio of E to W
EPerc = np.zeros((nyears,nseas,nlats,nlons),np.int)

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('years', nyears)
ncfile.createDimension('season',nseas)
ncfile.createDimension('size', 4)
ncfile.createDimension('lon', nlons)
ncfile.createDimension('lat', nlats)


OutETotalPrecip = ncfile.createVariable('ETotalPrecip','f4',('years','season','lat','lon'),fill_value=-9999)
OutWTotalPrecip = ncfile.createVariable('WTotalPrecip','f4',('years','season','lat','lon'),fill_value=-9999)
OutSTotalPrecip = ncfile.createVariable('STotalPrecip','f4',('years','season','lat','lon'),fill_value=-9999)

OutEPrecip = ncfile.createVariable('EPrecip','f4',('years','season','lat','lon'),fill_value=-9999)
OutWPrecip = ncfile.createVariable('WPrecip','f4',('years','season','lat','lon'),fill_value=-9999)
OutSPrecip = ncfile.createVariable('SPrecip','f4',('years','season','lat','lon'),fill_value=-9999)


OutEDensity = ncfile.createVariable('EDensity','f4',('years','season','lat','lon'),fill_value=-9999)
OutWDensity = ncfile.createVariable('WDensity','f4',('years','season','lat','lon'),fill_value=-9999)
OutSDensity = ncfile.createVariable('SDensity','f4',('years','season','lat','lon'),fill_value=-9999)
OutTDensity = ncfile.createVariable('TDensity','f4',('years','season','lat','lon'),fill_value=-9999)

OutESpeed = ncfile.createVariable('ESpeed','f4',('years','season','lat','lon'),fill_value=-9999)
OutWSpeed = ncfile.createVariable('WSpeed','f4',('years','season','lat','lon'),fill_value=-9999)
OutSSpeed = ncfile.createVariable('SSpeed','f4',('years','season','lat','lon'),fill_value=-9999)

OutEDistance = ncfile.createVariable('EDistance','f4',('years','season','lat','lon'),fill_value=-9999)
OutWDistance = ncfile.createVariable('WDistance','f4',('years','season','lat','lon'),fill_value=-9999)
OutSDistance = ncfile.createVariable('SDistance','f4',('years','season','lat','lon'),fill_value=-9999)

OutESize = ncfile.createVariable('ESize','f4',('years','season','lat','lon'),fill_value=-9999)
OutWSize = ncfile.createVariable('WSize','f4',('years','season','lat','lon'),fill_value=-9999)
OutSSize = ncfile.createVariable('SSize','f4',('years','season','lat','lon'),fill_value=-9999)

OutEUSize = ncfile.createVariable('EUniSize','f4',('years','season','lat','lon'),fill_value=-9999)
OutWUSize = ncfile.createVariable('WUniSize','f4',('years','season','lat','lon'),fill_value=-9999)
OutSUSize = ncfile.createVariable('SUniSize','f4',('years','season','lat','lon'),fill_value=-9999)

OutETSpan = ncfile.createVariable('ETSpan','f4',('years','season','lat','lon'),fill_value=-9999)
OutWTSpan = ncfile.createVariable('WTSpan','f4',('years','season','lat','lon'),fill_value=-9999)
OutSTSpan = ncfile.createVariable('STSpan','f4',('years','season','lat','lon'),fill_value=-9999)

OutEPerc = ncfile.createVariable('EPerc','f4',('years','season','lat','lon'),fill_value=-9999)

OutLongitude = ncfile.createVariable('Longitude','f4',('lon'),fill_value=-9999)
OutLatitude = ncfile.createVariable('Latitude','f4',('lat'),fill_value=-9999)
OutYears = ncfile.createVariable('years','f4',('years'),fill_value=-9999)
OutSeasons = ncfile.createVariable('Seasons','a3',('season'))


setattr(OutESpeed,'Extra Info','Speed calculated in m/s')
setattr(OutEDistance,'Extra Info','Distance calculated in m')
setattr(OutEPrecip,'Extra Info','EPrecip is easterly objects, i.e. from the East')
setattr(OutWPrecip,'Extra Info','WPrecip is westerly objects, i.e. from the West')
setattr(OutEPerc,'Extra Info','E Perc in percentage of Easterly; westerly includes stationary')


OutLongitude[:] = lons
OutLatitude[:] = lats

print DirI + FileI1

datain = xray.open_dataset(DirI + FileI1)

eventid = datain['eventid'].values

tspan=datain['timespan'].values
tstart = datain['tstart'].values
tend = tstart + datain['timespan'].values

xstart = datain['xcenterstart']
xend = datain['xcenterend']
ystart = datain['ycenterstart']
yend = datain['ycenterend']

xcenter = datain['xcentermean'].values
ycenter = datain['ycentermean'].values

totprecip = datain['totalprecipSA'].values
totsize = datain['gridboxspanSA'].values
totunisize = datain['uniquegridboxspanSA'].values


nevents = tspan.shape[0]
print nevents

ievent = 1

curtime = tstart[starttsteps]
n = 0
iyear = 0
iseas = 0
# find first timestep as defined above. 0 if beginning of file
for n in range(0,nevents):
	curtime = tstart[n]
	if curtime >= starttsteps:
		mints[iyear,iseas] = n
		print('mints found',n)
		break

# Loop through years and months to find the first and last timestep for each month
totaldays = starttsteps		# keep track of days of simulation

Seasons = ["MAM","JJA","SON","DJF"]
OutSeasons[0] = Seasons[0]
OutSeasons[1] = Seasons[1]
OutSeasons[2] = Seasons[2]
OutSeasons[3] = Seasons[3]


for iyear in range(0,nyears):
        OutYears[iyear] = anstartyr + iyear

	for iseas in range(0,nseas):
		totaldays = totaldays + seaststeps[iseas]
		for n in range(int(mints[iyear,iseas]),nevents):	# Start loop from beginning of last
			curtime = tstart[n]
			if curtime >= totaldays:		# End of this month
				if iseas < 3:
					mints[iyear,iseas+1] = n
				else:
					if iyear < nyears-1:
						mints[iyear+1,0] = n					

				maxts[iyear,iseas] = n  #switched to n instead of n-1 because python doesn't use the last index in a range
				break

# Calculate distance using the Haversine formula, so we can do this using numpy
lats1 = np.radians(lats[ystart[:]])
lats2 = np.radians(lats[yend[:]])
lons1 = np.radians(lons[xstart[:]])
lons2 = np.radians(lons[xend[:]]) 

a = (np.power(np.sin((lats2-lats1)/2),2) + np.cos(lats2) * np.cos(lats1) * np.power(np.sin((lons2-lons1)/2),2))
c = 2.0 * np.arctan(np.sqrt(a),np.sqrt(1-a))
distance = R * c
pi = 3.14

angle = np.arctan2((lats2-lats1),(lons2-lons1))

#print angle[4:6]
#print distance[4:6]
#print tspan[4:6]

if test == 1:
	nloopyears = 1
	nloopseas = 1
else:
	nloopyears = nyears
	nloopseas = nseas
tminchunk = 0
tmaxchunk = 0

# Loop through years and months to count the number
for iyear in range(0,nloopyears):
	print 'year ', iyear
	for iseas in range(0,nloopseas):
		print 'season ', iseas
		print mints[iyear,iseas]
		print maxts[iyear,iseas]
		if test == 1:
			toevent = mints[iyear,iseas] + 10
		else:
			toevent = maxts[iyear,iseas]	

		for ievent in range(mints[iyear,iseas],toevent):
			
			if (ievent % 1000 == 0):
				print "ievent: " + str(ievent)	# keep track of program progress	
			if mapping == 'genesis':
				yidx = int(round(ystart[ievent]))
				xidx = int(round(xstart[ievent]))
			elif mapping == 'centre':
                                yidx = int(round(ycenter[ievent]))
                                xidx = int(round(xcenter[ievent]))
			
	
			if lons2[ievent] == lons1[ievent]:
			        SDensity[iyear,iseas,yidx,xidx] += 1
                                STotalPrecip[iyear,iseas,yidx,xidx] += totprecip[ievent]
                                SSpeed[iyear,iseas,yidx,xidx] += distance[ievent]/(tspan[ievent]*3.0*60.0*60.0)
                                SDistance[iyear,iseas,yidx,xidx] += distance[ievent]
                                SSize[iyear,iseas,yidx,xidx] += totsize[ievent]
                                SUSize[iyear,iseas,yidx,xidx] += totunisize[ievent]
                                STSpan[iyear,iseas,yidx,xidx] += tspan[ievent]

			elif lons2[ievent] > lons1[ievent]:	# If end lon > start lon it's moving to the east, i.e. it's Westerly

                                WDensity[iyear,iseas,yidx,xidx] += 1
				WTotalPrecip[iyear,iseas,yidx,xidx] += totprecip[ievent]
				WSpeed[iyear,iseas,yidx,xidx] += distance[ievent]/(tspan[ievent]*3.0*60.0*60.0) 
				WDistance[iyear,iseas,yidx,xidx] += distance[ievent]
				WSize[iyear,iseas,yidx,xidx] += totsize[ievent]
				WUSize[iyear,iseas,yidx,xidx] += totunisize[ievent]
				WTSpan[iyear,iseas,yidx,xidx] += tspan[ievent]
			else:
                                EDensity[iyear,iseas,yidx,xidx] += 1
                                ETotalPrecip[iyear,iseas,yidx,xidx] += totprecip[ievent]
                                ESpeed[iyear,iseas,yidx,xidx] += distance[ievent]/(tspan[ievent]*3.0*60.0*60.0)
                                EDistance[iyear,iseas,yidx,xidx] += distance[ievent]
                                ESize[iyear,iseas,yidx,xidx] += totsize[ievent]
                                EUSize[iyear,iseas,yidx,xidx] += totunisize[ievent]
                                ETSpan[iyear,iseas,yidx,xidx] += tspan[ievent]



		        index = eventid[ievent] + minevent
		        tmin = int(max(0,tstart[ievent]-1))
        		tmax = int(min(ntimes,tend[ievent]+1))
        		if (tmax > tmaxchunk):
                		#print "tmax is " + str(tmax) + " whilst tmaxchunk is " + str(tmaxchunk)
                		tminchunk = max(tmin - 10,0)
                		tmaxchunk = min(tmax + chunksize + 1,ntimes)
                		#print "tminchunk is now " + str(tminchunk) + "whilst tmin is " + str(tmin)
                		#print "tmaxchunk is now " + str(tmaxchunk) + "whilst tmax is " + str(tmax)
                		eventschunk = eventsin.isel(time=slice(tminchunk,tmaxchunk)).values
                		# Here we include a multiplier to get from mm/hr (or mm/day) to mm/timestep (usually 3 hrs)
                		precipchunk = mult * precipin.isel(time=slice(tminchunk,tmaxchunk)).values
				ntimeshere = eventschunk.shape[0]
        		tminsel = max((tmin-tminchunk)-1,0)
        		tmaxsel = min((tmax-tminchunk)+2,ntimeshere)
			
			for it in range(tminsel,tmaxsel):
			        eventsin_small = eventschunk[it,:,:]
        			data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == index))
        			precipin_small = precipchunk[it,:,:]
			
				if lons2[ievent] == lons1[ievent]:
					SPrecip[iyear,iseas,:,:] += np.squeeze(data_mask_small.mask * precipin_small)
				elif lons2[ievent] > lons1[ievent]:
                                        WPrecip[iyear,iseas,:,:] += np.squeeze(data_mask_small.mask * precipin_small)
				else:
					EPrecip[iyear,iseas,:,:] += np.squeeze(data_mask_small.mask * precipin_small)


TDensity = WDensity + EDensity + SDensity
EPerc = EDensity * 100.0 / (TDensity)

WSpeed = WSpeed / WDensity
ESpeed = ESpeed / EDensity
SSpeed = SSpeed / SDensity

WDistance = WDistance / WDensity
EDistance = EDistance / EDensity
SDistance = SDistance / SDensity

WSize = WSize / WDensity
ESize = ESize / EDensity
SSize = SSize / SDensity

WUSize = WUSize / WDensity
EUSize = EUSize / EDensity
SUSize = SUSize / SDensity

WTSpan = WTSpan / WDensity
ETSpan = ETSpan / EDensity
STSpan = STSpan / SDensity


WSpeed[np.where(WDensity == 0)] = np.nan
WDistance[np.where(WDensity == 0)] = np.nan
WSize[np.where(WDensity == 0)] = np.nan
WUSize[np.where(WDensity == 0)] = np.nan
WTSpan[np.where(WDensity == 0)] = np.nan

ESpeed[np.where(EDensity == 0)] = np.nan
EDistance[np.where(EDensity == 0)] = np.nan
ESize[np.where(EDensity == 0)] = np.nan
EUSize[np.where(EDensity == 0)] = np.nan
ETSpan[np.where(EDensity == 0)] = np.nan



OutETotalPrecip[...]  = ETotalPrecip
OutWTotalPrecip[...]  = WTotalPrecip
OutSTotalPrecip[...]  = STotalPrecip

OutEPrecip[...]  = EPrecip
OutWPrecip[...]  = WPrecip
OutSPrecip[...]  = SPrecip

OutEDensity[...]  = EDensity
OutWDensity[...]  = WDensity
OutSDensity[...]  = SDensity
OutTDensity[...]  = TDensity

OutESpeed[...]  = ESpeed
OutWSpeed[...]  = WSpeed
OutSSpeed[...]  = SSpeed

OutEDistance[...]  = EDistance
OutWDistance[...]  = WDistance
OutSDistance[...]  = SDistance

OutESize[...]  = ESize
OutWSize[...]  = WSize
OutSSize[...]  = SSize

OutEUSize[...]  = EUSize
OutWUSize[...]  = WUSize
OutSUSize[...]  = SUSize

OutETSpan[...] = ETSpan
OutWTSpan[...] = WTSpan
OutSTSpan[...] = STSpan

OutEPerc[...] = EPerc

datain.close()
ncfile.close()


