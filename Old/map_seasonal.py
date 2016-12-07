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


Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

print Version

R = 6371000	# radius of Earth in m

# Testing for haversine and angle formula
# Angle is defined relative to lat-lon, not spherical so that angle is 0 if travelling along a latitude circle
#lat1 = np.array([46.0,45.0,45.0,46.0])
#lat2 = np.array([45.0,45.0,46.0,45.0])
#lon1 = np.array([0.0,1.0,1.0,1.0])
#lon2 = np.array([0.0,0.0,0.0,0.0])
#print(lat1,lat2,lon1,lon2)
# convert to radians
#np.radians(lat1,lat1)
#np.radians(lat2,lat2)
#np.radians(lon1,lon1)
#np.radians(lon2,lon2)
#
#a = (np.power(np.sin((lat2-lat1)/2),2) + np.cos(lat2) * np.cos(lat1) * np.power(np.sin((lon2-lon1)/2),2))
#c = 2.0 * np.arcsin(np.minimum(1,np.sqrt(a)))
#d = R * c
#print d / 1000.0
#
#angle2 = np.arctan2((lat2-lat1),(lon2-lon1))
#
#quit()

Data = "CESM"

filetimespan = "3hrly"

startyr = 1990 # Don't change - tied to file names!
endyr = 2014

Pfilestartyr = 1979
Pfileendyr = 2012

anstartyr = 1990 #year for analysis start
anendyr = 2011 #year for analysis end
nyears = anendyr - anstartyr + 1
nmonths =12
nmonthsrun = 12
mints = np.zeros([nyears,nmonths],np.int)
maxts = np.zeros([nyears,nmonths],np.int)

starttsteps = 1		# can start not at the beginning of the file, e.g. skip the first year, but assumption is that first time is beginning of Jan!
endtsteps = 46752 # 16 years
anntsteps = 2920 # timesteps per year
montsteps = [248,224,248,240,248,240,248,248,240,248,240,248]
minevent = 100000

if Data == "TRMM":
	DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
	FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'

	if Version in ["Standard","7thresh","6th_from6","5th_from48"]:
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
		DirO = DirI

		FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

	elif Version == '5th_nanto25':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'

		FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_5th_nanto25.nc'

	elif Version == '5th_nantozero':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'

		FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_5th_n2zero.nc'
	else:
		sys.exit('unexpected Version')
elif Data == "ERAI":
        DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
        FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_preprocess.nc'

	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	DirO = DirI
	FileI1 = 'Precip_Sizes_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

elif Data == "CESM":
        DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
        FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '.nc'

        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
	DirO = DirI
        FileI1 = 'Precip_Sizes_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'


FileO = 'DenDirSpd_Map_monthly_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'


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


#Create new datasets
RDenMap = np.zeros((nyears,nmonths,4,nlats,nlons),np.int)
RDenEWMap = np.zeros((nyears,nmonths,4,nlats,nlons),np.int)
RDenWWMap = np.zeros((nyears,nmonths,4,nlats,nlons),np.int)
RDenStatMap = np.zeros((nyears,nmonths,4,nlats,nlons),np.int)

RDirMap = np.zeros((nyears,nmonths,4,nlats,nlons),np.float64)
RDisMap = np.zeros((nyears,nmonths,4,nlats,nlons),np.float64)
RSpdMap = np.zeros((nyears,nmonths,4,nlats,nlons),np.float64)

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('years', nyears)
ncfile.createDimension('months',nmonths)
ncfile.createDimension('size', 4)
ncfile.createDimension('lon', nlons)
ncfile.createDimension('lat', nlats)


DensityMap = ncfile.createVariable('DensityMap','f4',('years','months','size','lat','lon'),fill_value=-9999)
WestwardMap = ncfile.createVariable('WestwardDensityMap','f4',('years','months','size','lat','lon'),fill_value=-9999)
EastwardMap = ncfile.createVariable('EastwardDensityMap','f4',('years','months','size','lat','lon'),fill_value=-9999)
StationaryMap = ncfile.createVariable('StationaryDensityMap','f4',('years','months','size','lat','lon'),fill_value=-9999)


DirectionMap = ncfile.createVariable('DirectionMap','f4',('years','months','size','lat','lon'),fill_value=-9999)
DistanceMap = ncfile.createVariable('DistanceMap','f4',('years','months','size','lat','lon'),fill_value=-9999)
SpeedMap = ncfile.createVariable('SpeedMap','f4',('years','months','size','lat','lon'),fill_value=-9999)



Longitude = ncfile.createVariable('Longitude','f4',('lon'),fill_value=-9999)
Latitude = ncfile.createVariable('Latitude','f4',('lat'),fill_value=-9999)
Years = ncfile.createVariable('years','f4',('years'),fill_value=-9999)
Months = ncfile.createVariable('months','f4',('months'),fill_value=-9999)


setattr(DensityMap,'Extra Info','Size based on timespan: 0 is < 8 (1 day), 1: < 16 (2 days) < 48 (6 days), 2: > 48)')
 
Longitude[:] = lons
Latitude[:] = lats

print DirI + FileI1

datain = xray.open_dataset(DirI + FileI1)

tspan=datain['timespan'].values
tmean = datain['tmean']
tstart = datain['tstart'].values
xmin = datain['xmin']
xmax = datain['xmax']
ymin = datain['ymin']
ymax = datain['ymax']

xstart = datain['xcenterstart']
xend = datain['xcenterend']
ystart = datain['ycenterstart']
yend = datain['ycenterend']

xcenter = datain['xcentermean'].values
ycenter = datain['ycentermean'].values

nevents = tspan.shape[0]
print nevents

ievent = 1

curtime = tstart[starttsteps]
n = 0
iyear = 0
imonth = 0
# find first timestep as defined above. 0 if beginning of file
for n in range(0,nevents):
	curtime = tstart[n]
	if curtime == starttsteps:
		mints[iyear,imonth] = n
		print('mints',iyear,imonth,n)
		break

# Loop through years and months to find the first and last timestep for each month
totaldays = starttsteps		# keep track of days of simulation

for iyear in range(0,nyears):
        Years[iyear] = anstartyr + iyear

	for imonth in range(0,nmonthsrun):
		Years[iyear] = anstartyr + iyear
		Months[imonth] = imonth + 1
		totaldays = totaldays + montsteps[imonth]
		for n in range(int(mints[iyear,imonth]),nevents):	# Start loop from beginning of last
			curtime = tstart[n]
			if curtime == totaldays:		# End of this month
				if imonth < 11:
					mints[iyear,imonth+1] = n
				else:
					if iyear < nyears-1:
						mints[iyear+1,0] = n					

				maxts[iyear,imonth] = n  #switched to n instead of n-1 because python doesn't use the last index in a range
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

print angle[4:6]
print distance[4:6]
print tspan[4:6]

# Loop through years and months to count the number
for iyear in range(0,nyears):
	print iyear
	for imonth in range(0,nmonthsrun):
		print imonth
		print mints[iyear,imonth]
		print maxts[iyear,imonth]
	
		for ievent in range(mints[iyear,imonth],maxts[iyear,imonth]):
			
			if (ievent % 1000 == 0):
				print "ievent: " + str(ievent)	# keep track of program progress	
			if tspan[ievent] <= 8:
				sizeindex = 0
                        elif tspan[ievent] < 16:
				sizeindex = 1
                        elif tspan[ievent] < 48:
				sizeindex = 2
			else:
				sizeindex = 3

			RDenMap[iyear,imonth,sizeindex,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
	
			if distance[ievent] == 0.0:
				RDenStatMap[iyear,imonth,sizeindex,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
			elif angle[ievent] < 1.57 and angle[ievent] > -1.57:
	                        RDenEWMap[iyear,imonth,sizeindex,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
			else:	
				RDenWWMap[iyear,imonth,sizeindex,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1

			RDirMap[iyear,imonth,sizeindex,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += angle[ievent]
			RDisMap[iyear,imonth,sizeindex,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += distance[ievent]
			RSpdMap[iyear,imonth,sizeindex,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += distance[ievent]/(tspan[ievent]*3.0*60.0*60.0)	# calculate speed in m/s	

DensityMap[:,:,:,:,:] = RDenMap
#RDenMap[np.where(RDenMap == 0)] = np.nan


RDisMap = RDisMap / RDenMap
RDirMap = RDirMap / RDenMap
RSpdMap = RSpdMap / RDenMap

RDisMap[np.where(RDenMap == 0)] = np.nan
RDirMap[np.where(RDenMap == 0)] = np.nan
RSpdMap[np.where(RDenMap == 0)] = np.nan


DistanceMap[:,:,:,:,:] = RDisMap
DirectionMap[:,:,:,:,:] = RDirMap
SpeedMap[:,:,:,:,:] = RSpdMap
EastwardMap[:,:,:,:,:] = RDenEWMap
WestwardMap[:,:,:,:,:] = RDenWWMap



datain.close()

Ngl.end()









