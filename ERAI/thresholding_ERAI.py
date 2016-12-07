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
import time
import xray
import pandas
from shutil import copyfile
from rhwhitepackages.query import query_yes_no

#Thresholds are in mm/day, data is in mm/hr
# Switched code so now thresholds should be from smallest to largest
# Values lower than the 1st threshold are set to 0
#Threshold = [120, 80, 56, 40, 24, -9999]
# Threshold = [24, 40, 56, 80, 120]
#Threshold = [6,12,24,48,72,96]
#Threshold = [24, 40, 56, 80, 120] "Standard"
#Threshold = [48,72,96,120,168] "5th_from48_"
nThresh = 5
name = "ERAvlow"
chunksize = 1000


if name == "ERAlow":
        Threshold = [12, 18, 26, 36, 48]
elif name == "ERAvlow":
	Threshold = [2, 4, 8, 16, 32, 48] 
else:
	error("unknown name")
print Threshold
styear = 1980
edyear = 2015

Thresh = len(Threshold)
Thresh = [x/24.0 for x in Threshold]

Dir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
DirO = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/Thresholds/'

Filein2 = 'ERAI_Totalprecip_' + str(styear) + '-' + str(edyear) + '.nc'
Filein = 'ERAI_Totalprecip_' + str(styear) + '-' + str(edyear) + '_preprocess.nc'

try:
	precipin = xray.open_dataset(Dir + Filein)

except RuntimeError:
	print "need to process files first"

	copyfile(Dir + Filein2, Dir + Filein)
	precipin = xray.open_dataset(Dir + Filein2)

	latsin = precipin['latitude'][:]
	lonsin = precipin['longitude'][:]
	timesin = precipin['time'][:]
	precip = precipin['tp'][:]

	times = pandas.to_datetime(precip['time'].values)
	
	hours = times.hour
	years = times.year
	print years

	nlons = len(lonsin)
	nlats = len(latsin)
	ntimes = len(hours)
	#Loop through each timestamp, and subtract previous value if necessary, to make all of them 3hrly
	precipnew = np.zeros_like(precip)
	
	#check that the first time doesn't need to be processed - can't do this!
	if hours[0] % 12 != 3:
		sys.error("need to preprocess so that the first timestep in a file is either hour 03 or 15, so it doesn't neeed to be processed")


	for itime in range(0,ntimes):
		if hours[itime] % 12 == 3:
			precipnew[itime,:,:] = precip[itime,:,:]
		else:	
			precipnew[itime,:,:] = precip[itime,:,:] - precip[itime-1,:,:]
	precipin.close()
	print('done here!')

	# convert from m (in 3 hours) to mm/hr
	precipnew = precipnew * 1000.0/3.0
	writefile = Dataset(Dir + Filein,'r+')
	tpnew = writefile.createVariable('tpnew','f8',('time','latitude','longitude'))
	tpnew[...] = precipnew[...]
	writefile.close()

	del precip
	del precipin
        precipin = xray.open_dataset(Dir + Filein)


latsin = precipin['latitude'][:]
lonsin = precipin['longitude'][:]
precipnew = precipin['tpnew'] #[:]
#timesin = precipin['time'][:]
#times = pandas.to_datetime(precipnew['time'].values)
#print times
#years = times.year
#hours = times.hour

nlons = len(lonsin)
nlats = len(latsin)
ntimes = len(precipnew.time)
print ntimes

print "Hold on, I'm calculating the total amount of precip missed by your chosen thresholds, estimating the value using the first year"
check = query_yes_no("Do you want to confirm whether or not it's acceptable?")

lowestTH = Thresh[0]
totals = []
missed = []

startTime = time.time()

nits = np.ceil(ntimes/chunksize)
print str(nits) + "iterations required"
for i in range(0,np.int(nits)):
        print i
	for attempt in range(10):
		try:
			precipchunk = np.float32(precipnew.isel(time=slice(i*chunksize,np.amin([(i+1)*chunksize,ntimes]))))
			total = np.sum(precipchunk)
			missedvalues = np.ma.sum(np.ma.masked_greater_equal(precipchunk,Thresh[0]))
			totals.append(total)
			missed.append(missedvalues)
		

print totals
totalmissed = 100.0* sum(missed)/sum(totals)

endTime = time.time()
workTime = endTime - startTime
print "job took " + str(workTime) + " seconds"

print missed
print "this threshold misses " + str(sum(missed)) + " which is " + str(totalmissed) + " % of the total rain"

if check:
        accept = query_yes_no("is this acceptable?")
        if not accept: exit("you didn't accept the thresholds")
exit()

for year in range(styear,edyear + 2):    #(1998,2015)
	print year
	stryear = str(year)

	# find start and end points for this year
#	if years[0] == year:
#		startitime = 0
#	if years[-1] == year:
#		enditime = ntimes-1
#	for itime in range(1,ntimes-1):
#		if years[itime] == year and years[itime-1] < year:
#			startitime = itime
#
#		if years[itime] == year and years[itime+1] > year:
#			enditime = itime

#	print startitime, enditime

	precipchunk = precipnew.sel(time=str(year))
	timeschunk = pandas.to_datetime(precipchunk['time'].values)
	hours = timeschunk.hour
	# reset counter for each year, but don't start at 0 if the first hour isn't 0!
	counter = hours[0]/3
	ntimes = len(hours)
	for x in np.nditer(precipchunk, op_flags=['readwrite']):
		if x < Thresh[0]:	# Deal with values below 1st threshold (0)
			x[...] = 0
		elif x >= Thresh[nThresh-1]:	# Deal with values above last threshold (nThresh - 1)
			x[...] = nThresh
		else:
			for ith in range(0,nThresh-1):	# iterate through thresholds 0 to nThresh - 2
				if x >= Thresh[ith] and x < Thresh[ith+1]:
					x[...] = ith + 1
					break

#		for y in np.nditer(precip_thrsh):
#		    if y > 5:
#			print 'missing value'
#			print y
#		    elif y < 0:
#			print 'weird value'
#			print y

	maxcount = ntimes
	print maxcount
	for count in range(0,maxcount):
		filecount = str(counter)
		FileO = 'ERAI_3hrly_' + name + "_" + stryear + '_' + filecount.zfill(5) + '.nc'
		try:
		    os.remove(DirO + FileO)
		except OSError:
		    pass       

		ncfile = Dataset(DirO + FileO, 'w')
		ncfile.createDimension('z', 1)
		ncfile.createDimension('y',len(latsin))
		ncfile.createDimension('x',len(lonsin))
		ncfile.createDimension('time',1)

		precipTH = ncfile.createVariable('value','f4',('z','y','x'),fill_value=0)

		precipTH[:,:,:] = precipchunk[count,:,:].values
		precipTH.grid_type = 'linear'
		precipTH.grid_name = 'grid-1'


		#global attribute
		ncfile.Conventions = 'COARDS'
		ncfile.calendar = 'standard'
		ncfile.center = 'gsfc'
		ncfile.History = 'Thresholds: ' + str(Threshold)

		ncfile.close()
		counter = counter + 1		   
	 
