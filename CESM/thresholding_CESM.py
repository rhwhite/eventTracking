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
import xray
import pandas
from shutil import copyfile

#Thresholds are in mm/day, data is in mm/hr
# Switched code so now thresholds should be from smallest to largest
# Values lower than the 1st threshold are set to 0
#Threshold = [120, 80, 56, 40, 24, -9999]
# Threshold = [24, 40, 56, 80, 120]
#Threshold = [6,12,24,48,72,96]
Threshold = [24, 40, 56, 80, 120]
#Threshold = [48,72,96,120,168] "5th_from48_"
nThresh = 5
name = "Standard_"

styear = 1980
edyear = 2012
# Convert thresholds from mm/day to m/s
Thresh = [x/(1000.0 * 60.0 * 60.0 * 24.0) for x in Threshold]

print Thresh

Dir = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
DirO = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/Thresholds/'

for year in range(styear,edyear+1):
	print year
	stryear = str(year)

	Filein = 'f.e13.FAMIPC5.ne120_ne120.1979_2012.001.cam.h4.PRECT.' + stryear + '010100Z-' + stryear + '123121Z.nc'

	print Dir + Filein
	precipin = xray.open_dataset(Dir + Filein)

	latsin = precipin['lat'][:]
	lonsin = precipin['lon'][:]
	precipnew = precipin['PRECT'][:]
	timesin = precipin['time'][:]
	times = pandas.to_datetime(precipnew['time'].values)
	print times
	years = times.year
	hours = times.hour

	nlons = len(lonsin)
	nlats = len(latsin)
	ntimes = len(years)

	precipchunk = precipnew.sel(time=str(year))
	timeschunk = pandas.to_datetime(precipchunk['time'].values)
	hours = timeschunk.hour
	print hours
	# reset counter for each year, but don't start at 0 if the first hour isn't 0!
	counter = hours[0]/3
	print counter
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
		FileO = 'CESM_f.e13.FAMPIC5.ne120_ne120.001.3hrly_' + name + stryear + '_' + filecount.zfill(5) + '.nc'
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
	 
