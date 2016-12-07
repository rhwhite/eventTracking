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
import psutil
from rhwhitepackages.query import query_yes_no

name = "Standard_"

styears = [1980]
edyears = [2011]

firsttimestep = 9	# 0 forecast hour

nsets = len(styears)

Dir = '/home/disk/eos4/rachel/Obs/ERA_20C/'
DirO = '/home/disk/eos4/rachel/Obs/ERA_20C/Thresholds/'

for idec in range(0,nsets):

	Filein2 = 'ERA_20C_Totalprecip_' + str(styears[idec]) + '-' + str(edyears[idec]) + '_up.nc'
	Filein = 'ERA_20C_Totalprecip_' + str(styears[idec]) + '-' + str(edyears[idec]) + '_preprocess.nc'

	try:
		precipin = xray.open_dataset(Dir + Filein)
		check = query_yes_no("The input file you specified already exists. Are you sure you want to overwrite it?")
		if not check:
			exit('Ok, exiting')

	
	except RuntimeError:
		pass

	copyfile(Dir + Filein2, Dir + Filein)
	precipin = xray.open_dataset(Dir + Filein2)

	latsin = precipin['latitude'][:]
	lonsin = precipin['longitude'][:]
	timesin = precipin['time'].values
	precip = precipin['tp'][:,:,:]

	times = pandas.to_datetime(timesin)
	#times = pandas.to_datetime(precip['time'].values)
	
	hours = times.hour
	years = times.year
	print years

	nlons = len(lonsin)
	nlats = len(latsin)
	ntimes = len(hours)
	print "ntimes = ", ntimes
	#Loop through each timestamp, and subtract previous value if necessary, to make all of them 3hrly
	
	writefile = Dataset(Dir + Filein,'r+')
	tpnew = writefile.createVariable('tpnew','f8',('time','latitude','longitude'))
	#check that the first time doesn't need to be processed
	if hours[0] % 24 != firsttimestep:
		sys.error("need to preprocess so that the first timestep in a file is hour " + str(firsttimestep) + ", so it doesn't neeed to be processed")	
	# Calc first for averaging:
	tpnew[0,:,:] = (1000.0/3.0) * precip[0,:,:].values

	for itime in range(1,ntimes):
		if hours[itime] % 24 == firsttimestep:
			# If it's a new forecast step, then average between this and last to get rid of discontinuity every 24 hours
			precipnew = (1000.0/3.0) * 0.5 * (precip[itime,:,:]+precip[itime-1,:,:]-precip[itime-2,:,:])
		else:	
			precipnew = (1000.0/3.0) * (precip[itime,:,:] - precip[itime-1,:,:])	# calculate 3 hourly accumulations
		if itime % 1000 == 0:
			print 'itime update:', itime
		tpnew[itime,:,:] = precipnew.values	# average of previous and next 3 hourly accumulations
		# Update precip prev:


	precipin.close()
	print('done here!')

	# convert from m (in 3 hours) to mm/hr
	tpnew.units = "mm/hr"
	writefile.close()

	del precip
	del precipin


