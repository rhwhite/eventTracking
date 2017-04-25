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

name = "Standard_"

styears = [1980,1990,2000,2010]
edyears = [1989,1999,2009,2015]

nsets = len(styears)

Dir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
DirO = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/Thresholds/'

for idec in range(0,nsets):

	Filein2 = 'ERAI_Totalprecip_' + str(styears[idec]) + '-' + str(edyears[idec]) + '.nc'
	Filein = 'ERAI_Totalprecip_' + str(styears[idec]) + '-' + str(edyears[idec]) + '_preprocess.nc'
	try:
		precipin = xray.open_dataset(Dir + Filein)
	
	except RuntimeError:

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

		#check that the first time doesn't need to be processed
		if hours[0] % 12 != 3:
			sys.error("need to preprocess so that the first timestep in a file is either hour 03 or 15, so it doesn't neeed to be processed")
		
		for itime in range(0,5): #ntimes):
			if hours[itime] % 12 == 3:
				precipnew[itime,:,:] = precip[itime,:,:]
			else:	
				precipnew[itime,:,:] = precip[itime,:,:] - precip[itime-1,:,:]
			if itime % 1000 == 0:
				print 'itime update:', itime
		precipin.close()
		print('done here!')

		# convert from m (in 3 hours) to mm/hr
		precipnew = precipnew * 1000.0/3.0
		writefile = Dataset(Dir + Filein,'r+')
		tpnew = writefile.createVariable('tpnew','f8',('time','latitude','longitude'))
		tpnew[...] = precipnew[...]
		tpnew.units = "mm/hr"
		writefile.close()

		del precip
		del precipin


