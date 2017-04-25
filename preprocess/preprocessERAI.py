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

styear = [1998]
edyear = [2015]

Dir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'

nfiles = len(styear)

#Get dimensions and dimension sizes for new file
Filein = 'ERAI_Totalprecip_' + str(styear[0]) + '-' + str(edyear[0]) + '.nc'
Fileout = 'ERAI_Totalprecip_' + str(styear[0]) + '-' + str(edyear[0]) + 'preprocess.nc'
precipin = xray.open_dataset(Dir + Filein)

latsin = precipin['latitude'][:]
lonsin = precipin['longitude'][:]

nlats = len(latsin)
nlons = len(lonsin)

copyfile(Dir + Filein, Dir + Fileout)

writefile = Dataset(Dir + Fileout,'r+')
tpout = writefile.createVariable('pcp','f8',('time','latitude','longitude'))

timesin = precipin['time'][:]
precip = precipin['tp'][:]
precipvalues = np.asfarray(precipin['tp'][:],dtype='float') 
times = pandas.to_datetime(precip['time'].values)

hours = times.hour
years = times.year
print years
print precip.shape
nlons = len(lonsin)
nlats = len(latsin)
ntimes = len(hours)
#Loop through each timestamp, and subtract previous value if necessary, to make all of them 3hrly
precipnew = np.zeros_like(precipvalues)


for itime in range(1,ntimes):
	if hours[itime] % 12 == 3:
		precipnew[itime,:,:] = precipvalues[itime,:,:]
	else:	
		precipnew[itime,:,:] = precipvalues[itime,:,:] - precipvalues[itime-1,:,:]
	# write out to file and convert to mm/hour from m / 3 hours
	tpout[itime,:,:] = precipnew[itime,:,:] * 1000.0/3.0
	if itime % 1000 = 0:
		print 'itime now ', itime
# For now just set the first time to its original value/4, as we don't have the correct values in the file to calculate the absolute
tpout[0,:,:] = precipvalues[0,:,:] *1000.0 / (3.0 * 4.0)

precipin.close()
writefile.close()
print('done here!')

