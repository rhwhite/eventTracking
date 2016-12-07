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

#Thresholds are in mm/day, data is in mm/hr
# Switched code so now thresholds should be from smallest to largest
# Values lower than the 1st threshold are set to 0
#Threshold = [120, 80, 56, 40, 24, -9999]
Threshold = [24, 40, 56, 80, 120]
#Threshold = [6,12,24,48,72,96]
#Threshold = [48,72,96,120,168]
nThresh = 5
name = "Std_ERAIgrid"

Thresh = [x/24.0 for x in Threshold]

print Thresh

Dir = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
DirO = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/Thresholds/'

Filein = 'regrid2ERAI_TRMM_3B42_1998-2014.nc'

precipin = xray.open_dataset(Dir + Filein)

latsin = precipin['latitude'][:]
lonsin = precipin['longitude'][:]
timesin = precipin['time'][:]

nlons = len(lonsin)
nlats = len(latsin)
for year in range(1998,2015):    #(1998,2015)
	print year
	stryear = str(year)
	counter = 0	#reset counter for each year

	precip_thrsh=precipin['pcp']

	precipchunk = precip_thrsh.sel(time=str(year))
	precipthresh = np.zeros(precipchunk.shape,np.float)
	timeschunk = timesin.sel(time=str(year))
	ntimes = len(timeschunk)
	
	for x,y in np.nditer([precipchunk,precipthresh],op_flags=['readwrite']):
		if x < Thresh[0]:	# Deal with values below 1st threshold (0)
			y[...] = 0
		elif x >= Thresh[nThresh-1]:	# Deal with values above last threshold (nThresh - 1)
			y[...] = nThresh
		else:
			for ith in range(0,nThresh-1):	# iterate through thresholds 0 to nThresh - 2
				if x >= Thresh[ith] and x < Thresh[ith+1]:
					y[...] = ith + 1
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
		FileO = 'TRMMregrid_3B42_3hrly_' + name + stryear + '_' + filecount.zfill(5) + '.nc'
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

		precipTH[:,:,:] = precipthresh[count,:,:]
		precipTH.grid_type = 'linear'
		precipTH.grid_name = 'grid-1'


		#global attribute
		ncfile.Conventions = 'COARDS'
		ncfile.calendar = 'standard'
		ncfile.center = 'gsfc'
		ncfile.History = 'Thresholds: ' + str(Threshold)

		ncfile.close()
		counter = counter + 1		   
	 
