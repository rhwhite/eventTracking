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
#Thresholds are in mm/day, data is in mm/hr
#Threshold = [140, 110, 80, 40, 30, 20, 10, 1]
Threshold = [120, 80, 56, 40, 24, -9999]
nThresh = 5

Thresh = [x/24.0 for x in Threshold]

Dir = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/raw/'
DirO = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/Thresholds/'
startyr = 1998
endyr = 2015
nyears = endyr-startyr

nummissing = np.zeros(nyears*12)
nummissing2zero = np.zeros(nyears*12)

imonth = 0
for year in range(1998,2015):    #(1998,2015)
	stryear = str(year)
	counter = 0	#reset counter for each year
	for month in range(1,13):
	        nmissing = 0 # reset missing data counter for each month
		nmissing2zero = 0

		monthpre = month-1
		monthpost = month+1
		yearpre = year
		yearpost = year
		if month == 1:
			monthpre = 12
			yearpre = year-1 
		elif month == 12:
			monthpost = 1
			yearpost = year+1

                if (year < 2000 or year > 2010 or (year == 2010 and month > 9)):    
                        File1 = 'TRMM_3B42_3hrly.' + str(year) + str(month).zfill(2) + '.7.nc'
                else:
                        File1 = 'TRMM_3B42_3hrly.' + str(year) + str(month).zfill(2) + '.7A.nc'

		print File1

		nc_fid = Dataset(Dir + File1, 'r')

		latsin = nc_fid.variables['latitude'][:]
		lonsin = nc_fid.variables['longitude'][:]
		timesin = nc_fid.variables['time'][:]
		precipin = ma.getdata(nc_fid.variables['pcp'][:,:,:])

		nlons = len(lonsin)
		nlats = len(latsin)
		ntimes = len(timesin)
		print ntimes

		precip_thrsh= precipin[:,:,:]

		nc_fid.close()

		print "start looking at missing values"

		mvs = np.where(precipin < -1)
		nvalues = len(mvs[0])
		print nvalues
		nmissing = nmissing + nvalues
	
                precipin[np.where(precipin < -1)] = 0
	
		mvs2 = np.where(precipin < -1)
		nvalues2 = len(mvs2[0])
		if nvalues2 != 0:
			sys.err("error with filling missing data; missing data remains")


		print 'finished filling in missing data'
				
		for x in np.nditer(precip_thrsh, op_flags=['readwrite']):
			for ith in range(0,nThresh+1):
				if x >= Thresh[ith]:
					x[...] = nThresh-ith
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
			FileO = 'nan2zero_TRMM_3B42_3hrly_' + str(nThresh) + 'thresh_' + stryear + '_' + filecount.zfill(5) + '.nc'
			try:
			    os.remove(DirO + FileO)
			except OSError:
			    pass       

			ncfile = Dataset(DirO + FileO, 'w')
			ncfile.createDimension('z', 1)
			ncfile.createDimension('y',len(latsin))
			ncfile.createDimension('x',len(lonsin))

			precipTH = ncfile.createVariable('value','f4',('z','y','x'),fill_value=0)

			precipTH[:,:,:] = precip_thrsh[count,:,:]
			precipTH.grid_type = 'linear'
			precipTH.grid_name = 'grid-1'
			#global attribute
			ncfile.Conventions = 'COARDS'
			ncfile.calendar = 'standard'
			ncfile.center = 'gsfc'
			ncfile.History = 'Thresholds: ' + str(Threshold)
			ncfile.close()
			counter = counter + 1		   

	imonth = imonth + 1
	print imonth


 
