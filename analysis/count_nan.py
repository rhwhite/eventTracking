# -*- coding: utf-8 -*-
"""
Script file to count the number of NaNs in the raw precipitation data to check
how these change over time
"""
import os, errno
import numpy as np
import numpy.ma as ma
import netCDF4
from netCDF4 import Dataset
import datetime as dt
#Thresholds are in mm/day, data is in mm/hr
Threshold = [140, 110, 80, 40, 30, 20, 10, 1]
#Threshold = [120, 80, 56, 40, 24, -9999]
nThresh = 7

Thresh = [x/24.0 for x in Threshold]

Dir = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/raw/'
DirO = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/Thresholds/'
startyr = 1998
endyr = 2015
nyears = endyr-startyr

nummissing = np.zeros(nyears*12)
nummissing2zero = np.zeros(nyears*12)

FileONM = 'TRMM_3B42_3hrly_' + str(startyr) + '_' + str(endyr) + '_num_missing_fast.nc'
try:
    os.remove(Dir + FileONM)
except OSError:
    pass
ncfileNM = Dataset(Dir + FileONM, 'w')
ncfileNM.createDimension('months', nyears*12)
ncfileNM.createDimension('years', nyears)

FileNumMissing = ncfileNM.createVariable('NumMissing','int32',('months'))
FileNumMissing2Zero = ncfileNM.createVariable('NumMissing2Zero','int32',('months'))
FileNumMissingY = ncfileNM.createVariable('NumMissingY','int32',('years'))
FileNumMissing2ZeroY = ncfileNM.createVariable('NumMissing2ZeroY','int32',('years'))


imonth = 0
iyear = 0
for year in range(1998,2015):    #(1998,2015)
	nmissingY = 0
	nmissing2zeroY = 0
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

		if (yearpre < 2000 or yearpre > 2010 or (yearpre == 2010 and monthpre > 9)):	
			Filepre = 'TRMM_3B42_3hrly.' + str(yearpre) + str(monthpre).zfill(2) + '.7.nc'
		else:
			Filepre = 'TRMM_3B42_3hrly.' + str(yearpre) + str(monthpre).zfill(2) + '.7A.nc'
		

                if (year < 2000 or year > 2010 or (year == 2010 and month > 9)):    
                        File1 = 'TRMM_3B42_3hrly.' + str(year) + str(month).zfill(2) + '.7.nc'
                else:
                        File1 = 'TRMM_3B42_3hrly.' + str(year) + str(month).zfill(2) + '.7A.nc'

                if (yearpost < 2000 or yearpost > 2010 or (yearpost == 2010 and monthpost > 9)):   
                        Filepost = 'TRMM_3B42_3hrly.' + str(yearpost) + str(monthpost).zfill(2) + '.7.nc'
                else:   
                        Filepost = 'TRMM_3B42_3hrly.' + str(yearpost) + str(monthpost).zfill(2) + '.7A.nc'

#		print Filepre
		print File1
#		print Filepost

		nc_fid = Dataset(Dir + File1, 'r')
		
		if (year != 1998 or month != 1):
			nc_pre = Dataset(Dir + Filepre,'r')
		if (year != 2014 or month != 12):
			nc_post = Dataset(Dir + Filepost,'r')
		
		latsin = nc_fid.variables['latitude'][:]
		lonsin = nc_fid.variables['longitude'][:]
		timesin = nc_fid.variables['time'][:]
		precipin = ma.getdata(nc_fid.variables['pcp'][:,:,:])

		nlons = len(lonsin)
		nlats = len(latsin)
		ntimes = len(timesin)

		precip_thrsh= precipin[:,:,:]


                if (year == 1998 and month == 1):
			precip_pre = precipin[0,:,:]
			precip_pre[:,:] = -9999
			print precip_pre.shape
		else:
			ntimespre = len(nc_pre.variables['time'])
			precip_pre = ma.getdata(nc_pre.variables['pcp'][ntimespre-1,:,:])
               		print precip_pre.shape 
	                nc_pre.close()

		if (year == 2014 and month == 12):
			precip_post = precipin[0,:,:]
			precip_post[:,:] = -9999
			print precip_post.shape
		else:
			precip_post = ma.getdata(nc_post.variables['pcp'][0,:,:])
			print precip_post.shape
			nc_post.close()	
		nc_fid.close()

		print "start looking at missing values"

		mvs = np.where(precipin < -1)

		nvalues = len(mvs[0])
		print nvalues
		nmissing = nmissing + nvalues
                nmissingY = nmissingY + nvalues

		for nv in range(0,nvalues):
			if mvs[0][nv] == 0 or mvs[0][nv] == ntimes-1:
				if precip_post[mvs[1][nv],mvs[2][nv]] < -1 or precip_pre[mvs[1][nv],mvs[2][nv]] < -1:
					nmissing2zero = nmissing2zero + 1
					nmissing2zeroY = nmissing2zeroY + 1
					precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.0
				else:
					precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = \
						0.5 * precip_pre[mvs[1][nv],mvs[2][nv]] + \
						0.5 * precip_post[mvs[1][nv],mvs[2][nv]]
			elif (precipin[mvs[0][nv]-1,mvs[1][nv],mvs[2][nv]] > -1 and precipin[mvs[0][nv]+1,mvs[1][nv],mvs[2][nv]] > -1):
				precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.5 * precipin[mvs[0][nv]-1,mvs[1][nv],mvs[2][nv]] + \
					0.5 * precipin[mvs[0][nv]+1,mvs[1][nv],mvs[2][nv]]
			else:
				precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.0
                                nmissing2zero = nmissing2zero + 1
				nmissing2zeroY = nmissing2zeroY + 1
		#	if precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] != 0.0:
		#		print precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]]


		mvs2 = np.where(precipin < -1)
		nvalues2 = len(mvs2[0])
		if nvalues2 != 0:
			sys.err("error with filling missing data; missing data remains")


		print 'finished filling in missing data'
		"""		
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
			FileO = 'TRMM_3B42_3hrly_' + str(nThresh) + 'thresh_' + stryear + '_' + filecount.zfill(5) + '.nc'
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
		"""

		FileNumMissing[imonth] = nmissing
		FileNumMissing2Zero[imonth] = nmissing2zero
	
		imonth = imonth + 1
		print imonth

	FileNumMissingY[iyear] = nmissingY
        FileNumMissing2ZeroY[iyear] = nmissing2zeroY
	iyear = iyear + 1

ncfileNM.close()
 
