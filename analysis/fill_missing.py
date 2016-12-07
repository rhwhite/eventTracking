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

import xray

#Thresholds are in mm/day, data is in mm/hr


startyr = 1998
endyr = 2014


#FileP = 'TRMM_pcp_3hrly_nonan2.2005.nc'
FileP = 'TRMM_' + str(startyr) + '-' + str(endyr) + '_3B42_3hrly_nonan2.nc'

minevent = 100000

#DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/raw/'

precipdata = Dataset(DirP + FileP,'r+')

ntimes = len(precipdata['time'])

print(ntimes)
print "start looking at missing values"

# for first timestep, can't fill missing values
step = 4000

count = 0
for nt in range(0,ntimes,step):
	
	print nt
	upto = min(nt + step+2,ntimes)
	print upto
	precipin = ma.getdata(precipdata.variables['pcp'][nt:upto,:,:])

	mvs = np.where(precipin < -1)

	nvalues = len(mvs[0])
	print nvalues
	for nv in range(0,nvalues):
	        if mvs[0][nv]+count == 0 or mvs[0][nv]+count == ntimes-1:
#can't do anything for missing values at the very beginning or very end of the timeseries
                #print "can't deal with missing values sequentially in time, setting to zero"
        	        precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.0

		elif mvs[0][nv] < step+1:
        		if (precipin[mvs[0][nv]-1,mvs[1][nv],mvs[2][nv]] > -1 and precipin[mvs[0][nv]+1,mvs[1][nv],mvs[2][nv]] > -1):
        	        	precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.5 * precipin[mvs[0][nv]-1,mvs[1][nv],mvs[2][nv]] + \
        	                	0.5 * precipin[mvs[0][nv]+1,mvs[1][nv],mvs[2][nv]]
			else:
        	        	precipin[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.0 
		else:
			#don't do anything for the last timestamp in each segment - it will be included in the next unless it's the
			#last one in which case we've already dealt with it!
			print "not doing anything for time " + str(count + mvs[0][nv])
	                
	count = count + step

	mvs2 = np.where(precipin[0:step+1,:,:] < -1)
	nvalues2 = len(mvs2[0])
	if nvalues2 != 0:
	        sys.err("error with filling missing data; missing data remains")
	        
	precipdata.variables['pcp'][nt:upto,:,:] = precipin
	print 'finished filling in missing data'
	
			
