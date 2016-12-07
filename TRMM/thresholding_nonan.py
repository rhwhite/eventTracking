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
from rhwhitepackages.query import query_yes_no
from rhwhitepackages.parallel import calcmissed
from multiprocessing import Process,Queue
from sys import getsizeof
import psutil

#Thresholds are in mm/day, data is in mm/hr
# Switched code so now thresholds should be from smallest to largest
# Values lower than the 1st threshold are set to 0
#Threshold = [120, 80, 56, 40, 24, -9999]
# Threshold = [24, 40, 56, 80, 120]
#Threshold = [6,12,24,48,72,96]
#Threshold = [48,72,96,120,168]

lowmem = True
name = "6th_from6" #"Standard"
startyr = 1998
endyr = 1999

if name == "5th_from48_":
        Threshold = [48,72,96,120,168]
elif name == "Standard":
        Threshold = [24, 40, 56, 80, 120]
elif name == "6th_from6":
	Threshold = [6,12,24,48,72,96]
else:
        error("unknown name")

nThresh = len(Threshold)
Thresh = [x/24.0 for x in Threshold]

print Thresh
threshlow = Thresh[0]
Dir = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
DirO = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/Thresholds/'

Filein = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'

precipin = xray.open_dataset(Dir + Filein,slice=time(0,100)chunks={'time':10})

latsin = precipin['latitude'][:]
lonsin = precipin['longitude'][:]
precipnew=precipin['pcp']
#print precipnew.time

nlons = len(lonsin)
nlats = len(latsin)

print "Hold on, I'm calculating the total amount of precip missed by your chosen thresholds, estimating the value using the first year"
check = query_yes_no("Do you want to confirm whether or not it's acceptable?")

print(psutil.phymem_usage())

#	result = 100.0* missedvalues/np.sum(precipchunk)
#        q.put(np.ma.sum(np.ma.masked_greater(precipnew.sel(time=str(yearin)+'-' + '{:02d}'.format(monthin)),Thresh[0])))
#	q2.put(np.nansum(precipnew.sel(time=str(yearin)+'-' + '{:02d}'.format(monthin))))
results = []
results2 = []

startTime = time.time()

total = np.sum(precipin)
caughtvalues = np.ma.sum(precipin[precipin >=Thresh[0]])
	
totalmissed = (total - caughtvalues)
percmissed = 100.0* (totalmissed)/total

endTime = time.time()
workTime = endTime - startTime
print "job took " + str(workTime) + " seconds"
print "this threshold misses "  + str(totalmissed) + " which is " + str(percmissed) + " % of the total rain"
if check:
	accept = query_yes_no("is this acceptable?")
	if not accept: exit("you didn't accept the thresholds")
exit()


timesin = precipin['time'][:]

for year in range(startyr,endyr):    #(1998,2015)
	print year
	stryear = str(year)
	counter = 0	#reset counter for each year

	precipchunk = precipnew.sel(time=str(year))
	timeschunk = timesin.sel(time=str(year))
	ntimes = len(timeschunk)
	
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
		FileO = 'TRMM_3B42_3hrly_' + name + stryear + '_' + filecount.zfill(5) + '.nc'
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
	 
