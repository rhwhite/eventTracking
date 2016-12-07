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

Version = "Standard"

filetimespan = "3hrly"

if filetimespan == "3hrly":
	mult = 3.0

chunksize = 2000
latlonaddsize = 2

startyr = 1998
endyr = 2014 


#DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
#FileP = 'TRMM_pcp_3hrly_nonan.2005.nc'
DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
FileP = 'TRMM_' + str(startyr) + '-' + str(endyr) + '_3B42_3hrly_nonan.nc'

minevent = 100000

if Version == "Standard":
	Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/'
	DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
	File1 = 'ts_TRMM1998-2014_final_4Dobjects.nc'
	FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_Standard.nc'
	TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_final_4Dobject_tree.txt'
elif Version == "7thresh":
	Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/'
	DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip'
	File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_7thresh_4Dobjects.nc'
	TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_7thresh_4Dobject_tree.txt'
	FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_7thresh.nc'
elif Version == "5thresh_n2z":
	Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/'
	DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
	File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_5th_n2z_4Dobjects.nc'
	TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_5th_n2z_4Dobject_tree.txt'
	FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_5th_n2zero.nc'
elif Version == "6th_from6" or Version == "5th_from48":
        Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
        File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
        TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
        FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
else:
	sys.exit("Version " + Version + "doesn't exist here!")

# If output directory doesn't exist, create it!
if not os.path.exists(DirO):
        os.makedirs(DirO)

textdata = np.loadtxt(Dir + TxtFileIn,skiprows = 1, usecols = (0,))
print textdata.shape
print textdata[0]
print textdata[-1]

event1 = textdata[0]
nevents = textdata[-1] - textdata[0] + 1
print("nevents: ", nevents)


minx = np.zeros(nevents)
maxx = np.zeros(nevents)
miny = np.zeros(nevents)
maxy = np.zeros(nevents)
mint = np.zeros(nevents)
maxt = np.zeros(nevents)
startx = np.zeros(nevents)
endx = np.zeros(nevents)
starty = np.zeros(nevents)
endy = np.zeros(nevents)

#centerxstart = np.zeros(nevents)
#centerystart = np.zeros(nevents)

centerx = np.zeros(nevents)
centery = np.zeros(nevents)
meant = np.zeros(nevents)

nlines = len(textdata)
print("nlines: ", nlines)

preeventnum = -100

with open(Dir + TxtFileIn,"r") as textFile:
	next(textFile)	#Skip header line

	for lines in textFile:	#loop through all lines
		line = lines.split('\t')
	        eventnum = int(float(line[0]))
        	if (eventnum % 10000 == 0):
                	print "eventnum: " + str(eventnum)
		if eventnum == preeventnum:
			listcenterx.append(int(re.findall(r'\d+',line[6])[0]))
			listcentery.append(int(re.findall(r'\d+',line[6])[1]))
		else:
			if (preeventnum > 0): #then it's not the very first event
				#new event: take mean of center x, center y and t and put into array before updating index
				startx[index] = listcenterx[0]
				starty[index] = listcentery[0]
				endx[index] = listcenterx[-1]
				endy[index] = listcentery[-1]
							
			#refresh list
			#update index for new event
			index = eventnum-event1
			preeventnum = eventnum

                        listcenterx = [int(re.findall(r'\d+',line[6])[0])]
                        listcentery = [int(re.findall(r'\d+',line[6])[1])]

textFile.close()

print 'centerx'
print centerx[0:10]
print 'centery'
print centery[0:10]
print 'meant'
print meant[0:10]
print 'startx'
print startx[0:10]
print 'endx'
print endx[0:10]
print 'starty'
print starty[0:10]
print 'endy'
print endy[0:10]


eventsdata = xray.open_dataset(Dir + File1)

print(eventsdata.coords)

eventsin=eventsdata['value']

eventsin = eventsin.squeeze(dim="z")

ntimes = eventsin.shape[0]
nlats = eventsin.shape[1]
nlons = eventsin.shape[2]
maxevent = int(np.amax(eventsin[ntimes-10:ntimes,:,:]))

nevents2 = maxevent - minevent + 1


if (nevents2 != nevents):
	print nevents2, nevents
	sys.exit("event numbers not equal in netcdf and text file")


ncfile = Dataset(DirO + FileO, 'r+')

xstartmean = ncfile.createVariable('xcenterstart','f4',('events'),fill_value=-9999)
xendmean = ncfile.createVariable('xcenterend','f4',('events'),fill_value=-9999)
ystartmean = ncfile.createVariable('ycenterstart','f4',('events'),fill_value=-9999)
yendmean = ncfile.createVariable('ycenterend','f4',('events'),fill_value=-9999)

print "starting now"


#for ievent in range(minevent, minevent + 3): ### For testing!!!
for ievent in range(minevent,maxevent + 1):

        index = ievent-minevent
	if (ievent % 5000 == 0):
		print "ievent: " + str(ievent)

	xstartmean[ievent-minevent] = startx[index]
        xendmean[ievent-minevent] = endx[index]

        ystartmean[ievent-minevent] = starty[index]
        yendmean[ievent-minevent] = endy[index]

#results = results2
ncfile.close()

