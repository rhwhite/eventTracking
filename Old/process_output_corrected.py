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
#Thresholds are in mm/day, data is in mm/hr


startyr = 1998
endyr = 2014 

DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
FileP1 = 'TRMM_pcp_3hrly_nonan.2005.nc'

minevent = 100000

Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Old/'
DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Old/Precip/'
File1 = 'ts_TRMMtest_5th_2005_4Dobjects.nc'
FileO = 'Precip_Sizes_2005_test_old_process_corrected.nc'

TxtFileIn = 'TRMMtest_5th_2005_4Dobject_tree.txt'

with open(Dir + TxtFileIn,"r") as textFile:
	lines = [line.split('\t') for line in textFile]

nlines = len(lines)
print("nlines: ", nlines)
event1 = int(lines[1][0])
#read first as float so scientific notation is read, then convert to int

nevents = int(float(lines[nlines-1][0]) - float(lines[1][0])) + 1

print("nevents: ", nevents)

preeventnum = int(lines[1][0])
index = preeventnum-event1

minx = np.zeros(nevents)
maxx = np.zeros(nevents)
miny = np.zeros(nevents)
maxy = np.zeros(nevents)
mint = np.zeros(nevents)
maxt = np.zeros(nevents)

centerxstart = np.zeros(nevents)
centerystart = np.zeros(nevents)

mint[0] = int(re.findall(r'\d+',lines[1][2])[0])
maxt[0] = int(re.findall(r'\d+',lines[1][2])[0])

minx[0] = int(re.findall(r'\d+',lines[1][3])[0])
maxx[0] = int(re.findall(r'\d+',lines[1][3])[1])

miny[0] = int(re.findall(r'\d+',lines[1][4])[0])
maxy[0] = int(re.findall(r'\d+',lines[1][4])[0])

centerxstart[0] = int(re.findall(r'\d+',lines[1][6])[0])
centerystart[0] = int(re.findall(r'\d+',lines[1][6])[1])

for event in range(1,nlines):
	eventnum = int(float(lines[event][0]))
	if eventnum == preeventnum:
                maxt[index] = max(maxt[index],int(re.findall(r'\d+',lines[event][2])[0]))

		minx[index] = min(minx[index],int(re.findall(r'\d+',lines[event][3])[0]))
		maxx[index] = max(maxx[index],int(re.findall(r'\d+',lines[event][3])[1]))

                miny[index] = min(miny[index],int(re.findall(r'\d+',lines[event][4])[0]))
                maxy[index] = max(maxy[index],int(re.findall(r'\d+',lines[event][4])[1]))
	else:
	        index = eventnum-event1
		preeventnum = eventnum

		centerxstart[index] = int(re.findall(r'\d+',lines[event][6])[0])
		centerystart[index] = int(re.findall(r'\d+',lines[event][6])[1])
          
		mint[index] = int(re.findall(r'\d+',lines[event][2])[0])
                maxt[index] = int(re.findall(r'\d+',lines[event][2])[0])
 
		minx[index] = int(re.findall(r'\d+',lines[event][3])[0])
                maxx[index] = int(re.findall(r'\d+',lines[event][3])[1])

                miny[index] = int(re.findall(r'\d+',lines[event][4])[0])
                maxy[index] = int(re.findall(r'\d+',lines[event][4])[1])


nc_fid = Dataset(Dir + File1, 'r')



eventsin = ma.getdata(nc_fid.variables['value'][:,0,:,:])

ntimes = eventsin.shape[0]
nlats = eventsin.shape[1]
nlons = eventsin.shape[2]
maxevent = int(np.amax(eventsin))
nevents2 = maxevent - minevent + 1

if (nevents2 != nevents):
	print nevents2, nevents
	sys.exit("event numbers not equal in netcdf and text file")

pre_nc = Dataset(DirP + FileP1,'r')

precip_in = ma.getdata(pre_nc.variables['pcp'][:,:,:])


#precip_in = ma.getdata(Pr1_fid.variables['pcp'][:,:,:])
ntimespre = len(pre_nc.variables['time'])

print ntimes
print ntimespre

try:
    os.remove(DirO + FileO)
except OSError:
    pass

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('events', nevents)

results = ncfile.createVariable('span','f4',('events'),fill_value=-9999)
results2 = ncfile.createVariable('totalprecip','f4',('events'),fill_value=-9999)
timestart = ncfile.createVariable('tstart','f4',('events'),fill_value=-9999)
xcentstart = ncfile.createVariable('xcenstart','f4',('events'),fill_value=-9999)
ycentstart = ncfile.createVariable('ycenstart','f4',('events'),fill_value=-9999)
xmins = ncfile.createVariable('xmin','f4',('events'),fill_value=-9999)
xmaxs = ncfile.createVariable('xmax','f4',('events'),fill_value=-9999)
ymins = ncfile.createVariable('ymin','f4',('events'),fill_value=-9999)
ymaxs = ncfile.createVariable('ymax','f4',('events'),fill_value=-9999)

#with open(DirO + FileO, "a") as text_file:
#	text_file.write("Event number   Num gridcells   Total precip \n")

print "starting now"

for ievent in range(minevent,maxevent):
#		print ievent
        index = ievent-minevent

        tmin = max(0,mint[index]-1)
        tmax = min(ntimes,maxt[index]+1)

        ymin = max(0,miny[index]-1)
        ymax = min(nlats,maxy[index]+1)

        xmin = max(0,minx[index]-1)
        xmax = min(nlons,maxx[index]+1)

        if (xmin > xmax):
                print "wrapping around lons"
                xmin = 0
                xmax = nlons

        if (ymin > ymax):
                sys.exit("wrapping around lats - something is wrong!")


	eventsin_small = eventsin[tmin:tmax,ymin:ymax,xmin:xmax]
	data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == ievent))
	precip_in_small = precip_in[tmin:tmax,ymin:ymax,xmin:xmax]

#		text_file.write('%s , %s \n' % (np.sum([data_mask.mask]), np.sum([data_mask.mask * precip_in])))

	results[ievent-minevent] = (np.sum([data_mask_small.mask]))
	results2[ievent-minevent] = (np.sum([data_mask_small.mask * precip_in_small]))
	timestart[ievent-minevent] = mint[index]
	xmins[ievent-minevent] = minx[index]
        xmaxs[ievent-minevent] = maxx[index]
        ymins[ievent-minevent] = miny[index]
        ymaxs[ievent-minevent] = maxy[index]
        xcentstart[ievent-minevent] = centerxstart[index]
        ycentstart[ievent-minevent] = centerystart[index]



#results = results2
ncfile.close()

