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

FileP = 'TRMM_' + str(startyr) + '-' + str(endyr) + '_3B42_3hrly.nc'

minevent = 100000
Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/'
DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Precip/'
File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_4Dobjects.nc'

TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_4Dobject_tree.txt'

with open(Dir + TxtFileIn,"r") as textFile:
	lines = [line.split('\t') for line in textFile]

nlines = len(lines)
print("nlimesn: ", nlines)
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

for event in range(1,12):    #nlines):
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

Pr1_fid = Dataset(DirP + FileP1,'r')
nc_pre = Dataset(DirP + FilePpre,'r')
nc_post = Dataset(DirP + FilePpost,'r')

#Really need to figure out why output data is shorter than input!!
precip_in = ma.getdata(Pr1_fid.variables['pcp'][0:ntimes,:,:])


#precip_in = ma.getdata(Pr1_fid.variables['pcp'][:,:,:])
ntimespre = len(nc_pre.variables['time'])
precip_pre = ma.getdata(nc_pre.variables['pcp'][ntimespre-1,:,:])
precip_post = ma.getdata(nc_post.variables['pcp'][0,:,:])

nc_fid.close()
nc_pre.close()
nc_post.close()

#Fill in missing precip data first, if possible. In same manner as for thresholding.
# May want to consider doing this in a more efficient manner? Maybe using iter command? Currently takes

"""
print "start looking at missing values"

mvs = np.where(precip_in < -1)

nvalues = len(mvs[0])
print nvalues
for nv in range(0,nvalues):
	if mvs[0][nv] == 0 or mvs[0][nv] == ntimes-1:
#can't do anything for missing values at the very beginning or very end of the timeseries
		print "can't deal with missing values sequentially in time, setting to zero"
		precip_in[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.0
	elif (precip_in[mvs[0][nv]-1,mvs[1][nv],mvs[2][nv]] > -1 and precip_in[mvs[0][nv]+1,mvs[1][nv],mvs[2][nv]] > -1):
        	precip_in[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.5 * precip_in[mvs[0][nv]-1,mvs[1][nv],mvs[2][nv]] + \
			0.5 * precip_in[mvs[0][nv]+1,mvs[1][nv],mvs[2][nv]]
	else:
		print "can't deal  missing values sequentially in time, setting to zero"
		precip_in[mvs[0][nv],mvs[1][nv],mvs[2][nv]] = 0.0 

#	if precip_in[mvs[0][nv],mvs[1][nv],mvs[2][nv]] != 0.0:
#		print precip_in[mvs[0][nv],mvs[1][nv],mvs[2][nv]]

		
mvs2 = np.where(precip_in < -1)
nvalues2 = len(mvs2[0])
if nvalues2 != 0:
	sys.err("error with filling missing data; missing data remains")
	

print 'finished filling in missing data'

"""

FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '.nc'
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

for ievent in range(minevent,minevent+10):    #maxevent):
#		print ievent
        index = ievent-minevent

	eventsin_small = eventsin[mint[index]:maxt[index]+1,miny[index]:maxy[index]+1,minx[index]:maxx[index]+1]
	data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == ievent))
	precip_in_small = precip_in[mint[index]:maxt[index]+1,miny[index]:maxy[index]+1,minx[index]:maxx[index]+1]

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

