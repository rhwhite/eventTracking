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

#Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Old/'
#DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Old/Precip2/'
#File1 = 'ts_TRMMtest_5th_2005_4Dobjects.nc'
#FileO = 'Precip_Sizes_2005_test_new_process_mm.nc'
#TxtFileIn = 'TRMMtest_5th_2005_4Dobject_tree.txt'

Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/final/'
DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Precip/'
File1 = 'ts_TRMM' + str(startyr) + '-' + str(endyr) + '_final_4Dobjects.nc'
TxtFileIn = 'TRMM' + str(startyr) + '-' + str(endyr) + '_final_4Dobject_tree.txt'
FileO = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + 'tspan.nc'

# for testing

with open(Dir + TxtFileIn,"r") as textFile:
        lines = [line.split('\t') for line in textFile]

textFile.close()

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

centerx = np.zeros(nevents)
centery = np.zeros(nevents)
meant = np.zeros(nevents)

mint[0] = int(re.findall(r'\d+',lines[1][2])[0])
maxt[0] = int(re.findall(r'\d+',lines[1][2])[0])

minx[0] = int(re.findall(r'\d+',lines[1][3])[0])
maxx[0] = int(re.findall(r'\d+',lines[1][3])[1])

miny[0] = int(re.findall(r'\d+',lines[1][4])[0])
maxy[0] = int(re.findall(r'\d+',lines[1][4])[0])

listcenterx = int(re.findall(r'\d+',lines[1][6])[0])
listcentery = int(re.findall(r'\d+',lines[1][6])[1])
listt = int(re.findall(r'\d+',lines[1][2])[0])

for event in range(1,nlines):
        eventnum = int(float(lines[event][0]))
        if eventnum == preeventnum:
                maxt[index] = max(maxt[index],int(re.findall(r'\d+',lines[event][2])[0]))

                minx[index] = min(minx[index],int(re.findall(r'\d+',lines[event][3])[0]))
                maxx[index] = max(maxx[index],int(re.findall(r'\d+',lines[event][3])[1]))

                miny[index] = min(miny[index],int(re.findall(r'\d+',lines[event][4])[0]))
                maxy[index] = max(maxy[index],int(re.findall(r'\d+',lines[event][4])[1]))
		listcenterx.append(int(re.findall(r'\d+',lines[event][6])[0]))
                listcenterx.append(int(re.findall(r'\d+',lines[event][6])[1]))
                listt.append(int(re.findall(r'\d+',lines[event][2])[0]))
        else:

		centerx[index] = np.mean(listcenterx)
		centery[index] = np.mean(listcentery)
		meant[index] = np.mean(listt)

                index = eventnum-event1
                preeventnum = eventnum

                mint[index] = int(re.findall(r'\d+',lines[event][2])[0])
                maxt[index] = int(re.findall(r'\d+',lines[event][2])[0])

                minx[index] = int(re.findall(r'\d+',lines[event][3])[0])
                maxx[index] = int(re.findall(r'\d+',lines[event][3])[1])

                miny[index] = int(re.findall(r'\d+',lines[event][4])[0])
                maxy[index] = int(re.findall(r'\d+',lines[event][4])[1])

		listcenterx = int(re.findall(r'\d+',lines[event][6])[0])
		listcentery = int(re.findall(r'\d+',lines[event][6])[1])
		listt = int(re.findall(r'\d+',lines[event][2])[0])


print 'centerx'
print centerx[0:10]
print 'centery'
print centery[0:10]
print 'meant'
print meant[0:10]


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


#precipdata = xray.open_dataset(DirP + FileP)

#precipin = precipdata['pcp']
#ntimespre = len(precipdata['time'])

#xray.Dataset.close(eventsdata)
#xray.Dataset.close(precipdata)

try:
    os.remove(DirO + FileO)
except OSError:
    pass

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('events', nevents)

#results = ncfile.createVariable('gridboxspan','f4',('events'),fill_value=-9999)
#results2 = ncfile.createVariable('totalprecip','f4',('events'),fill_value=-9999)
timespan = ncfile.createVariable('timespan','f4',('events'),fill_value=-9999)
timestart = ncfile.createVariable('tstart','f4',('events'),fill_value=-9999)
timemean = ncfile.createVariable('tmean','f4',('events'),fill_value=-9999)
#xcentstart = ncfile.createVariable('ycenterstart','f4',('events'),fill_value=-9999)
#ycentstart = ncfile.createVariable('xcenterstart','f4',('events'),fill_value=-9999)
xcentmean = ncfile.createVariable('ycentermean','f4',('events'),fill_value=-9999)
ycentmean = ncfile.createVariable('xcentermean','f4',('events'),fill_value=-9999)
xmins = ncfile.createVariable('xmin','f4',('events'),fill_value=-9999)
xmaxs = ncfile.createVariable('xmax','f4',('events'),fill_value=-9999)
ymins = ncfile.createVariable('ymin','f4',('events'),fill_value=-9999)
ymaxs = ncfile.createVariable('ymax','f4',('events'),fill_value=-9999)

#with open(DirO + FileO, "a") as text_file:
#	text_file.write("Event number   Num gridcells   Total precip \n")

print "starting now"

tminchunk = 0
tmaxchunk = 0

for ievent in range(minevent,maxevent + 1):
        index = ievent-minevent
	tmin = max(0,mint[index]-1)
	tmax = min(ntimes,maxt[index]+1)

	ymin = max(0,miny[index]-latlonaddsize)
	ymax = min(nlats,maxy[index]+latlonaddsize)

	xmin = max(0,minx[index]-latlonaddsize)
	xmax = min(nlons,maxx[index]+latlonaddsize)

	if (xmin > xmax):
		print "wrapping around lons"
		xmin = 0
		xmax = nlons

        if (ymin > ymax):
		sys.exit("wrapping around lats - something is wrong!")

#	if (tmax > tmaxchunk):
#		print "tmax is " + str(tmax) + " whilst tmaxchunk is " + str(tmaxchunk)
		tminchunk = max(tmin - 10,0)
#		tmaxchunk = min(tmax + chunksize + 1,ntimes)
#		print "tminchunk is now" + str(tminchunk) + "whilst tmin is " + str(tmin)
#		print "tmaxchunk is " + str(tmaxchunk) + "whilst tmax is " + str(tmax)
#		try:
#			del eventschunk
#		except NameError:
#			pass
#		try:
#			del precipchunk
#		except NameError:
#			pass
#		eventschunk = eventsin.isel(time=slice(tminchunk,tmaxchunk)).values
#		precipchunk = mult * precipin.isel(time=slice(tminchunk,tmaxchunk)).values
		
	if (ievent % 10000 == 0):
		print "ievent: " + str(ievent)

#	tminsel = max((tmin-tminchunk)-1,0) 
#	tmaxsel = min((tmax-tminchunk)+2,ntimes)

#	eventsin_small = eventschunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
#	data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == ievent))

#        precipin_small = precipchunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
#		text_file.write('%s , %s \n' % (np.sum([data_mask.mask]), np.sum([data_mask.mask * precipin])))

#	results[ievent-minevent] = (np.sum([data_mask_small.mask]))
#	results2[ievent-minevent] = (np.sum([data_mask_small.mask * precipin_small]))
	timestart[ievent-minevent] = mint[index]
	xmins[ievent-minevent] = minx[index]
        xmaxs[ievent-minevent] = maxx[index]
        ymins[ievent-minevent] = miny[index]
        ymaxs[ievent-minevent] = maxy[index]
#        xcentstart[ievent-minevent] = centerxstart[index]
#        ycentstart[ievent-minevent] = centerystart[index]

	timespan[ievent-minevent] = 1 + maxt[index]-mint[index]
	timemean[ievent-minevent] = meant[index]
	xcentmean[ievent-minevent] = centerx[index]
	ycentmean[ievent-minevent] = centery[index]


#results = results2
ncfile.close()

