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

Version = "5th_from48"

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
	FileO = 'Precip_Sizes' + str(startyr) + '-' + str(endyr) + '_Standard.nc'
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
			maxt[index] = max(maxt[index],int(re.findall(r'\d+',line[2])[0]))

			minx[index] = min(minx[index],int(re.findall(r'\d+',line[3])[0]))
			maxx[index] = max(maxx[index],int(re.findall(r'\d+',line[3])[1]))

			miny[index] = min(miny[index],int(re.findall(r'\d+',line[4])[0]))
			maxy[index] = max(maxy[index],int(re.findall(r'\d+',line[4])[1]))
			listt.append(int(re.findall(r'\d+',line[2])[0]))
			listcenterx.append(int(re.findall(r'\d+',line[6])[0]))
			listcentery.append(int(re.findall(r'\d+',line[6])[1]))
		else:
			if (preeventnum > 0): #then it's not the very first event
				#new event: take mean of center x, center y and t and put into array before updating index
				centerx[index] = np.mean(listcenterx)
				centery[index] = np.mean(listcentery)
				meant[index] = np.mean(listt)
			
			#refresh list
			#update index for new event
			index = eventnum-event1
			preeventnum = eventnum

	#               centerxstart[index] = int(re.findall(r'\d+',line[6])[0])
	#               centerystart[index] = int(re.findall(r'\d+',line[6])[1])

			mint[index] = int(re.findall(r'\d+',line[2])[0])
			maxt[index] = int(re.findall(r'\d+',line[2])[0])

			minx[index] = int(re.findall(r'\d+',line[3])[0])
			maxx[index] = int(re.findall(r'\d+',line[3])[1])

			miny[index] = int(re.findall(r'\d+',line[4])[0])
			maxy[index] = int(re.findall(r'\d+',line[4])[1])

			listt = [int(re.findall(r'\d+',line[2])[0])]
			listcenterx = [int(re.findall(r'\d+',line[6])[0])]
			listcentery = [int(re.findall(r'\d+',line[6])[1])]


textFile.close()

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


precipdata = xray.open_dataset(DirP + FileP)

precipin = precipdata['pcp']
ntimespre = len(precipdata['time'])

#xray.Dataset.close(eventsdata)
#xray.Dataset.close(precipdata)

try:
    os.remove(DirO + FileO)
except OSError:
    pass

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('events', nevents)

results = ncfile.createVariable('gridboxspan','f4',('events'),fill_value=-9999)
results2 = ncfile.createVariable('totalprecip','f4',('events'),fill_value=-9999)
results3 = ncfile.createVariable('uniquegridboxspan','f4',('events'),fill_value=-9999)
timespan = ncfile.createVariable('timespan','f4',('events'),fill_value=-9999)
timestart = ncfile.createVariable('tstart','f4',('events'),fill_value=-9999)
timemean = ncfile.createVariable('tmean','f4',('events'),fill_value=-9999)
#xcentstart = ncfile.createVariable('ycenterstart','f4',('events'),fill_value=-9999)
#ycentstart = ncfile.createVariable('xcenterstart','f4',('events'),fill_value=-9999)
xcentmean = ncfile.createVariable('xcentermean','f4',('events'),fill_value=-9999)
ycentmean = ncfile.createVariable('ycentermean','f4',('events'),fill_value=-9999)
xmins = ncfile.createVariable('xmin','f4',('events'),fill_value=-9999)
xmaxs = ncfile.createVariable('xmax','f4',('events'),fill_value=-9999)
ymins = ncfile.createVariable('ymin','f4',('events'),fill_value=-9999)
ymaxs = ncfile.createVariable('ymax','f4',('events'),fill_value=-9999)

#with open(DirO + FileO, "a") as text_file:
#	text_file.write("Event number   Num gridcells   Total precip \n")

print "starting now"

tminchunk = 0
tmaxchunk = 0

#for ievent in range(minevent, minevent + 3): ### For testing!!!
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

	if (tmax > tmaxchunk):
		print "tmax is " + str(tmax) + " whilst tmaxchunk is " + str(tmaxchunk)
		tminchunk = max(tmin - 10,0)
		tmaxchunk = min(tmax + chunksize + 1,ntimes)
		print "tminchunk is now " + str(tminchunk) + "whilst tmin is " + str(tmin)
		print "tmaxchunk is now " + str(tmaxchunk) + "whilst tmax is " + str(tmax)
		eventschunk = eventsin.isel(time=slice(tminchunk,tmaxchunk)).values
		precipchunk = mult * precipin.isel(time=slice(tminchunk,tmaxchunk)).values
		
	if (ievent % 5000 == 0):
		print "ievent: " + str(ievent)

	tminsel = max((tmin-tminchunk)-1,0) 
	tmaxsel = min((tmax-tminchunk)+2,ntimes)

	eventsin_small = eventschunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
	data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == ievent))

        precipin_small = precipchunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
#		text_file.write('%s , %s \n' % (np.sum([data_mask.mask]), np.sum([data_mask.mask * precipin])))

	
        data_mask_max = np.amax(data_mask_small.mask,axis=0)
#               text_file.write('%s , %s \n' % (np.sum([data_mask.mask]), np.sum([data_mask.mask * precipin])))

        results3[ievent-minevent] = (np.sum([data_mask_max]))
	results[ievent-minevent] = (np.sum([data_mask_small.mask]))
	results2[ievent-minevent] = (np.sum([data_mask_small.mask * precipin_small]))
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

