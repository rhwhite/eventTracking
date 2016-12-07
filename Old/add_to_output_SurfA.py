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

test = 0

Data = "ERAI"

Version = "Standard"

filetimespan = "3hrly"

if filetimespan == "3hrly":
        if Data == "ERAI" or Data == "TRMM":
                # convert from mm/hour to mm/3 hours to get total rain over event for 3-hourly data
                mult = 3.0
        elif Data == "CESM":
                # convert from m/s to mm/3 hours to get total rain over event for 3-hourly data
                mult =  1000.0 * 60.0 * 60.0 * 3.0
if test == 1:
        chunksize = 50
else:
        chunksize = 2000

latlonaddsize = 2

startyr = 1998
endyr = 2014 


#DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
#FileP = 'TRMM_pcp_3hrly_nonan.2005.nc'
DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
FileP = 'TRMM_' + str(startyr) + '-' + str(endyr) + '_3B42_3hrly_nonan.nc'

minevent = 100000

if Data == "TRMM":
        startyr = 1998
	endyr = 2014
	Pfilestartyr = 1998
        Pfileendyr = 2014
        DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
        FileP = 'TRMM_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_3B42_3hrly_nonan.nc'
        Filegrid = "SurfaceArea.nc"
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
elif Data == "ERAI":
        Pfilestartyr = 1980
        Pfileendyr = 2015
        startyr = 1980
	endyr = 2014
	DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
        FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_preprocess.nc'
        Filegrid = "SurfaceArea.nc"

        Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
        File1 = 'ts_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
        TxtFileIn = 'ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
        FileO = 'Precip_Sizes_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
elif Data == "CESM":
        DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
        FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '.nc'
        Filegrid = "SurfaceArea.nc"
        Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) +'/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
        File1 = 'ts_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
        TxtFileIn = 'CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
        FileO = 'Precip_Sizes_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

# If output directory doesn't exist, create it!
if not os.path.exists(DirO):
        os.makedirs(DirO)

if test == 1:
        FileO = FileO + "_test.nc"

print FileO

# Open precip file (check it's there!)
print DirP + FileP
precipdata = xray.open_dataset(DirP + FileP)
# Open surface area grid:
print DirP + Filegrid
griddata = xray.open_dataset(DirP + Filegrid)
SurfA = griddata['SurfaceArea']


if test == 1:
        textdata = np.genfromtxt(Dir + TxtFileIn,skip_header = 1,usecols = (0,),max_rows = 100)
else:
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
                                startx[index] = listcenterx[0]
                                starty[index] = listcentery[0]
                                endx[index] = listcenterx[-1]
                                endy[index] = listcentery[-1]

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

                if test == 1:
                        if index == 10: break

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
        nevents = max(nevents,nevents2)
        print nevents2, nevents
        if test != 1:
                sys.exit("event numbers not equal in netcdf and text file")

if Data == "TRMM":
        precipin = precipdata['pcp']
elif Data == "ERAI":
        precipin = precipdata['tpnew']
elif Data == "CESM":
        # Convert from m/s to mm/day
        precipin = precipdata['PRECT']

ntimespre = len(precipdata['time'])

ncfile = Dataset(DirO + FileO, 'r+')
print DirO + FileO
try:
	Ogridbox = ncfile.createVariable('gridboxspanSA','f4',('events'),fill_value=-9999)
        setattr(Ogridbox, 'Extra Info', 'Total gridbox-span in m2')
except RuntimeError:
        print "gridboxspanSA already in file"
        Ogridbox = ncfile["gridboxspanSA"]

try:
	OtotalP = ncfile.createVariable('totalprecipSA','f4',('events'),fill_value=-9999)
        setattr(OtotalP, 'Extra Info', 'Total event precip in m3')
except RuntimeError:
        print "totalprecipSA already in file"
        OtotalP = ncfile["totalprecipSA"]

try:
        Oungridbox = ncfile.createVariable('uniquegridboxspanSA','f4',('events'),fill_value=-9999)
        setattr(OtotalP, 'Extra Info', 'Unique grid-box span in m3')
except RuntimeError:
        print "uniquegridboxspanSA already in file"
        Oungridbox = ncfile["uniquegridboxspanSA"]

print "starting now"

tminchunk = 0
tmaxchunk = 0

if test == 1:
        maxrun = minevent + 3
else:
        maxrun = maxevent + 1

#for ievent in range(minevent, minevent + 3): ### For testing!!!
for ievent in range(minevent,maxrun):

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
                # Here we include a multiplier to get from mm/hr (or mm/day) to mm/timestep (usually 3 hrs)

                precipchunk = mult * precipin.isel(time=slice(tminchunk,tmaxchunk)).values

        if (ievent % 5000 == 0):
                print "ievent: " + str(ievent)

        tminsel = max((tmin-tminchunk)-1,0)
        tmaxsel = min((tmax-tminchunk)+2,ntimes)

        eventsin_small = eventschunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
        data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == ievent))
        SurfA_small = SurfA[ymin:ymax,xmin:xmax].values
        precipin_small = precipchunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
        
	data_mask_max = np.amax(data_mask_small.mask,axis=0)

        Oungridbox[ievent-minevent] = (np.sum([data_mask_max*SurfA_small]))
        Ogridbox[ievent-minevent] = (np.sum([data_mask_small.mask*SurfA_small[None,:,:]]))
	#Factor of 0.001 to convert from mm to m
	# SurfS_small is in m2
	# So OtotalP is in m3
        OtotalP[ievent-minevent] = (0.001 * np.sum([data_mask_small.mask * precipin_small * SurfA_small[None,:,:]]))


#results = results2
ncfile.close()

