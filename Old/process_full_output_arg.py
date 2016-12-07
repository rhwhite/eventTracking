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
import argparse

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--test',type=int,nargs='?',default=0,help='if 1, run quick test version of code')
parser.add_argument('--filetspan',type=str,nargs='?',default='3hrly',help='string for file time resolution, 3hrly etc')
parser.add_argument('--speedtspan',type=int,nargs='?',default=4,help='how many timespans to avg speed over')

args = parser.parse_args()

print args.Data

if args.Data[0] not in ['TRMM','ERAI','CESM']:
        exit("incorrect Data option " + str(args.Data[0]) + " must be TRMM, ERAI, or CESM")

Data = args.Data[0]
Version = args.Version[0]

filetimespan = args.filetspan

if filetimespan == "3hrly":
	if Data in ["ERAI","TRMM","TRMM_ERAIgd"]:
		# convert from mm/hour to mm/3 hours to get total rain over event for 3-hourly data
		mult = 3.0
	elif Data == "CESM":
		# convert from m/s to mm/3 hours to get total rain over event for 3-hourly data
		mult =  1000.0 * 60.0 * 60.0 * 3.0
	else:
		sys.error(Data + " not defined")
		
if args.test == 1:
	chunksize = 50
else:
	chunksize = 2000

latlonaddsize = 2

startyr = 1998
endyr = 2014
minevent = 100000


if Data == "TRMM":
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
elif Data == "TRMM_ERAIgd":
        Pfilestartyr = 1998 # Don't change - tied to file names!
        Pfileendyr = 2014
        filetimespan = "3hrly"

        DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
        FileP = "regrid2ERAI_TRMM_3B42_" + str(Pfilestartyr) + '-' + str(Pfileendyr) + ".nc"

	Filegrid = Data + "_SurfaceArea.nc"
        Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/TRMM_output/ERAIgd/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/TRMM_output/ERAIgd/Precip/'
        File1 = 'ts_TRMM_ERAIgd_' + str(startyr) + '-' + str(endyr) + '_4Dobjects.nc'
        TxtFileIn = 'testTRMM_ERAIgd_' + str(startyr) + '-' + str(endyr) + '_4Dobject_tree.txt'
        FileO = 'Precip_Sizes_TRMM_ERAIgd_' + str(startyr) + '-' + str(endyr) + '.nc'
elif Data == "ERAI":
	startyr = 1980
	endyr = 2015
	Pfilestartyr = 1980
	Pfileendyr = 2015
        DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
        FileP = 'ERAI_Totalprecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '_preprocess.nc'
	Filegrid = "SurfaceArea.nc"
	
	Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/'
	DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	File1 = 'ts_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
	TxtFileIn = 'ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
	FileO = 'Precip_Sizes_ERAI' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
elif Data == "CESM":
	Pfilestartyr = 1979
	Pfileendyr = 2012
        DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
        FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(Pfilestartyr) + '-' + str(Pfileendyr) + '.nc'
	Filegrid = "SurfaceArea.nc"
        Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) +'/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
        File1 = 'ts_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc'
        TxtFileIn = 'CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobject_tree.txt'
        FileO = 'Precip_Sizes_CESM' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

else:
	sys.error("unknown Data type")

# If output directory doesn't exist, create it!
if not os.path.exists(DirO):
        os.makedirs(DirO)

if args.test == 1:
        FileO = FileO + "_test.nc"

print FileO

print DirP + FileP
precipdata = xray.open_dataset(DirP + FileP)
# Open surface area grid:
print DirP + Filegrid
griddata = xray.open_dataset(DirP + Filegrid)
SurfA = griddata['SurfaceArea']


if args.test == 1:
	textdata = np.genfromtxt(Dir + TxtFileIn,skip_header = 1,usecols = (0,),max_rows = 100)
else:
	textdata = np.loadtxt(Dir + TxtFileIn,skiprows = 1, usecols = (0,))

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
maxzonalspeed = np.zeros(nevents)

#centerxstart = np.zeros(nevents)
#centerystart = np.zeros(nevents)

centerx = np.zeros(nevents)
centery = np.zeros(nevents)
meant = np.zeros(nevents)

preeventnum = -100

with open(Dir + TxtFileIn,"r") as textFile:
	next(textFile)	#Skip header line

	for lines in textFile:	#loop through all lines
		line = lines.split('\t')
	        eventnum = int(float(line[0]))
        	if (eventnum % 5000000 == 0):
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
                        eventcount += 1

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
                                if eventcount == 1:
                                        maxzonalspeed[index] = 0.0      # Stationary event by definition of one timestep
                                elif eventcount <= speedtspan:
                                        diffs = listcenterx[eventcount-1] - listcenterx[0]
                                        maxzonalspeed[index] = np.cos(np.radians(lats[listcentery[int(eventcount/2 -1)]])) * R * 2.0 * math.pi * diffs/(nlons * 3.0 * eventcount * 60.0 * 60.0)   # speed in m/s
                                elif eventcount > speedtspan:
                                        diffs = np.array([x - listcenterx[i-speedtspan] for i,x in enumerate(listcenterx)][speedtspan:])
                                        coslat = np.cos(np.radians(lats[listcentery[speedtspan:]]))
                                        maxindex = np.argmax(abs(diffs)*coslat)
                                        maxzonalspeed[index] = np.cos(np.radians(lats[listcentery[maxindex]])) * R * 2.0 * math.pi * diffs[maxindex]/(nlons * 3.0 * speedtspan * 60.0 * 60.0)   # speed in m/s
                                else:
                                        sys.exit("eventcount is less than 0?")

							
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
			eventcount=1

		if args.test == 1:
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
	#if args.test != 1:
	#	sys.exit("event numbers not equal in netcdf and text file")

if Data in ["TRMM","TRMM_ERAIgd"]:
	precipin = precipdata['pcp']
elif Data == "ERAI":
        precipin = precipdata['tpnew']
elif Data == "CESM":
        # No conversion here as put into mult instead
	precipin = precipdata['PRECT']

ntimespre = len(precipdata['time'])

#xray.Dataset.close(eventsdata)
#xray.Dataset.close(precipdata)

try:
    os.remove(DirO + FileO)
except OSError:
    pass


ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('events', nevents)

Ogridboxspan = ncfile.createVariable('gridboxspan','f4',('events'),fill_value=-9999)
Ototalprecip = ncfile.createVariable('totalprecip','f4',('events'),fill_value=-9999)
Oungridboxspan = ncfile.createVariable('uniquegridboxspan','f4',('events'),fill_value=-9999)
OgridboxspanSA = ncfile.createVariable('gridboxspanSA','f4',('events'),fill_value=-9999)
OtotalprecipSA = ncfile.createVariable('totalprecipSA','f4',('events'),fill_value=-9999)
OungridboxspanSA = ncfile.createVariable('uniquegridboxspanSA','f4',('events'),fill_value=-9999)
Otimespan = ncfile.createVariable('timespan','f4',('events'),fill_value=-9999)
Otimestart = ncfile.createVariable('tstart','f4',('events'),fill_value=-9999)
Otimemean = ncfile.createVariable('tmean','f4',('events'),fill_value=-9999)
Oxstartmean = ncfile.createVariable('xcenterstart','f4',('events'),fill_value=-9999)
Oxendmean = ncfile.createVariable('xcenterend','f4',('events'),fill_value=-9999)
Oystartmean = ncfile.createVariable('ycenterstart','f4',('events'),fill_value=-9999)
Oyendmean = ncfile.createVariable('ycenterend','f4',('events'),fill_value=-9999)
Oxcentmean = ncfile.createVariable('xcentermean','f4',('events'),fill_value=-9999)
Oycentmean = ncfile.createVariable('ycentermean','f4',('events'),fill_value=-9999)
Oxmins = ncfile.createVariable('xmin','f4',('events'),fill_value=-9999)
Oxmaxs = ncfile.createVariable('xmax','f4',('events'),fill_value=-9999)
Oymins = ncfile.createVariable('ymin','f4',('events'),fill_value=-9999)
Oymaxs = ncfile.createVariable('ymax','f4',('events'),fill_value=-9999)

varname = 'xmaxspeed_' + str(speedtspan) + 'ts'
Omaxzonalspeed = ncfile.createVariable(varname,'f4',('events'),fill_value=-9999)
setattr(Omaxzonalspeed, 'Extra Info', 'maximum zonal speed during lifetime for 3 * ' + str(speedtspan) + 'hours, for non-zero values; minimum non-zero value is ~2.5m/s at equator based on spatial and temporal resolution')

setattr(gridboxspanSA,'Extra Info','spatial span in m2')
setattr(totalprecipSA,'Extra Info','total precip in m3')
setattr(ungridboxspanSA,'Extra Info','unique spatial span in m2')


#with open(DirO + FileO, "a") as text_file:
#	text_file.write("Event number   Num gridcells   Total precip \n")

print "starting now"

tminchunk = 0
tmaxchunk = 0

if args.test == 1:
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
		
	if (ievent % 5000000 == 0):
		print "ievent: " + str(ievent)

	tminsel = max((tmin-tminchunk)-1,0) 
	tmaxsel = min((tmax-tminchunk)+2,ntimes)

	eventsin_small = eventschunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
	data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == ievent))
	SurfA_small = SurfA[ymin:ymax,xmin:xmax].values
        precipin_small = precipchunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
	
        data_mask_max = np.amax(data_mask_small.mask,axis=0)

	ungridboxspan[ievent-minevent] = (np.sum([data_mask_max]))
        gridboxspan[ievent-minevent] = (np.sum([data_mask_small.mask]))
        totalprecip[ievent-minevent] = (0.001 * np.sum([data_mask_small.mask * precipin_small]))

	# SurfA is in m2
        ungridboxspanSA[ievent-minevent] = (np.sum([data_mask_max*SurfA_small]))
	gridboxspanSA[ievent-minevent] = (np.sum([data_mask_small.mask*SurfA_small[None,:,:]]))
        # Factor of 0.001 to convert from mm to m
        # SurfS_small is in m2
        # So OtotalP is in m3
	totalprecipSA[ievent-minevent] = (0.001 * np.sum([data_mask_small.mask * precipin_small * SurfA_small[None,:,:]]))


timespan[:] = 1 + maxt[:]-mint[:]

Otimestart[:] = mint[:]
Oxmins[:] = minx[:]
Oxmaxs[:] = maxx[:]
Oymins[:] = miny[:]
Oymaxs[:] = maxy[:]

Oxstartmean[:] = startx[:]
Oxendmean[:] = endx[:]
Oystartmean[:] = starty[:]
Oyendmean[:] = endy[:]

Otimemean[:] = meant[:]
Oxcentmean[:] = centerx[:]
Oycentmean[:] = centery[:]

ncfile.close()

