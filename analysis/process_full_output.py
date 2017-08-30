"""
Code to process raw output from the Forward-in-Time event tracking code.
Currently set up for TRMM, ERA-interim, CESM, ERA20C and GPCP data

Author: Rachel White, rachel.white@cantab.net

Created: Jan 2016

Example use: python process_full_output.py --Data TRMM --Version
Standard --startyr 1998 --endyr 2014

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
import math
import argparse
import resource
import time
import warnings

from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.readwrite import getunitsdesc

np.seterr(all='warn')
warnings.simplefilter("error",UserWarning)

rsrcV = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit starts as  :', soft
print 'Hard limit starts as :', hard
resource.setrlimit(rsrcV, (120000000000, hard)) #limit memory usage
#                          137438953472
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit changed to :', soft


parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, ERA20C, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--test',type=int,nargs='?',default=0,help='if 1, run quick test version of code')
parser.add_argument('--speedtspan',type=int,nargs='?',default=[4],help='how many timespans to avg speed over')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')

args = parser.parse_args()

print "here's what I have as arguments: ", args

speedtspan = args.speedtspan[0]
print speedtspan

if args.Data[0] not in ['TRMM','ERAI','CESM','TRMMERAIgd','ERA20C','GPCP']:
    exit("incorrect Data option " + str(args.Data[0]) + 
         " must be TRMM, TRMMERAIgd,ERAI,ERA20C, CESM, or GPCP")

Data = args.Data[0]
Version = args.Version[0]

startyr = args.startyr[0]
endyr = args.endyr[0]


if Data in ["ERAI","TRMM","TRMMERAIgd","ERA20C","CESM"]:
    filetimespan = "3hrly"
    tmult = 3   # multiplier for timespans to get to hours
    if Data in ["ERAI","TRMM","TRMMERAIgd","ERA20C"]:
        # convert from mm/hour to mm/3 hours to get total rain over event for
        # 3-hourly data
        mult = 3.0

    elif Data == "CESM":
        # convert from m/s to mm/3 hours to get total rain over event for 3-hourly
        # data
        mult =  1000.0 * 60.0 * 60.0 * 3.0

elif Data == 'GPCP':
    filetimespan = "daily"
    tmult = 24   # multiplier for timespans to get to hours
    # no need to convert: data in mm/day, and resolution is daily
    mult = 1.0
else:
    sys.error(Data + " not defined")


if args.test == 1:
    chunksize = 50
else:
    chunksize = 20000

latlonaddsize = 2

minevent = 100000
R = 6371000     # radius of Earth in m


def geteventdata(ievent):
    maxt[index] = eventtimes[iline-1]
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


Filegrid = Data + "_SurfaceArea.nc"

if Data == "TRMM":
    DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
    FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'
elif Data == "TRMMERAIgd":
    DirP = '/home/disk/eos4/rachel/Obs/TRMM/'
    FileP = "regrid2ERAI_TRMM_3B42_1998-2014.nc"
elif Data == "ERAI":
    DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
    FileP = 'ERAI_Totalprecip_1980-2015_preprocess.nc'
elif Data == "ERA20C":
    DirP = '/home/disk/eos4/rachel/Obs/ERA_20C/'
    FileP = 'ERA_20C_Totalprecip_' + str(startyr) + '-' + str(endyr) + '_preprocess.nc'
elif Data == "CESM":
    DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
    FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'
elif Data == "GPCP":
    DirP = '/home/disk/eos4/rachel/Obs/GPCP/Daily/'
    FileP = 'GPCP_1DD_v1.2_199610-201510.nc'


else:
    sys.error("unknown Data type")

if Data == 'TRMM' and Version in ["Standard","6th_from6"]:
    Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW/' + Data + '_output/' + Version + str(startyr) + '/raw/'
else:
    Dir = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' + Version + str(startyr) + '/raw/'

File1 = 'ts_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_4Dobjects.nc'
TxtFileIn = Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_4Dobject_tree.txt'

DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' + Version + str(startyr) + '/proc/'
FileO = 'All_Precip_' + str(startyr) + '-' + str(endyr) + '_' + Data + '_' + Version + '.nc'

# If output directory doesn't exist, create it!
if not os.path.exists(DirO):
    os.makedirs(DirO)

if args.test == 1:
    FileO = FileO + "_test.nc"

print FileO

precipdata = xrayOpen(DirP + FileP)
# Open surface area grid:
griddata = xrayOpen(DirP + Filegrid)
SurfA = griddata['SurfaceArea']

# Open precip file for latitudes
if Data in ["CESM",'GPCP']:
    lons = precipdata['lon'].values
    lats = precipdata['lat'].values
else:
    lons = precipdata['longitude'].values
    lats = precipdata['latitude'].values

nlats = lats.shape[0]
nlons = lons.shape[0]

with open(Dir + TxtFileIn, 'rb') as fh:
    header = str.split(str(next(fh).decode()))
    try:
        first = str.split(str(next(fh).decode()))
        fh.seek(-2048,2)
        last = str.split(str(fh.readlines()[-1].decode()))
    except StopIteration:
        print "found end of file!"
event1 = int(first[0])
nevents = int(last[0]) - event1 + 1

if args.test == 1:
    nevents = 500
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

print first
timestart = time.time()

maxlines = 3000000
skips = 1
eventnum = 0

while eventnum < nevents:
    for attempt in range(10):
        try:
            print "skips,",skips
            try:
                textdata1 = np.genfromtxt(Dir + TxtFileIn,dtype='int',skip_header = skips, max_rows = maxlines,usecols = (0,2))
                textdata2 = np.genfromtxt(Dir + TxtFileIn,dtype='string',skip_header = skips, max_rows = maxlines,usecols = (3,4,6))
            except (UserWarning,StopIteration):
                print "We've got to the end of the file!"
                eventnum += 1   # to get out of outer loop!
                break

            nlines = textdata1.shape[0]
            print "nlines ", nlines
            # Check that the next read will have more than one line
            try:
                test = np.genfromtxt(Dir + TxtFileIn,dtype='int',skip_header = skips+nlines, max_rows = maxlines,usecols = (0))
            except (UserWarning,StopIteration):
                print "couldn't read in test; must be end of file"
            else:
                if test.shape[0] == 1:
                    print test.shape
                    maxlines = maxlines + 5
                    textdata1 = np.genfromtxt(Dir + TxtFileIn,dtype='int',skip_header = skips, max_rows = maxlines,usecols = (0,2))
                    textdata2 = np.genfromtxt(Dir + TxtFileIn,dtype='string',skip_header = skips, max_rows = maxlines,usecols = (3,4,6))

                    nlines = textdata1.shape[0]

            print "textdata1 shape: ", textdata1.shape
            eventnums = textdata1[:,0] - event1
            eventtimes = textdata1[:,1]

            timeend = time.time()
            #print "number of events: " + str(eventnums[nlines-1])
            #print "time to read " + str(timeend - timestart)

            xs = np.zeros([nlines,2],int)
            ys = np.zeros([nlines,2],int)
            cxy = np.zeros([nlines,3],int)

            # Get int data out of strings:
            preeventnum = -100
            print "eventnum first in chunk", eventnums[0]

            # Start loop through every line in this chunk
            for iline in range(0,nlines):
                eventnum = eventnums[iline]
                xs[iline,:] = (re.findall(r'\d+',textdata2[iline,0]))
                ys[iline,:] = (re.findall(r'\d+',textdata2[iline,1]))
                cxy[iline,:] = (re.findall(r'\d+',textdata2[iline,2]))

                # if line is a continuation of a previous event
                if eventnum == preeventnum:
                    minx[index] = min(minx[index],xs[iline,0])
                    maxx[index] = max(maxx[index],xs[iline,1])
                    miny[index] = min(miny[index],ys[iline,0])
                    maxy[index] = max(maxy[index],ys[iline,1])

                    listt.append(eventtimes[iline])

                    listcenterx.append(cxy[iline,0])
                    listcentery.append(cxy[iline,1])
                    eventcount += 1
                    if iline == nlines -1:  #i.e. we're at the end of our chunk
                        #print "in the middle of event ", eventnum

                        if nlines < maxlines:   # then we reached the end of the file 
                            print "end of file: ", nlines, maxlines
                            skips += nlines
                            print "skips now ",skips
                            eventnum += 1   # to get out of outer loop!
                            geteventdata(index)

                            break   # break out of loop
                        else:
                            print "eventnum last in chunk ", eventnum
                            print "eventcount", eventcount
                            print "skips before", skips
                            skips += nlines-eventcount      # Skip back to beginning of this event
                            print 'skips here', skips
                        # End if for continued event as last in chunk

                else:       # if this line is a new event
                    if (preeventnum >= 0): #then it's not the very first event
                        geteventdata(index)

                    #update index for new event
                    index = eventnum

                    if index >= nevents:     # then we're testing and have got
                                            # all the events we need.
                        break
                    else:
                        preeventnum = eventnum

                        mint[index] = eventtimes[iline]
                        minx[index] = xs[iline,0]
                        maxx[index] = xs[iline,1]
                        miny[index] = ys[iline,0]
                        maxy[index] = ys[iline,1]

                        listt = [(eventtimes[iline])]

                        listcenterx = [(cxy[iline,0])]
                        listcentery = [(cxy[iline,1])]

                        eventcount=1

                    # If last event was single, now we've read in the line
                    # about it, get data, and then add skips to the filecount
                    # for the next chunk of the file
                    if iline == nlines-1:   # i.e. we're at the end of our chunk
                        print "last event was single"
                        print "eventnum last in chunk ", eventnum
                        geteventdata(index)

                        if nlines < maxlines:   # then we reached the end of the file 
                            print "add to skips so we read end of file", skips, nlines
                            skips += nlines
                            print "skips now", skips
                        else:
                            skips += nlines-1 # go back to start of this event.
                # End if for checking if event is a single line or not
            # End loop over every line in chunk
        except MemoryError:
            print "Memory error, trying to chunk"
            maxlines = int(maxlines/2)
            print "new maxlines: ", maxlines
            try:
                del(textdata1,textdata2,eventnums,eventtimes,xs,ys,cxy)
            except NameError:
                continue
            continue
        break

#print lines[0]
timeend = time.time()
print "time to read and process = " + str(timeend - timestart) + " seconds"

#print(np.genfromtxt(Dir + TxtFileIn,dtype='int',skip_header = 1, max_rows = 6,usecols = (0,2)))
                        #textdata = np.loadtxt(Dir + TxtFileIn,skiprows = 1, usecols = (0,2,3,4,6))
#print(np.genfromtxt(Dir + TxtFileIn,dtype='string',skip_header = 1, max_rows = 6,usecols = (3,4,6)))


timeend = time.time()
print "this took " + str(timeend - timestart) + " seconds"

eventsdata = xray.open_dataset(Dir + File1)

eventsin=eventsdata['value']

eventsin = eventsin.squeeze(dim="z")

ntimes = eventsin.shape[0]
nlats = eventsin.shape[1]
nlons = eventsin.shape[2]
maxevent = int(np.amax(eventsin[ntimes-10:ntimes,:,:]))

nevents2 = maxevent - minevent + 1

print maxevent
if (nevents2 != nevents):
    print nevents2, nevents
#   if args.test != 1:
#       sys.exit("event numbers not equal in netcdf and text file")

if Data in ["TRMM","TRMMERAIgd"]:
    precipin = precipdata['pcp']
elif Data in ["ERAI","ERA20C"]:
    precipin = precipdata['tpnew'].sel(time=slice(str(startyr),str(endyr)))
elif Data == "CESM":
    # No conversion here as put into mult instead
    precipin = precipdata['PRECT'].sel(time=slice(str(startyr),str(endyr)))
elif Data == 'GPCP':
    precipin = precipdata['PREC'].sel(time=slice(str(startyr),str(endyr)))

ntimespre = len(precipdata['time'])

#xray.Dataset.close(eventsdata)
#xray.Dataset.close(precipdata)

try:
    os.remove(DirO + FileO)
except OSError:
    pass

print "nevents here:", nevents
ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('events', nevents)
ncfile.createDimension('lats', nlats)
ncfile.createDimension('lons', nlons)

today = dt.date.today()
ncfile.history=('Created using process_full_output.py on ' +
            today.strftime('%d/%m/%y')  + 'with inputs: ' + Data +
            '; ' + Version + '; ' + str(startyr) + '-' + str(endyr))
varlistin = []

newvar = 'xmaxspeed_' + str(speedtspan) + 'ts'
varlistin.append((str(newvar)))

varlistin.extend(('gridboxspanSA','totalprecipSA','uniquegridboxspanSA','avg_intensity'))

varlistin.extend(('gridboxspan','totalprecip','uniquegridboxspan','timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax'))

f4vars = np.array(varlistin)
nf4vars = len(f4vars)

Ofilevars = []

# Add lats and lons

Ofilevars.append(ncfile.createVariable('lats','f4',('lats'),fill_value=-9999))
setattr(Ofilevars[-1],'units','degrees latitude')
setattr(Ofilevars[-1],'description',
                'latitude of indexes used in start, centre and end points')

Ofilevars.append(ncfile.createVariable('lons','f4',('lons'),fill_value=-9999))
setattr(Ofilevars[-1],'units','degrees longitude')
setattr(Ofilevars[-1],'description',
                'longitude of indexes used in start, centre and end points')


# Add all other variables
for ivar in range(0,nf4vars):
    Ofilevars.append(ncfile.createVariable(f4vars[ivar],'f4',('events'),fill_value=-9999))

    units,desc = getunitsdesc(f4vars[ivar])
    setattr(Ofilevars[-1],'units',units)
    setattr(Ofilevars[-1],'description',desc)

tminchunk = 0
tmaxchunk = 0

maxrun = minevent + nevents


ncfile['timespan'][:] = (maxt[:]-mint[:]+ 1) * tmult

ncfile['tstart'][:] = mint[:]
ncfile['xmin'][:] = minx[:]
ncfile['xmax'][:] = maxx[:]
ncfile['ymin'][:] = miny[:]
ncfile['ymax'][:] = maxy[:]

ncfile['xcenterstart'][:] = startx[:]
ncfile['xcenterend'][:] = endx[:]
ncfile['ycenterstart'][:] = starty[:]
ncfile['ycenterend'][:] = endy[:]

ncfile['tmean'][:] = meant[:]
ncfile['xcentermean'][:] = centerx[:]
ncfile['ycentermean'][:] = centery[:]
ncfile[newvar][:] = maxzonalspeed[:]

ncfile['lats'][:] = lats[:]
ncfile['lons'][:] = lons[:]

#for ievent in range(minevent, minevent + 3): ### For testing!!!

for ievent in range(minevent,maxrun):
    index = ievent-minevent
    tmin = int(max(0,mint[index]-1))
    tmax = int(min(ntimes,maxt[index]+1))

    ymin = int(max(0,miny[index]-latlonaddsize))
    ymax = int(min(nlats,maxy[index]+latlonaddsize))

    xmin = int(max(0,minx[index]-latlonaddsize))
    xmax = int(min(nlons,maxx[index]+latlonaddsize))

    if (xmin > xmax):
        print "wrapping around lons"
        xmin = 0
        xmax = nlons

    if (ymin > ymax):
        sys.exit("wrapping around lats - something is wrong!")

    if (tmax > tmaxchunk):
        print "tmax is " + str(tmax) + " whilst tmaxchunk is " + str(tmaxchunk)

        for attempt in range(10):
            try:
                tminchunk = int(max(tmin - 10,0))
                tmaxchunk = int(min(tmax + chunksize + 1,ntimes))
                print "tminchunk is now " + str(tminchunk) + "whilst tmin is " + str(tmin)
                print "tmaxchunk is now " + str(tmaxchunk) + "whilst tmax is " + str(tmax)
                eventschunk = eventsin.isel(time=slice(tminchunk,tmaxchunk)).values
                # Here we include a multiplier to get from mm/hr (or mm/day) to mm/timestep (usually 3 hrs)
                precipchunk = mult * precipin.isel(time=slice(tminchunk,tmaxchunk)).values
            except MemoryError:
                chunksize = chunksize/2
                print "memory error, reducing chunksize to ", chunksize
                continue
            break
        else:
            exit("couldn't find enough memory")

    if (ievent % 5000000 == 0):
        print "ievent: " + str(ievent)

    tminsel = int(max((tmin-tminchunk)-1,0))
    tmaxsel = int(min((tmax-tminchunk)+2,ntimes))

    eventsin_small = eventschunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]
    data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == ievent))
    SurfA_small = SurfA[ymin:ymax,xmin:xmax].values

    precipin_small = precipchunk[tminsel:tmaxsel,ymin:ymax,xmin:xmax]

    try:
        data_mask_max = np.amax(data_mask_small.mask,axis=0)
    except ValueError:
        print ievent
        print mint[ievent-minevent]
        print maxt[ievent-minevent]
        print tmin, tmax
        print tminchunk
        print tminsel, tmaxsel
        print ymin,ymax
        print xmin,xmax
        print ievent
        print data_mask_small
        print data_mask_small.shape
        exit("ValueError")

    ncfile['gridboxspan'][ievent-minevent] = (np.sum([data_mask_small.mask]))
    ncfile['gridboxspanSA'][ievent-minevent] = (np.sum([data_mask_small.mask*SurfA_small[None,:,:]]))

    # there used to be a factor of 0.001 here to convert from mm to m.
    # but we want to keep it in mm!
    #        ncfile['totalprecip'][ievent-minevent] = (0.001 * np.nansum([data_mask_small.mask * precipin_small]))
    ncfile['totalprecip'][ievent-minevent] = (np.nansum([data_mask_small.mask * precipin_small]))

    # SurfA is in m2
    ncfile['uniquegridboxspan'][ievent-minevent] = (np.sum(data_mask_max))
    ncfile['uniquegridboxspanSA'][ievent-minevent] = (np.sum([data_mask_max*SurfA_small]))
    # Factor of 0.001 to convert from mm to m
    # SurfS_small is in m2
    # So OtotalP is in m3
    ncfile['totalprecipSA'][ievent-minevent] = (0.001 * np.nansum([data_mask_small.mask * precipin_small * SurfA_small[None,:,:]]))

    ncfile['avg_intensity'][ievent-minevent] = (np.nansum([data_mask_small.mask * precipin_small]) / 
                                                  ( 3.0 * np.sum([data_mask_small.mask])))

ncfile.close()

timeend = time.time()
print "whole thing took " + str(timeend - timestart) + " seconds"


