# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os, errno
import numpy as np
import netCDF4
from netCDF4 import Dataset
import datetime as dt
import gc
import re
import sys
import Ngl
import xray as xr
import math
import resource
import argparse
import calendar
from rhwhitepackages.readwrite import getunitsdesc
from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.readwrite import getdirectory

#rsrcV = resource.RLIMIT_AS
#soft, hard = resource.getrlimit(rsrcV)
#print 'Soft limit starts as  :', soft
#print 'Hard limit starts as :', hard
#resource.setrlimit(rsrcV, (137000000000, hard)) #limit memory usage
#                          137438953472
#soft, hard = resource.getrlimit(rsrcV)
#print 'Soft limit changed to :', soft

def memory_usage_psutil():
    # return the memory usage in MB
    import psutil
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem

#print memory_usage_psutil()

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs='?',default=4,help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound',type=float,nargs='+',help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=float,nargs="+",help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI,ERA20C, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--filetspan',type=str,nargs='?',default=['3hrly'],help='string for file time resolution, 3hrly etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--minGB',type=int,nargs='?',default=-1,help='minimum number of unique grid boxes to count as an event')


args = parser.parse_args()
print "here's what I have as arguments: ", args

if args.splittype[0] not in ["day","speed","maxspeed"]:
        exit("incorrect splittype " + str(args.splittype[0]) + " must be day, speed, or maxspeed")
if args.speedtspan not in [0,1,4]:
        exit("incorrect speedtspan " + str(args.speedtspan[0]) + " must be 0 or 1")
if args.Data[0] not in ['TRMM','TRMMERAIgd','ERAI','ERA20C','CESM']:
        exit("incorrect Data option " + str(args.Data[0]) + " must be TRMM, TRMMERAIgd, ERAI,ERA20C or CESM")

splittype = args.splittype[0]
speedtspan = args.speedtspan
unit = args.unit[0]
Data = args.Data[0]
Version = args.Version[0]
filetimespan = args.filetspan[0]
startyr = args.startyr[0]
endyr = args.endyr[0]
minGB = args.minGB

if filetimespan == '3hrly':
    daymult = 8     # 8 timesteps per day
elif filetimespan == 'hrly':
    daymult = 24    # 24 timesteps per day
else:
    print(filetimespan)
    exit('unknown filetimespan')


if splittype == 'day' and unit == 'day':
    tbound1 = np.array(args.tbound1) * 24.0    # multiple by 24 to get hours, not days
    tbound2 = np.array(args.tbound2) * 24.0


print tbound1
diradd = getdirectory(splittype)

nbounds = len(tbound1)

nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

nmonths = 12    # months per year
dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
minevent = 100000

if Data == 'TRMM':
    DirRaw = ('/home/disk/eos4/rachel/EventTracking/FiT_RW/' + Data + '_output/' +
            Version + str(startyr) + '/raw/')
else:
    DirRaw = ('/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data +
                '_output/' + Version + str(startyr) + '/raw/')
DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' + Version + str(startyr) + '/proc/'


if Data == "TRMM":
    if Version == '6th_from6' or Version == '5th_from48':
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/proc/'
    FileInLats = '/home/disk/eos4/rachel/Obs/TRMM/SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc'
elif Data == "TRMMERAIgd":
    FileInLats = '/home/disk/eos4/rachel/Obs/TRMM/regrid2ERAI_TRMM_3B42_1998-2014.nc'

elif Data == "ERAI":
    FileInLats = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/SeasAnn_ERAI_Totalprecip_' + str(startyr) + '-' + str(endyr) + '_preprocess.nc' 

elif Data == "ERA20C":
    FileInLats = '/home/disk/eos4/rachel/Obs/ERA_20C/ERA_20C_LatLon.nc'

elif Data == "CESM":
    DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '/proc/' 
    FileInLats = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'
else:
    print("unexpected data type")
    exit()

DirO = DirI + diradd + '/'
# open files: 
# dataIn is list of all precip events
# eventsIn is maps of each timestep with event id number
#eventsIn = xrayOpen(DirRaw + '/ts_' + Data +  str(startyr) + '-' + str(endyr) + '_' + Version + '_4Dobjects.nc')
eventsIn = xrayOpen(DirRaw + '/ts_' + Data + '_' + Version + '_' + 
                    str(startyr) + '-' + str(endyr)
                    + '_4Dobjects.nc')

dataIn = xrayOpen(DirI + '/All_Precip_' + str(startyr) + '-' + str(endyr) + '_' + Data + '_' + Version + '.nc',decodetimes=False)

#Get lons and lats
#print FileInLats
FileIn = xrayOpen(FileInLats)

if Data == "CESM":
    lats = FileIn['lat'].values
    lons = FileIn['lon'].values
elif Data in ["ERA20C","TRMMERAIgd"]:
        lats = FileIn['latitude'].values
        lons = FileIn['longitude'].values
else:
    lats = FileIn['Latitude'].values
    lons = FileIn['Longitude'].values

print lats

nlats = len(lats)
nlons = len(lons)


# get data for all events
nevents = len(dataIn.events)

# test here

evDurations = dataIn['timespan'].values
evGridboxes = dataIn['gridboxspan'].values
evTstarts = dataIn['tstart'].values
evYmins = dataIn['ymin'].values
evYmaxs = dataIn['ymax'].values
evXmins = dataIn['xmin'].values
evXmaxs = dataIn['xmax'].values

maxtimes = len(evDurations)
# start on first year, first month
iyear = 0
imonth = 0

# assuming first month is January for nextmonth initialization
print('this code assumes the first month is January')
ntimes = nyears * 12

dtime = 0   # initialize days since as 0
timestamp = np.zeros(ntimes,'float')
for iyear in range(0,nyears):
    for imonth in range(0,nmonths):
        itime = (iyear*12) + imonth
        timestamp[itime] = dtime
        dtime = dtime + dayspermonth[imonth]
        if (calendar.isleap(iyear+startyr) and imonth == 1):
            dtime = dtime + 1

for iday in range(0,nbounds):

    nextmonth = 1 + 31 * daymult     # number of timesteps until next month
    # Read in events one month at a time
    # Add extra 100 to catch tails of events in different months
    events = eventsIn['value'][0:min(nextmonth+100,maxtimes),0,:,:].values


    locDensity = None
    locFracDensity = None
    # Set up arrays
    locDensity = np.zeros([ntimes,nlats,nlons])
    locFracDensity = np.zeros([ntimes,nlats,nlons])

    imonth = 0
    calmonth = 0
    Tmin = 0
    for ievent in range(0,nevents):
        # if Tstart is in next month, increment imonth and nextmonth
        if evTstarts[ievent] >= nextmonth:
            imonth += 1
            calmonth +=1
            if calmonth >= 12:  # if reach end of year, reset calmonth
                calmonth = 0
            # increment nextmonth
            if (calendar.isleap(iyear + startyr) and calmonth == 1):
                # calmonth = 1 is february
                # if leap year, and Feb, add one day
                nextmonth = nextmonth + ((dayspermonth[calmonth] + 1) * daymult)
            else:
                nextmonth = nextmonth + (dayspermonth[calmonth] * daymult)

            # update events map
            Tmin = evTstarts[ievent]    # to get correct index later
            events = None   # clear memory
            print 'memory usage', memory_usage_psutil()
            events = (eventsIn['value']
                        [evTstarts[ievent]:min(nextmonth+50,maxtimes),0,:,:].values)
        if imonth >= ntimes:
            print 'nextmonth ',nextmonth
            print 'imonth ', imonth
            break
        # if event is big enough to be counted
        if evGridboxes[ievent] > minGB: 
            # get duration bound index
            evindex = ievent + minevent
            if splittype == 'day':
                evduration = evDurations[ievent]
                # Is event correct duration for this iday
                if evduration >= tbound1[iday] and evduration < tbound2[iday]:
                    # Old code
                    #    durindex = 0    # start with smallest events
                    #while evduration >= tbound2[durindex]:
                    #    durindex += 1       # increment by one until find index

                    tminsel = max(0,evTstarts[ievent]-1 - Tmin)
                    tmaxsel = min((tminsel + evduration)+2,maxtimes)

                    yminsel = max(0,evYmins[ievent]-1)
                    ymaxsel = min(nlats,evYmaxs[ievent]+2)

                    xminsel = max(0,evXmins[ievent]-1)
                    xmaxsel = min(nlons,evXmaxs[ievent]+2)

                    eventsin_small = (events[tminsel:tmaxsel,
                                            yminsel:ymaxsel,
                                            xminsel:xmaxsel])

                    data_mask_small = np.ma.array(eventsin_small,mask=(eventsin_small == evindex))
                    # take max to not count multiple timesteps in one location as more than one
                    # event
                    try:
                        data_mask_max = np.amax(data_mask_small.mask,axis=0)
                        data_mask_sum = np.sum(data_mask_small.mask,axis=0)
                    except ValueError:
                        print tminsel, tmaxsel
                        print evTstarts[ievent]
                        print evduration
                        print eventsin_small.shape
                        print eventsin_small
                        print data_mask_small
                        print evindex
                        exit()
                    # Add to array
                    try:
                        locDensity[imonth,yminsel:ymaxsel,xminsel:xmaxsel] += data_mask_max
                        locFracDensity[imonth,yminsel:ymaxsel,xminsel:xmaxsel] +=(
                                                   data_mask_sum/evGridboxes[ievent])
                    except ValueError:
                        print(locDensity[imonth,yminsel:ymaxsel,xminsel:xmaxsel].shape)
                        print tminsel, tmaxsel
                        print evTstarts[ievent]
                        print evduration
                        print eventsin_small.shape
                        print eventsin_small
                        print data_mask_small
                        print evindex
                        exit()

            else:
                exit('not set up for anything except duration splits' + 
                        ' you will need to change the code')

    datimestamp = xr.DataArray(timestamp,coords=[('time',range(0,ntimes))],
                    attrs=[('units','days since ' + str(startyr) + '-01-01')])
    daout = xr.DataArray(locDensity,
                            coords=[
                                    ('time',datimestamp),
                                    ('lat',lats),
                                    ('lon',lons)],
                            attrs=[('description','local density of events'),
                                    ('units','event passes/gridbox')])

    daout2 = xr.DataArray(locFracDensity,
                            coords=[
                                    ('time',datimestamp),
                                    ('lat',lats),
                                    ('lon',lons)],
                            attrs=[('description','local fractional event density'
                                        + ' where event density is spread among'
                                        + ' all gridboxes the event passed ' 
                                        + ' weighted by time spent at each gridbox'),
                                    ('units','event presence/gridbox')])

    newDataset = xr.Dataset({'LocalDensity':daout}
                            ,{'LocalFractionalDensity':daout2})
    FileO = ('LocFracDensity_' + str(startyr) + '-' + str(endyr) + '_' + Version +
                   '_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '.nc')

    newDataset.to_netcdf(DirO + FileO,mode='w')

eventsIn.close()
dataIn.close()








