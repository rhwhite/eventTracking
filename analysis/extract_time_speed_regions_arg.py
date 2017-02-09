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
import re
import sys
import Ngl
import xray
import math
import resource
import argparse
from rhwhitepackages.readwrite import getunitsdesc
from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.readwrite import getdirectory

rsrcV = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit starts as  :', soft
print 'Hard limit starts as :', hard
resource.setrlimit(rsrcV, (120000000000, hard)) #limit memory usage
#                          137438953472
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit changed to :', soft

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
parser.add_argument('--SA',type=int,nargs='?',default=1,help='do surface area variables (1) or not (0)?')
parser.add_argument('--minGB',type=int,nargs='?',default=-1,help='minimum number of unique grid boxes to count as an event')

parser.add_argument('--minlat',metavar='minlat',type=int,nargs=1,
                    help='minimum latitude')
parser.add_argument('--maxlat',metavar='maxlat',type=int,nargs=1,
                    help='maximum latitude')
 

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
tbound1 = args.tbound1
tbound2 = args.tbound2
unit = args.unit[0]
Data = args.Data[0]
Version = args.Version[0]
filetimespan = args.filetspan[0]
startyr = args.startyr[0]
endyr = args.endyr[0]
SA = args.SA
minGB = args.minGB

diradd = getdirectory(splittype)

nbounds = len(tbound1)

R = 6371000     # radius of Earth in m

nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

minevent = 100000

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
FileI1 = 'All_Precip_' + str(startyr) + '-' + str(endyr) + '_' + Data + '_' + Version + '.nc'


#Get lons and lats
iday = 0
print FileInLats
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

nlats = len(lats)
nlons = len(lons)


print DirI + FileI1
datain = xrayOpen(DirI + FileI1,decodetimes=False)

varlistin = []

newvar = "xmaxspeed_" + str(speedtspan) + "ts"
varlistin.append((str(newvar)))

if SA == 1:
        varlistin.extend(('gridboxspanSA','totalprecipSA','uniquegridboxspanSA'))

varlistin.extend(('timespan','tstart','tmean','xcenterstart','xcenterend','ycenterstart','ycenterend','xcentermean','ycentermean','xmin','xmax','ymin','ymax','ydist'))

print varlistin

#if speedtspan == 0:    # have updated process script to not use this variable name now
#       newvar = "xmaxtravel"
#       newvarout = "xmaxspeed"

f4vars = np.array(varlistin)
f4varsin = np.array(varlistin)

f8vars = np.array(['eventid'])

nf4vars = len(f4vars)
nf8vars = len(f8vars)

nevents = len(datain.events)
chunksize = nevents

def runchunk(nevents,chunksize,writeidxchunk):
    for chunkevent in range(0,nevents,chunksize):
        chunkmax = np.amin([chunkevent+chunksize,nevents])
        print("chunkmax is: ", chunkmax)
        variablesin = np.zeros([nf4vars,chunksize],np.float)

        splitvar,SAgridboxes = getsplitvar(chunkevent,chunkmax)
        for ivar in range(0,nf4vars):
            variablesin[ivar,:] = datain[f4varsin[ivar]][chunkevent:chunkmax].values

        print "chunksize: ", chunksize, ", chunkevent: ",chunkevent

        varforwrite = np.zeros([nf4vars+nf8vars,chunksize*neventmult],np.float)         # saving memory by using neventmult

        writeidx = 0
 
        for ievent in range(chunkevent,chunkmax):
            if latstart[ievent-chunkevent] < minlat or
                latstart[ievent-chunkevent] > maxlat:
                continue
            if latend[ievent-chunkevent] < minlat or
                latstart[ievent-chunkevent] > maxlat:
                continue
            if SAgridboxes[ievent-chunkevent] > minGB:      # if unique gridboxes exceed specified threshold
                if tbound1[ibound] < 0:
                    if (splitvar[ievent-chunkevent] >= tbound1[ibound]) and (splitvar[ievent-chunkevent] < tbound2[ibound]):
                        for ivar in range(0,nf4vars):
                            varforwrite[ivar,writeidx] = variablesin[ivar][ievent-chunkevent]
                        varforwrite[nf4vars,writeidx] = ievent
                        writeidx += 1
                else:
                    if (splitvar[ievent-chunkevent] > tbound1[ibound]) and (splitvar[ievent-chunkevent] <= tbound2[ibound]):
                        for ivar in range(0,nf4vars):
                            varforwrite[ivar,writeidx] = variablesin[ivar][ievent-chunkevent]
                        varforwrite[nf4vars,writeidx] = ievent
                        writeidx += 1

        if writeidx > 0:
            if writeidx == 1:
                for ifvar in range(0,nf4vars+1):
                    Ofilevars[ifvar][writeidxchunk] = varforwrite[ifvar,0]
            else:
                for ifvar in range(0,nf4vars+1):
                    Ofilevars[ifvar][writeidxchunk:writeidxchunk+writeidx] = varforwrite[ifvar,0:writeidx]
        writeidxchunk += writeidx

def getsplitvar(chunkfirst,chunkmax):
    if splittype == "day":
        try:
            if datain['timespan'].units == "hours":
                # check that it really is hours, and not 3hrly units!
                if np.any(datain['timespan'][0:1000].values < 3):
                    exit("Timespan units says hours, but it doesn't seem to be! Unless your data has hourly resolution, in which case you need to update this script")
                splitvar = datain['timespan'][chunkfirst:chunkmax].values/24.0       # Divide by 24 to get values in days instead of hours 
            else:
                exit("not set up for anything other than hourly right now")
        except AttributeError:
            print("timespan has no units. You should really rerun process_output and write out units for your variables. For now, I will assume hours, but on your head be it if that's wrong")
            splitvar = datain['timespan'][chunkfirst:chunkmax].values/24.0

        latstarts = lats[datain['ycenterstart'][chunkfirst:chunkmax]]
        latends = lats[datain['ycenterend'][chunkfirst:chunkmax]]


    elif splittype == "speed":
        #xcenterstart = datain['xcenterstart'][chunkfirst:chunkmax].values
        #xcenterend = datain['xcenterend'][chunkfirst:chunkmax].values
        #ycenterstart = datain['ycenterstart'][chunkfirst:chunkmax].values
        #ycenterend = datain['ycenterend'][chunkfirst:chunkmax].values
        # Calculate zonal speeds
        lats1 = np.radians(lats[datain['ycenterstart'][chunkfirst:chunkmax]])
        lats2 = np.radians(lats[datain['ycenterend'][chunkfirst:chunkmax]])
        lons1 = np.radians(lons[datain['xcenterstart'][chunkfirst:chunkmax]])
        lons2 = np.radians(lons[datain['xcenterend'][chunkfirst:chunkmax]])

        latsmean = (lats1 + lats2) / 2.0

        az = (np.cos(latsmean) * np.cos(latsmean) * np.power(np.sin((lons2-lons1)/2),2))
        cz = 2.0 * np.arctan(np.sqrt(az),np.sqrt(1-az))

        distancez = R * cz

        pi = 3.14
        # Calculate speed
        angle = np.arctan2((lats2-lats1),(lons2-lons1))
        ones = np.zeros(lons1.shape)
        ones[...] = 1.0
        negones = ones * -1.0

        direction = np.where(lons2>=lons1,ones,negones) # True where 
        splitvar = direction * distancez/(timespan*3.0*60.0*60.0)

    elif splittype == "maxspeed":
        #if speedtspan == 0:
        #   splitvar = datain['xmaxtravel'].values
        #else:
        splitvar = datain['xmaxspeed_' + str(speedtspan) +  'ts'][chunkfirst:chunkmax].values
    else:
        exit("unacceptable splittype: " + str(splittype))

    # for all, get surface area in grid boxes:
    SAgridboxes = datain['uniquegridboxspan'][chunkfirst:chunkmax].values

    return splitvar,SAgridboxes


for ibound in range(0,nbounds):
    if splittype == 'day' and tbound1[ibound] > 0:
        neventmult = 0.25        # Assume fewer than neventmult of events will be caught by whatever bound we have
    elif splittype == 'MaxSpeeds' and (tbound1[ibound]/tbound2[ibound] > 0) and tbound1[ibound] != 0 and tbound2[ibound] != 0:
        neventmult = 0.25
    else:
        neventmult = 1.0

    varforwrite = np.zeros([nf4vars+nf8vars,nevents*neventmult],np.float)
    #print memory_usage_psutil()


    if splittype == "maxspeed":
        fileadd = "MaxSpeeds_" + str(speedtspan) + "ts_"
    elif splittype == "speed":
        fileadd = "Speeds_"
    elif splittype == "day":
        fileadd = "Sizes_"

    if minGB > 0:
        fileaddGB = '_min' + str(minGB) + 'GB'
    else:
        fileaddGB = ''

    if tbound1[ibound] < 0:
        tboundtitle = str(tbound1[ibound]) + '-' + str(tbound2[ibound])
    else:
        tboundtitle = str(tbound1[ibound]) + '-' + str(tbound2[ibound])


    FileO = ('testPrecip_' + fileadd + tboundtitle + unit + "_" + str(startyr) +
            '-' + str(endyr) + '_' + Version + fileaddGB + '_' + str(minlat) +
            'N-' + str(maxlat) + 'N.nc')

    print FileO
    # Create new file 
    ncfile = Dataset(DirO + FileO, 'w')
    ncfile.createDimension('events', None)

    Ofilevars = []
    for ivar in range(0,nf4vars):
        Ofilevars.append(ncfile.createVariable(f4vars[ivar],'f4',('events'),fill_value=-9999))

        units,desc = getunitsdesc(f4vars[ivar])
        setattr(Ofilevars[-1],'units',units)
        setattr(Ofilevars[-1],'description',desc)

    # create last variable 
    Ofilevars.append(ncfile.createVariable('eventid','f8',('events'),fill_value=-9999))

    #Extract events of startlen-endlen days
    # Do in chunks

    writeidxchunk = 0
    for attempt in range(10):
        try:
            print chunksize
            runchunk(nevents,chunksize,writeidxchunk)
            print("ran chunk")
        except MemoryError:
            #print 'memory usage', memory_usage_psutil()
            try:
                del(variablesin,splitvar,varforwrite)
            except NameError:
                pass

            chunksize = int(chunksize/4)
            print "failed on memory, quartering chunksize"
            print chunksize
            continue

        print 'success'
        break

    ncfile.close()
datain.close()






