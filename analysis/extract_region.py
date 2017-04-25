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
#import Ngl
import xray
import math
import resource
import argparse
from rhwhitepackages.readwrite import getunitsdesc
from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.readwrite import getdirectory

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--Name',type=str,nargs=1,help='name of region')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI,ERA20C, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--filetspan',type=str,nargs='?',default=['3hrly'],help='string for file time resolution, 3hrly etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--minlon',type=int,nargs='?',default=-180,
                    help='minimum longitude for extraction')
parser.add_argument('--maxlon',type=int,nargs='?',default=180,
                    help='maximum longitude for extraction')
parser.add_argument('--minlat',type=int,nargs='?',default=-90,
                    help='minimum latiitude for extraction')
parser.add_argument('--maxlat',type=int,nargs='?',default=90,
                    help='maximum latitude for extraction')
## Define functions

def createfile(filenameM, f4varsin, f8varsin):

    filenameM.createDimension('events', None)

    filevarsM = []
    for ivar in range(0,len(f4varsin)):
            filevarsM.append(filenameM.createVariable(f4varsin[ivar],'f4',('events'),fill_value=-9999))
            units,desc = getunitsdesc(f4vars[ivar])
            setattr(filevarsM[-1],'units',units)
            setattr(filevarsM[-1],'description',desc)

    for ivar in range(0,len(f8varsin)):
            filevarsM.append(filenameM.createVariable(f8varsin[ivar],'f8',('events'),fill_value=-9999))
            setattr(filevarsM[-1],'description','event id number')

    return (filevarsM)



def isinregion(ilat,ilon,minlat,maxlat,minlon,maxlon,checklons):
    if ilat < minlat or ilat > maxlat:
        return(False)
    else:
        if checklons:
            if ilon < minlon or ilon > maxlon:
                return False
            else:
                return(True)
        else:
            return(True)

args = parser.parse_args()
print "here's what I have as arguments: ", args

if args.Data[0] not in ['TRMM','TRMMERAIgd','ERAI','ERA20C','CESM']:
        exit("incorrect Data option " + str(args.Data[0]) + " must be TRMM, TRMMERAIgd, ERAI,ERA20C or CESM")

filename = args.Name[0]
Data = args.Data[0]
Version = args.Version[0]
filetimespan = args.filetspan[0]
startyr = args.startyr[0]
endyr = args.endyr[0]
minlon = args.minlon
maxlon = args.maxlon
minlat = args.minlat
maxlat = args.maxlat

R = 6371000     # radius of Earth in m

nyears = endyr - startyr + 1

checklons = True
if minlon == 0 and maxlon == 360:
    checklons = False
elif minlon ==-180 and maxlon == 180:
    checklons = False
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

FileI1 = 'All_Precip_' + str(startyr) + '-' + str(endyr) + '_' + Data + '_' + Version + '.nc'

#Get lons and lats
iday = 0
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

nevents = len(datain.events)

f4vars = np.array(['gridboxspanSA','totalprecipSA','uniquegridboxspanSA',
                   'gridboxspan','totalprecip','uniquegridboxspan','timespan',
                   'tstart','tmean','xcenterstart','xcenterend','ycenterstart',
                   'ycenterend','loncentermean','latcentermean',
                   'xmin','xmax','ymin','ymax'])

f8vars = np.array(['eventid'])

nf4vars = len(f4vars) 

startlats = lats[datain.ycenterstart[0:nevents].astype(int)]
endlats = lats[datain.ycenterend[0:nevents].astype(int)]
startlons = lons[datain.xcenterstart[0:nevents].astype(int)]
endlons = lons[datain.xcenterend[0:nevents].astype(int)]

midlats = lats[datain.ycentermean[0:nevents].astype(int)]
midlons = lons[datain.xcentermean[0:nevents].astype(int)]


invars = np.zeros([len(f4vars),nevents],np.float)
#read in variables
for ivar in range(0,len(f4vars)):
    if f4vars[ivar] not in ['loncentermean','latcentermean']:
        invars[ivar,:] = datain[f4vars[ivar]].values


## Create new files and put these file handles in arrays

filevars = []

ncfiles = []

FileO = filename + '_Precip_Events_' + Data  + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

ncfileO = Dataset(DirI + FileO, 'w')
filevars = createfile(ncfileO,f4vars,f8vars)

print len(filevars)
writeidx = 0
for ievent in range(0,nevents):
    # check if in region

    # center HAS to be in region
    if (isinregion(midlats[ievent],midlons[ievent],minlat,maxlat,minlon,maxlon,checklons)):

        # Plus either start OR end
        if (isinregion(startlats[ievent],startlons[ievent],minlat,maxlat,minlon,maxlon,checklons)
            or
            isinregion(endlats[ievent],endlons[ievent],minlat,maxlat,minlon,maxlon,checklons)):


            # print to file
            for ivar in range(0,len(f4vars)):
                if f4vars[ivar] not in ['loncentermean','latcentermean']:
                    filevars[ivar][writeidx] = invars[ivar,ievent]
                else:
                    if f4vars[ivar] == 'loncentermean':
                        filevars[ivar][writeidx] = midlons[ievent]
                    elif f4vars[ivar] == 'latcentermean':
                        filevars[ivar][writeidx] = midlats[ievent]

            # print last variable, ievent
            filevars[nf4vars][writeidx] = ievent + minevent
            writeidx += 1

    if ievent % 10000 == 0:
        print ievent

datain.close()






