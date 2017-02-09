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

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,
                    help = 'the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs='?',default=4,
                    help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound',type=float,nargs='+',
                    help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=float,nargs="+",
                    help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')

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

def isinregion(ilat,ilon,minlat,maxlat,minlon,maxlon,checklons):
    if ilat < minlat or ilat > maxlat:
        return(False)
    else:
        if checklons:
            if ilon < minlon or ilon > maxlon:
                return False
        else:
            return(True)

args = parser.parse_args()
print "here's what I have as arguments: ", args

if args.splittype[0] not in ["day","speed","maxspeed"]:
        exit("incorrect splittype " + str(args.splittype[0]) + 
            " must be day, speed, or maxspeed")
if args.speedtspan not in [0,1,4]:
        exit("incorrect speedtspan " + str(args.speedtspan[0]) + 
            " must be 0 or 1")
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
minlon = args.minlon
maxlon = args.maxlon
minlat = args.minlat
maxlat = args.maxlat
diradd = getdirectory(splittype)

nbounds = len(tbound1)

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

nevents = len(datain.events)

averageydist = np.zeros(nbounds)
averageyspeed = np.zeros(nbounds)
averageprecipperhr = np.zeros(nbounds)
averageprecipperareahr = np.zeros(nbounds)
averagetime = np.zeros(nbounds)

count = np.zeros(nbounds)

nevents = 100

timespan = datain.timespan[0:nevents].values
ycenterstart = datain.ycenterstart[0:nevents].values
ycenterend = datain.ycenterend[0:nevents].values
xcenterstart = datain.xcenterstart[0:nevents].values
xcenterend = datain.xcenterend[0:nevents].values

timespan = datain.timespan[0:nevents].values
totalprecip = datain.totalprecip[0:nevents].values
gridboxspan = datain.gridboxspan[0:nevents].values

startlats = lats[datain.ycenterstart[0:nevents].astype(int)]
endlats = lats[datain.ycenterend[0:nevents].astype(int)]

print startlats

for ievent in range(0,nevents):
    # check if in region
    ilat = startlats[ievent]
    if ilat < minlat or ilat > maxlat:
        continue
    else:
        ilat = endlats[ievent]
        if ilat < minlat or ilat > maxlat:
            continue
        else:

            for ibound in range(0,nbounds):
                if timespan[ievent] < tbound2[ibound]*24:
                    print (startlats[ievent], '-',
                            endlats[ievent])

                    # multiply by 24 to get hours
                    ydist = (ycenterend[ievent] - ycenterstart[ievent])
                    time = timespan[ievent]
                    precip = totalprecip[ievent]

                    averageydist[ibound] += ydist
                    averageyspeed[ibound] += ydist/time
                    # if negative then that's fine for NH, positive is fine
                    # for southern hemisphere, so get the average event
                    # distance travelled

                    averagetime[ibound] += time
                    averageprecipperhr[ibound] += (precip/time)
                    # Include factor of 3 to convert to hours, not timesteps
                    averageprecipperareahr[ibound] += (precip/(3.0 *
                                                        gridboxspan[ievent]))
                    count[ibound] += 1
                    break

print count
averageydist = averageydist / count
averageyspeed = averageyspeed / count
averagetime = averagetime / count
averageprecipperhr = averageprecipperhr / count
averageprecipperareahr = averageprecipperareahr / count


# Write out to a text file

with open(DirI + '/testDomainAverages_' + '{:d}'.format(minlat) + 'N-' +
            '{:d}'.format(maxlat) + 'N.txt', 'w') as text_file:
    text_file.write('Domain averages for ' + '{:d}'.format(minlat) + 'N-' +
            '{:d}'.format(maxlat) + 'N and ' + '{:d}'.format(minlon) + 'E-' +
            '{:d}'.format(maxlon) + 'E \n')
    text_file.write('Timespan (hours), averageydistance (degs), averageyspeed'+
            '(degs/hr), averagepreciphr (mm/hr), averagepreciphr (mm/gridbox/hr) \n')
    for ibound in range(0,nbounds):
        text_file.write('{0:d}'.format(tbound1[ibound]) + '-' +
                        '{0:d}'.format(tbound2[ibound]) + ' days,   ' +
                        '{0:f}'.format(averageydist[ibound]) + ' degrees,    ' +
                        '{0:f}'.format(averageyspeed[ibound]) + ' degrees/hr,    ' +
                        '{0:f}'.format(averagetime[ibound]) + ' hours,    ' +
                        '{0:f}'.format(averageprecipperhr[ibound]) + ' mm/hr,    ' +
                        '{0:f}'.format(averageprecipperareahr[ibound]) + 'mm/hr/gridbox \n')


datain.close()






