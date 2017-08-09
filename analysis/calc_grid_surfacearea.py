# -*- coding: utf-8 -*-
"""
Script to calculate the surface area of gridded data.
The output from this script is used when summing up total precipitation and
total area of precipitation

Created: Oct 2016
Author: Rachel H White rhwhite@uw.edu
"""

import os, errno
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray
import Ngl
from scipy import stats
import math
import argparse
import resource

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')

args = parser.parse_args()

print args.Data

pi = 3.14159    # pi
rE = 6.371E6    # radius of earth in m

Data = args.Data[0]
startyr = args.startyr[0]
endyr = args.endyr[0]


if Data == "TRMM":
    filetimespan = "3hrly"
    DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
    FileP = 'TRMM_' + str(startyr) + '-' + str(endyr) + '_3B42_3hrly_nonan.nc'
elif Data == "TRMM_ERAIgd":
    filetimespan = "3hrly"
    DirP = '/home/disk/eos4/rachel/Obs/TRMM/' + filetimespan + '/'
    FileP = "regrid2ERAI_TRMM_3B42_" + str(startyr) + '-' + str(endyr) + ".nc"
elif Data == "ERAI":
    filetimespan = "3hrly"
    DirP = '/home/disk/eos4/rachel/Obs/ERAI/Precip_' + filetimespan + '/'
    FileP = 'ERAI_Totalprecip_' + str(startyr) + '-' + str(endyr) + '_preprocess.nc'
elif Data == "ERA20C":
    DirP = '/home/disk/eos4/rachel/Obs/ERA_20C/'
    FileP = 'ERA_20C_Totalprecip_' + str(startyr) + '-' + str(endyr) + '_preprocess.nc'
elif Data == "CESM":
    DirP = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
    FileP = 'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_' + str(startyr) + '-' + str(endyr) + '.nc'
elif Data == "GPCP":
    DirP = '/home/disk/eos4/rachel/Obs/GPCP/Daily/'
    FileP = 'GPCP_1DD_v1.2_199610-201510.nc'


print DirP + FileP

#Get lons and lats
FileIn = xray.open_dataset(DirP + FileP)

if Data in ["CESM",'GPCP']:
    lats = FileIn['lat']
    lons = FileIn['lon']
else:
    lats = FileIn['latitude']
    lons = FileIn['longitude']

#convert to radians
latsr = lats * pi / 180.0
lonsr = lons * pi / 180.0

nlats = len(lats)
nlons = len(lons)

area = np.zeros([nlats,nlons],np.float)
lonvalue = np.zeros(nlons,np.float)

# Almost all grids will have equal longitude spacing, but just in case:
for ilon in range(0,nlons-1):
    lonvalue[ilon] = abs(lonsr[ilon+1] - lonsr[ilon])

lonvalue[nlons-1] = abs(lonsr[nlons-1] - lonsr[nlons-2])

for ilat in range(0,nlats):
    print ilat
    # Based on: area above a latitude lat = 2piR^2(1 - sin(lat)
    # Thus area between two latitudes: 2piR^2(sin(lat1) - sin(lat2))
    # Break into 2pi and multiply by difference between lons: R^2(sin(lat1)-sin(lat2)) * (lon1 - lon2)
    if ilat == 0:
        latvalue = abs(np.sin(0.5*(latsr[ilat+1] + latsr[ilat]) - np.sin(latsr[ilat])))
    elif ilat == nlats-1:
        latvalue = abs(np.sin(latsr[ilat]) - np.sin(0.5 * (latsr[ilat-1] + latsr[ilat])))

    else:
        latvalue = abs(np.sin(0.5*(latsr[ilat] + latsr[ilat+1])) - np.sin(0.5 * (latsr[ilat]+ latsr[ilat-1])))

    for ilon in range(0,nlons): 
        area[ilat,ilon] = abs(rE * rE * latvalue * lonvalue[ilon])

        #old version:
        """
        if ilat == nlats-1:
               xlength = (rE) * (latsr[ilat] - latsr[ilat-1])
        else:
            xlength = (rE) * (latsr[ilat+1] - latsr[ilat])

        if ilon == nlons-1:
            # Assuming constant longitude spacing!
            ylength = (rE * math.cos(latsr[ilat])) * (lonsr[ilon] - lonsr[ilon-1])
        else:
                        ylength = (rE * math.cos(latsr[ilat])) * (lonsr[ilon-1] - lonsr[ilon])

        area[ilat,ilon] = xlength * ylength
        """

ncfile = Dataset(DirP + Data + "_SurfaceArea.nc", 'w')
ncfile.createDimension('lon', nlons)
ncfile.createDimension('lat', nlats)


SurfA = ncfile.createVariable('SurfaceArea','f4',('lat','lon'),fill_value=-9999)

SurfA[...] = area[...]

print np.sum(SurfA)

ncfile.close()

Ngl.end()




