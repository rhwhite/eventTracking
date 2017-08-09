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
import time
import xray
from rhwhitepackages.query import query_yes_no
from rhwhitepackages.parallel import calcmissed
from rhwhitepackages.readwrite import xrayOpen
from multiprocessing import Process,Queue
from sys import getsizeof
import psutil
import argparse

#Thresholds are in mm/day, data is in mm/hr
# Switched code so now thresholds should be from smallest to largest
# Values lower than the 1st threshold are set to 0
#Threshold = [120, 80, 56, 40, 24, -9999]
# Threshold = [24, 40, 56, 80, 120]
#Threshold = [6,12,24,48,72,96]
#Threshold = [48,72,96,120,168]


parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')

args = parser.parse_args()

print "here's what I have as arguments: ", args


lowmem = True
Data = args.Data[0]
name = args.Version[0] #"Standard"
startyr = args.startyr[0]
endyr = args.endyr[0]

nyears = endyr-startyr + 1

#Defaults
latvar = 'latitude'
lonvar = 'longitude'
Dir = '/home/disk/eos4/rachel/Obs/' + Data + '/'
DirO = Dir + 'Thresholds/'

# Set dataset specific variables:
if Data == "TRMM":
    Filein = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'
    precipvar = 'pcp'
    precipunit = 'mm/hr'
    if name == "5th_from48_":
        Threshold = [48,72,96,120,168]
    elif name == "Standard":
        Threshold = [24, 40, 56, 80, 120]
    elif name == "6th_from6":
        Threshold = [6,12,24,48,72,96]
elif Data == "TRMMERAIgd":
    Dir = '/home/disk/eos4/rachel/Obs/TRMM/'
    DirO = '/home/disk/eos4/rachel/Obs/TRMM/Thresholds/ERAIgd/'
    Filein = "regrid2ERAI_TRMM_3B42_1998-2014.nc"
    precipvar = 'pcp'
    precipunit = 'mm/hr'

    if name == "Standard": 
        Threshold = [24, 40, 56, 80, 120]
    elif name == "2mmhr":
        Threshold = [48,72,96,120,148]
    else:
        exit("unknown name " + name)
elif Data == "ERAI":
    Dir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
    DirO = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/Thresholds/'
    Filein = 'ERAI_Totalprecip_1980-2015_preprocess.nc'
    precipvar = 'tpnew'
    precipunit = 'mm/hr'

    if name == "Standard":
        Threshold = [24, 40, 56, 80, 120]
    elif name == "low":
        Threshold = [12, 18, 26, 36, 48]
    elif name == "2mmhr":
        Threshold = [48,72,96,120,148]
    elif name == "vlow":
        Threshold = [2, 4, 8, 16, 32, 48]
    else:
        exit("unknown name " + name)
elif Data == "ERA20C":
    Filein = 'ERA_20C_Totalprecip_' + str(startyr) + '-' + str(endyr) + '_preprocess.nc'
    precipvar = 'tpnew'
    precipunit = 'mm/hr'
    if name == "20Cstd":
        Threshold = [4,8,16,32,64,128]
    elif name == "std":
        Threshold = [12,24,48,72,96]

    else:
        exit("unknown name" + name)
elif Data == "GPCP":
    Dir = '/home/disk/eos4/rachel/Obs/GPCP/Daily/'
    DirO = '/home/disk/eos4/rachel/Obs/GPCP/Daily/Thresholds/'
    Filein = 'GPCP_1DD_v1.2_199610-201510.nc'
    precipvar = 'PREC'
    precipunit = 'mm/day'
    lonvar = 'lon'
    latvar = 'lat'

    if name == "Daily12":
        Threshold = [12, 20, 28, 40, 60]
    else:
        exit("unknown name " + name)
else:
        exit("unknown data type " + Data)

# Make output directory if it doesn't already exist
try:
    os.mkdir(DirO + '/' + name + '_' + str(startyr) + '-' + str(endyr))
except OSError:
    pass

DirO = DirO + '/' + name + '_' + str(startyr) + '-' + str(endyr) + '/'

# Convert thresholds if necessary

nThresh = len(Threshold)
if precipunit == 'mm/day':
    Thresh = Threshold
elif precipunit == 'mm/hr':
    Thresh = [x/24.0 for x in Threshold]    # convert to mm/hr
else:
    exit("unexpected precipunit " + precipunit)


# Read in precip data
print Dir + Filein
precipin = xrayOpen(Dir + Filein)

latsin = precipin[latvar][:]
lonsin = precipin[lonvar][:]
precipnew=precipin[precipvar]

nlons = len(lonsin)
nlats = len(latsin)
ntimes = len(precipnew.time)
print ntimes
chunksize = 10000

print "Hold on, I'm calculating the total amount of precip missed by your chosen thresholds, estimating the value using the first year"
check = query_yes_no("Do you want to confirm whether or not it's acceptable?")

lowestTH = Thresh[0]
totals = []
missed = []

startTime = time.time()

try:
    precipchunk = precipnew.values
    total = np.nansum(precipchunk)
    missedvalues = np.ma.sum(np.ma.masked_greater_equal(precipchunk[np.isfinite(precipchunk)],Thresh[0]))
    totals.append(total)
    missed.append(missedvalues)

except MemoryError:
    print 'memory error, doing in chunks'
    for attempt in range(10):
        try:
            nits = np.ceil(ntimes/chunksize)
            print str(nits) + "iterations required"
            for i in range(0,np.int(nits)):
                print i
                precipchunk = np.float32(precipnew.isel(time=slice(i*chunksize,np.amin([(i+1)*chunksize,ntimes]))))
                total = np.nansum(precipchunk)
                missedvalues = np.ma.sum(np.ma.masked_greater_equal(precipchunk[np.isfinite(precipchunk)],Thresh[0]))
                totals.append(total)
                missed.append(missedvalues)
        except MemoryError:
            chunksize = chunksize/2
            print chunksize
            continue
        break



print totals
totalmissed = 100.0* sum(missed)/sum(totals)

endTime = time.time()
workTime = endTime - startTime
print "job took " + str(workTime) + " seconds"

print missed
print "this threshold misses " + str(sum(missed)) + " which is " + str(totalmissed) + " % of the total rain"

if check:
    accept = query_yes_no("is this acceptable?")
    if not accept: exit("you didn't accept the thresholds")

timesin = precipin['time'][:]

for year in range(startyr,endyr+1):    #(1998,2015)
    print year
    stryear = str(year)
    counter = 0 #reset counter for each year

    precipchunk = precipnew.sel(time=str(year))
    timeschunk = timesin.sel(time=str(year))
    ntimes = len(timeschunk)

    for x in np.nditer(precipchunk, op_flags=['readwrite']):
        if x < Thresh[0]:   # Deal with values below 1st threshold (0)
            x[...] = 0
        elif x >= Thresh[nThresh-1]:    # Deal with values above last threshold (nThresh - 1)
            x[...] = nThresh
        else:
            for ith in range(0,nThresh-1):  # iterate through thresholds 0 to nThresh - 2
                if x >= Thresh[ith] and x < Thresh[ith+1]:
                    x[...] = ith + 1
                    break

#       for y in np.nditer(precip_thrsh):
#           if y > 5:
#           print 'missing value'
#           print y
#           elif y < 0:
#           print 'weird value'
#           print y

    maxcount = ntimes
    print maxcount
    for count in range(0,maxcount):
        filecount = str(counter)
        if Data == "TRMM":
            FileO = 'TRMM_3B42_3hrly_' + name + stryear + '_' + filecount.zfill(5) + '.nc'
        elif Data == "TRMMERAIgd":
            FileO = 'TRMMERAIgd_3B42_3hrly_' + name + stryear + '_' + filecount.zfill(5) + '.nc'
        elif Data == "ERAI":
            FileO = 'ERAI_3hrly_' + name + "_" + stryear + '_' + filecount.zfill(5) + '.nc'
        elif Data == "ERA20C":
            FileO = 'ERA20C_3hrly_' + name + "_" + stryear + '_' + filecount.zfill(5) + '.nc'
        elif Data == "GPCP":
            FileO = 'GPCP_' + name + "_" + stryear + '_' + filecount.zfill(5) + '.nc'

        try:
            os.remove(DirO + FileO)
        except OSError:
           pass

        ncfile = Dataset(DirO + FileO, 'w')
        ncfile.createDimension('z', 1)
        ncfile.createDimension('y',len(latsin))
        ncfile.createDimension('x',len(lonsin))
        ncfile.createDimension('time',1)

        precipTH = ncfile.createVariable('value','f4',('z','y','x'),fill_value=0)

        precipTH[:,:,:] = precipchunk[count,:,:].values
        precipTH.grid_type = 'linear'
        precipTH.grid_name = 'grid-1'


        #global attribute
        ncfile.Conventions = 'COARDS'
        ncfile.calendar = 'standard'
        ncfile.center = 'gsfc'
        ncfile.Thresholds = 'Thresholds: ' + str(Threshold)
        ncfile.Missed = 'missed by thresholds: ' + str(totalmissed) + '% of total rain'
        ncfile.close()
        counter = counter + 1
