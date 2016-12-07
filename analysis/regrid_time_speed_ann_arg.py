#i-*- coding: utf-8 -*-
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
import Ngl
import xray
import math
import argparse
import resource

from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.readwrite import getdenfilename

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs=1,help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound1',type=float,nargs='+',help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=float,nargs="+",help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--sumlats',type=int,nargs='+',help='number of lats to sum over')
parser.add_argument('--sumlons',type=int,nargs='+',help='number of lons to sum over')
parser.add_argument('--minGB',type=int,nargs='?',default=0,help='minimum number of gridboxes to count as an event')


args = parser.parse_args()

print args.Data

splittype = args.splittype[0]
speedtspan = args.speedtspan[0]
tbound1 = args.tbound1
tbound2 = args.tbound2
unit = args.unit[0]
Data = args.Data[0]
Version = args.Version[0]
startyr = args.startyr[0]
endyr = args.endyr[0]
sumlats = args.sumlats
sumlons = args.sumlons
minGB = args.minGB

nbounds = len(tbound1)
nsumlats = len(sumlats)
nyears = endyr-startyr + 1 

mapping = "center"

if splittype == "maxspeed":
    diradd = "MaxSpeeds"
elif splittype == "speed":
    diradd = "Speeds"
elif splittype == "day":
    diradd = "Sizes"
else:
    exit("unexpected splittype")

DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' + Version + str(startyr) + '/proc/' + diradd + '/'
DirO = DirI

if Data == "TRMM":
    if Version == '5th_nanto25':
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
    elif Version == '5th_nantozero':
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
    elif Version == '7thresh':
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
    elif Version == '5th_from48':
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'
        DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'
#elif Data == 'ERAI':   no longer need this correction - data are correct
#   nyears = nyears-1


sumvars = ["TDensity","EDensity","WDensity","SDensity"]

if minGB > 0:
    GBadd = '_min' + str(minGB) + 'GB'
else:
    GBadd = ''

for isumlat in range(0,nsumlats):
    sumlat = sumlats[isumlat]
    sumlon = sumlons[isumlat]
    for ibound in range(0,nbounds):

        FileI = getdenfilename(mapping, Data, Version, startyr, endyr, ibound, splittype, unit, speedtspan, minGB, tbound1, tbound2, 'Ann',-1,-1)

        FileO = FileI.rstrip('.nc') + '_regrid_' + str(sumlon) + 'lons_' + str(sumlat) + 'lats.nc'

        print FileO

        filetimespan = "3hrly"

        file = xrayOpen(DirI + FileI)
        lons = file['lon'].values
        lats = file['lat'].values
        years = file['years'].values

        nlons = len(lons)
        nlats = len(lats)
        #print "nlons: ",nlons
        #print "nlats: ",nlats

        nlonsnew = math.floor(nlons/sumlon)
        nlatsnew = math.floor(nlats/sumlat)

        nlats2 = np.int(nlatsnew * sumlat)
        nlons2 = np.int(nlonsnew * sumlon)

        Latsnew = np.zeros(nlatsnew,np.float)
        Lonsnew = np.zeros(nlonsnew,np.float)

        # Create Outfile structure
        ncfile = Dataset(DirO + FileO, 'w')
        ncfile.createDimension('year', nyears) 
        ncfile.createDimension('size', 4)
        ncfile.createDimension('lon', nlonsnew)
        ncfile.createDimension('lat', nlatsnew)
        OLon = ncfile.createVariable('lon','f4',('lon'),fill_value=-9999)
        OLat = ncfile.createVariable('lat','f4',('lat'),fill_value=-9999)
        OLongitude = ncfile.createVariable('longitude','f4',('lon'),fill_value=-9999)
        OLatitude = ncfile.createVariable('latitude','f4',('lat'),fill_value=-9999)
        OYears = ncfile.createVariable('year','f4',('year'),fill_value=-9999)
        OYears[:] = years
    
        for variable in file.variables.keys(): 
            if variable in ['time','years','lat','lon','season','Longitude','Latitude','latitude','longitude']:
                continue
                print 'not doing anything for dimensions variables'
            else:
                VarMap = file[variable].values

                # Catch all silly values!
                VarMap[np.where(VarMap > 10E30)] = np.nan
                VarMap[np.where(VarMap == -9999)] = np.nan
                VarMap[np.where(VarMap < -10E30)] = np.nan
                

                VarMapnew = np.zeros((nyears,nlatsnew,nlonsnew),np.float64)
                VarMapnew[...] = np.nan

                #'regrid' by summing over large boxes (averaging lons and lats)
                inlat = 0
                for ilats in range(0,nlats2,sumlat):
                    inlon = 0
                    Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlat])
                #   print np.mean(lats[ilats:ilats+sumlat])
                    for ilons in range(0,nlons2,sumlon):
                        Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlon])
                        if variable in sumvars:
                            VarMapnew[:,inlat,inlon] = np.nansum(VarMap[:,ilats:ilats+sumlat,ilons:ilons+sumlon],axis=(1,2),dtype=np.float64)   
                        else:
                            VarMapnew[:,inlat,inlon] = np.nanmean(VarMap[:,ilats:ilats+sumlat,ilons:ilons+sumlon],axis=(1,2),dtype=np.float64)         
        
                        inlon += 1
                    inlat += 1
                
                OVarMap = ncfile.createVariable(variable,'f4',('year','lat','lon'),fill_value=-9999)
                
                if Latsnew[0] > Latsnew[-1]:
                    OVarMap[:,:,:] = VarMapnew[:,::-1,:]
                else:
                    OVarMap[:,:,:] = VarMapnew

        # After last variable, write lats and lons, which will be identical for all variables and just need to be written once
        if Latsnew[0] > Latsnew[-1]:
            OLatitude[:] = Latsnew[::-1]
            OLat[:] = Latsnew[::-1]
        else:
            OLatitude[:] = Latsnew
            OLat[:] = Latsnew
        OLongitude[:] = Lonsnew
        OLon[:] = Lonsnew

        ncfile.close()

Ngl.end()



