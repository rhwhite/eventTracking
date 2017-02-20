# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:45:25 2015
Code to plot figures of precipitation events density and total rainfall using command line arguments
Example use:
python paperplot_compare_arg.py --Data1 TRMMERAIgd --Version1 Standard
--anstartyr 1998 --anendyr 2014 --Data2 ERAI --Version2 Standard
 --tbound1 0 1 2 5 --tbound2 1 2 5 100 --unit day --splittype day --minlat -40
 --maxlat 40 --sumlons 3 --sumlats 3 --minGB 9 --plotvar1 TPrecip --plotvar2
 TPrecip
@author: rachel
"""

import os, errno
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray
import Ngl
import math
from scipy import stats
from rhwhitepackages.readwrite import *
from rhwhitepackages.stats import regressmaps
from rhwhitepackages.plots import *
from rhwhitepackages.regrid import *
import argparse
import resource

rsrcV = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrcV)
#print 'Soft limit starts as  :', soft
#print 'Hard limit starts as :', hard
resource.setrlimit(rsrcV, (50000000000, hard)) #limit memory usage
#                          137438953472
soft, hard = resource.getrlimit(rsrcV)
#print 'Soft limit changed to :', soft

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs='?',default=4,help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound1',type=float,nargs='+',help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=float,nargs="+",help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')
parser.add_argument('--Data1',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version1',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--Data2',type=str,nargs=1,help='type of Data, TRMM, ERAI,or CESM')
parser.add_argument('--Version2',type=str,nargs=1,help='Version of Data,Standard, low, 6th_from6 etc')
parser.add_argument('--anstartyr',type=int,nargs=1,help='start year foranalysis')
parser.add_argument('--anendyr',type=int,nargs=1,help='end year for analysis')

parser.add_argument('--plotvar1',type=str,nargs=1,help='variable to plot')
parser.add_argument('--plotvar2',type=str,nargs=1,help='variable to plot')

parser.add_argument('--minlat',type=int,nargs='?',default=-45,help='min lat')
parser.add_argument('--maxlat',type=int,nargs='?',default=45,help='max lat')
parser.add_argument('--minlon',type=int,nargs='?',default=0,help='min lon')
parser.add_argument('--maxlon',type=int,nargs='?',default=360,help='max lon')
parser.add_argument('--sumlats',type=int,nargs='?',default=-1,help='number of lats to sum over')
parser.add_argument('--sumlons',type=int,nargs='?',default=-1,help='number of lons to sum over')
parser.add_argument('--minGB',type=int,nargs='?',default=0,help='minimum number of gridboxes to count as an event')

parser.add_argument('--totals',type=int,nargs='?',default=0,
                    help='plot totals or not')


args = parser.parse_args()

print('here\'s what I have as arguments: ', args)

splittype = args.splittype[0]
speedtspan = args.speedtspan
tbound1 = args.tbound1
tbound2 = args.tbound2
unit = args.unit[0]
Data1 = args.Data1[0]
Version1 = args.Version1[0]
Data2 = args.Data2[0]
Version2 = args.Version2[0]
anstartyr = args.anstartyr[0]
anendyr = args.anendyr[0]
MinLonF = args.minlon
MaxLonF = args.maxlon
MinLatF = args.minlat
MaxLatF = args.maxlat
sumlats = args.sumlats
sumlons = args.sumlons
minGB = args.minGB

plotvar1 = args.plotvar1[0]
plotvar2 = args.plotvar2[0]

if splittype=='day':
    splitname = 'Sizes'
elif splittype =='maxspeed':
    splitname = 'MaxSpeeds_' + str(speedtspan) + 'ts'

nbounds = len(tbound1)

FillValue = -9999

# Time period for analysis
nyears = anendyr - anstartyr + 1
mapping = "center"
FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

base = '/home/disk/eos4/rachel/EventTracking/'

if minGB > 0:
    fileadd = '_min' + str(minGB) + 'GB'
else:
    fileadd = ''

fileOadd = fileadd
if sumlons > 0:
    fileOadd = fileOadd + '_' + str(sumlons) + 'sumlons'

# Define functions

def getplotdata(plotvarin,Datain,Versionin,
            MinLatF,MaxLatF,
            anstartyrin,anendyrin):

    # make assumptions on file start years
    if Datain in ['TRMMERAIgd','TRMM']:
        Fstartyrin = 1998
        Fendyrin = 2014
    elif Datain in ['ERAI']:
        Fstartyrin = 1980
        Fendyrin = 2014

    # get precip data for dataset
    # Add +/- 5 lats for regridding
    PrecipIn = getrawPrecipAnn(Datain,Versionin,MinLatF-5,MaxLatF+5,
                                anstartyrin,anendyrin)

    # take time mean
    PrecipIn = PrecipIn.mean(dim='time')
    # get data to plot 
    iday = 0

    DirIn = (base + '/FiT_RW_ERA/' + Datain + '_output/' + Versionin +
                 str(Fstartyrin) + '/proc/')

    FileIn = ('DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-' +
                str(tbound2[iday]) + unit + '_' + mapping + '_' + Datain + '_' + 
                str(Fstartyrin) + '-' + str(Fendyrin) + '_' + Versionin + fileadd +  '.nc')

    if splittype == "maxspeed":
        diradd = "MaxSpeeds"
    elif splittype == "speed":
        diradd = "Speeds"
    elif splittype == "day":
        diradd = "Sizes"
    else:
        exit("unexpected splittype")

    FileIn = xrayOpen(DirIn + '/' + diradd + '/' + FileIn)

    # Always read in file to go from South to North
    # Add 5 lats in either direction for regridded averages
    if FileIn['lat'][0] > FileIn['lat'][1]:
        switchlats = True
        lats = FileIn['lat'][::-1].sel(lat=slice(MinLatF-5,MaxLatF+5)).values
    else:
        switchlats = False
        lats = FileIn['lat'].sel(lat=slice(MinLatF-5,MaxLatF+5)).values

    lons = FileIn['lon'].values

    yearstr = 0
    try:
        years = FileIn['years'].sel(years=slice(anstartyr,anendyr)).values
    except KeyError:    # years is called year
        try:
            years = FileIn['year'].sel(year=slice(anstartyr,anendyr)).values
        except KeyError:    # time is 'time'
            years = FileIn['time'].sel(time=slice(
                                        str(anstartyr),str(anendyr)))
            yearstr = 1

    nlats = len(lats)
    nlons = len(lons)
    nyears = len(years)

    difflat = lats[1]-lats[0]
    difflon = lons[1]-lons[0]

    data = np.zeros([nbounds,nyears,nlats,nlons])      # +1 for total sum as first plot

    # read in data for each category.
    for iday in range(0,nbounds):
        FileIn = ('DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-'
                    + str(tbound2[iday]) + unit + '_' + mapping + '_' + Datain + '_'
                    + str(Fstartyrin) + '-' + str(Fendyrin) + '_' + Versionin + fileadd
                    + '.nc')

        #Get lons and lats
        FileIn = xrayOpen(DirIn + '/' + diradd + '/' + FileIn)

        # Adding +/- 5 lats for regridding purposes
        try:
            if switchlats:
                data[iday,:,:,:] = (FileIn[plotvarin][:,::-1,:]
                                        .sel(lat=slice(MinLatF-5,MaxLatF+5),
                                             years=slice(anstartyr,anendyr)))
            else:
                data[iday,:,:,:] = (FileIn[plotvarin].sel
                                                    (lat=slice(MinLatF-5,MaxLatF+5),
                                                      years=slice(anstartyr,anendyr)))

        except KeyError:    # time is year, not years
            if yearstr == 0:
                if switchlats:
                    data[iday,:,:,:] = (FileIn[plotvarin][:,::-1,:]
                                            .sel(lat=slice(MinLatF-5,MaxLatF+5),
                                                 year=slice(anstartyr,anendyr)))
                else:
                    data[iday,:,:,:] = (FileIn[plotvarin]
                                            .sel(lat=slice(MinLatF-5,MaxLatF+5),
                                                 year=slice(anstartyr,anendyr)))
            elif yearstr == 1:
                if switchlats:
                    data[iday,:,:,:] = (FileIn[plotvarin][:,::-1,:]
                                            .sel(lat=slice(MinLatF-5,MaxLatF+5),
                                                 time=slice(str(anstartyr)
                                                            ,str(anendyr))))
                else:
                    data[iday,:,:,:] = (FileIn[plotvarin]
                                            .sel(lat=slice(MinLatF-5,MaxLatF+5),
                                                 time=slice(str(anstartyr)
                                                            ,str(anendyr))))
    if args.totals == 1:
        lifetimemin =([-100] + tbound1)
    else:
        lifetimemin = tbound1

    data = xray.DataArray(data,
                           coords=[('lifetimemin',tbound1),
                                   ('year',years),
                                   ('lat',lats),
                                   ('lon',lons)])

    # Make time means
    dataAnnAvg = data.mean(dim='year') #np,nanmean(data,axis=1)
    # Make event lifetime means
    dataAllAvg = dataAnnAvg.sum(dim='lifetimemin')

    # Calculate percentages
    if plotvarin == 'TDensity':
        # get values per square degree
        dataAllPercent = (dataAllAvg /
                          (difflat * difflon)) 
        # get total value per year per square degree
        dataSum = (np.nansum(dataAllAvg.sel(
                                    lat=slice(MinLatF,MaxLatF)))/
                          ((MaxLatF - MinLatF) * (MaxLonF - MinLonF))) 
    elif plotvarin == "TPrecip":
        # convert PrecipIn from mm/day to mm/yr before dividing
        dataAllPercent = 100.0 * np.divide(dataAllAvg,
                                           PrecipIn*365.0)
        try:
            dataSum = 100.0*(np.nansum(dataAllAvg.sel(lat=slice(MinLatF,MaxLatF))/
                        np.nansum(PrecipIn.sel(lat=slice(MinLatF,MaxLatF)) * 365.0)))
        except KeyError:
            dataSum = 100.0 * (
                      np.nansum(dataAllAvg.sel(lat=slice(MinLatF,MaxLatF))/
                      np.nansum(PrecipIn.sel(latitude=slice(MinLatF,MaxLatF)) *
                                365.0)))

    if args.totals == 1:
        arraynbounds = nbounds + 1
    else:
        arraynbounds = nbounds

    dataPercent = np.zeros([arraynbounds,nlats,nlons]) 
    titlesdata = []

    # Fill percent arrays: first is total
    arrayindex = 0
    if args.totals == 1:
        dataPercent[0,:,:] = dataAllPercent
        if dataSum > 100:
            strtotal = '{:4.3g}'.format(dataSum)
        elif dataSum > 10:
            strtotal = '{:2.2g}'.format(dataSum)
        else:
            strtotal = '{:2.1g}'.format(dataSum)

        if plotvarin == 'TDensity':
            titlesdata.append('Density, events/' + 
                              'yr/deg~S1~2; mean = ' +
                              strtotal)
            titlemain = 'density'
        elif plotvarin == 'TPrecip':
            titlesdata.append('Precipitation; ' + 
                              ' mean = ' +
                              strtotal + '%')
            titlemain = 'precip'
        arrayindex += 1

    for iday in range(0,nbounds):
        dataPercent[iday+arrayindex,:,:] = 100.0 * (np.where(dataAllAvg > 0,
                                                    np.divide(dataAnnAvg[iday,:,:],dataAllAvg),0.0)) 

        if np.nanmean(dataPercent[iday+arrayindex,:,:]) > 100:
            strmean = '{:4.3g}'.format(
                            np.nanmean(dataPercent[iday+arrayindex,:,:]))
        elif np.nanmean(dataPercent[iday+arrayindex,:,:]) > 10:
            strmean = '{:2.2g}'.format(
                            np.nanmean(dataPercent[iday+arrayindex,:,:]))
        else:
            strmean = '{:2.1g}'.format(
                            np.nanmean(dataPercent[iday+arrayindex,:,:]))

        if (abs(tbound2[iday]) < abs(tbound1[iday]*10) or 
            abs(tbound1[iday]) < abs(tbound2[iday]*10) or
            tbound1[iday] == 0 or tbound2[iday] == 0):
            titlesdata.append(str(int(tbound1[iday])) + 
                              '-' + str(int(tbound2[iday])) + ' ' + unit + 
                              '; mean = ' +
                              strmean + '%')
        else:
            titlesdata.append('>' + str(int(tbound1[iday])) + ' ' +
                              unit + '; mean = ' +
                              strmean + '%')
    newDA = xray.DataArray(dataPercent,
                           coords=[('lifetimemin',lifetimemin),
                                   ('lat',lats),
                                   ('lon',lons)])
    return(titlesdata,newDA)

# Get the data

titlescol1,datacol1 = getplotdata(plotvar1,
                                  Data1,Version1,
                                  MinLatF,MaxLatF,
                                  anstartyr,anendyr)

titlescol2,datacol2 = getplotdata(plotvar2,
                                  Data2,Version2,
                                  MinLatF,MaxLatF,
                                  anstartyr,anendyr)

# And now plot 
# Plot 1: Average density, and percentage easterly

figtitlein = (FigDir + 'Paper_' + plotvar1 + '_' + Data1 + '_' + Version1 +
             '_' + plotvar2 + '_' + Data2 +
             '_' + Version2 + str(anstartyr) + '-' + str(anendyr) + '_' +
              str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit +
              fileOadd)

titlein = ('Events in ' + Data1 + " " + Version1 + '(L) and ' + Data2 + ' ' +
            Version2 + '(R), years ' + str(anstartyr) + '-' + str(anendyr))


cb1,cb2 = getFITcolorbars(Data1,minGB,splittype,plotvar1)
cb3,cb4 = getFITcolorbars(Data2,minGB,splittype,plotvar2)

if plotvar1 in ['TDensity','TPrecip']:
    datacol1 = conserveRegrid(datacol1,
                                  'lat','lon',
                                  sumlats,sumlons)

if plotvar2 in ['TDensity','TPrecip']:
    datacol2 = conserveRegrid(datacol2,
                                  'lat','lon',
                                  sumlats,sumlons)
if args.totals == 1:
    print('plotting totals')
    plotmap(datacol1,datacol2,
            cb1,cb2,cb3,cb4,
            titlescol1,titlescol2,
            titlein,figtitlein,
            MinLonF,MaxLonF,MinLatF,MaxLatF,
            FillValue)
else:
    # if not plotting totals, take 1st cb values out of array
    plotmap(datacol1,datacol2,
            cb1[1:-1],cb2[1:-1],cb3[1:-1],cb4[1:-1],
            titlescol1,titlescol2,
            titlein,figtitlein,
            MinLonF,MaxLonF,MinLatF,MaxLatF,
            FillValue)




