# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:45:25 2015
Code to plot figures of precipitation events density and total rainfall using command line arguments
Example use:
python paperplot_compare_arg.py --Data1 TRMMERAIgd --Version1 Standard
--anstartyr 1998 --anendyr 2014 --Data2 ERAI --Version2 Standard --anstartyr2
--tbound1 0 1 2 5 --tbound2 1 2 5 100 --unit day --splittype day --minlat -40
--maxlat 40 --sumlons 3 --sumlats 3 --minGB 9
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
import argparse
import resource

rsrcV = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit starts as  :', soft
print 'Hard limit starts as :', hard
resource.setrlimit(rsrcV, (50000000000, hard)) #limit memory usage
#                          137438953472
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit changed to :', soft

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

parser.add_argument('--plotvar',type=str,nargs=1,help='variable to plot')

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

plotvar = args.plotvar[0]

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

def getplotdata(Datain,Versionin,
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
    PrecipIn = getrawPrecipAnn(Datain,Versionin,MinLatF,MaxLatF,
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
    if FileIn['lat'][0] > FileIn['lat'][1]:
        switchlats = True
        lats = FileIn['lat'][::-1].sel(lat=slice(MinLatF,MaxLatF)).values
    else:
        switchlats = False
        lats = FileIn['lat'].sel(lat=slice(MinLatF,MaxLatF)).values

    lons = FileIn['lon'].values

    try:
        years = FileIn['years'].sel(years=slice(anstartyr,anendyr))
    except KeyError:    # years is called year
        years = FileIn['year'].sel(year=slice(anstartyr,anendyr))

    nlats = len(lats)
    nlons = len(lons)
    nyears = len(years)

    difflat = lats[1]-lats[0]
    difflon = lons[1]-lons[0]

    if args.totals == 1:
        arraynbounds = nbounds + 1
    else:
        arraynbounds = nbounds

    data = np.zeros([arraynbounds,nyears,nlats,nlons])      # +1 for total sum as first plot

    # read in data for each category.
    for iday in range(0,nbounds):
        FileIn = ('DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-'
                    + str(tbound2[iday]) + unit + '_' + mapping + '_' + Datain + '_'
                    + str(Fstartyrin) + '-' + str(Fendyrin) + '_' + Versionin + fileadd
                    + '.nc')

        #Get lons and lats
        print DirIn + '/' + diradd + '/' + FileIn
        FileIn = xray.open_dataset(DirIn + '/' + diradd + '/' + FileIn)

        try:
            if switchlats:
                data[iday,:,:,:] = (FileIn[plotvar][:,::-1,:]
                                        .sel(lat=slice(MinLatF,MaxLatF),
                                             years=slice(anstartyr,anendyr)))
            else:
                data[iday,:,:,:] = (FileIn[plotvar].sel(lat=slice(MinLatF,MaxLatF),
                                                        years=slice(anstartyr,anendyr)))

        except KeyError:    # time is year, not years
            if switchlats:
                data[iday,:,:,:] = (FileIn[plotvar][:,::-1,:]
                                        .sel(lat=slice(MinLatF,MaxLatF),
                                             year=slice(anstartyr,anendyr)))
            else:
                data[iday,:,:,:] = (FileIn[plotvar]
                                        .sel(lat=slice(MinLatF,MaxLatF),
                                             year=slice(anstartyr,anendyr)))

        # Make time means
        dataAnnAvg = np.nanmean(data,axis=1)

        # Make event lifetime means
        dataAllAvg = np.nansum(dataAnnAvg,axis=0)

        # Calculate percentages
        if plotvar == 'TDensity':
            # get values per square degree
            dataAllPercent = dataAllAvg  / (difflat * difflon) 
            dataSum =  np.nansum(dataAllAvg)
        elif plotvar == "TPrecip":
            dataAllPercent = 100.0 * np.divide(dataAllAvg,PrecipIn*365.0)  # convert PrecipIn from mm/day to mm/yr before dividing
            dataSum = 100.0 * np.nansum(dataAllAvg)/np.nansum(PrecipIn * 365.0)

        dataPercent = np.zeros(dataAnnAvg.shape)
        titlesdata = []

        # Fill percent arrays: first is total
        arrayindex = 0
        if args.totals == 1:
            dataPercent[0,:,:] = dataAllPercent

            if plotvar == 'TDensity':
                titlesdata.append('Total annual event density, events/ ' + 
                                  '(yr/deg~S1~2); mean = ' +
                                  '{:2.3g}'.format(dataSum))
                titlemain = 'density'
            elif plotvar == 'TPrecip':
                titlesdata.append('Percentage of total precipitation ' + 
                                  'captured in events; mean = ' +
                                  '{:2.3g}'.format(dataSum) + '%')
                titlemain = 'precip'
            arrayindex += 1

        for iday in range(0,nbounds):
            dataPercent[iday+arrayindex,:,:] = 100.0 * (np.where(dataAllAvg > 0,
                                                        np.divide(dataAnnAvg[iday,:,:],dataAllAvg),0.0)) 

            if tbound2[iday] < tbound1[iday]*10 or tbound1[iday] == 0 or tbound2[iday] == 0:
                titlesdata.append(str(int(tbound1[iday])) + 
                                  ' to ' + str(int(tbound2[iday])) + ' ' + unit + 
                                  ' events; mean = ' +
                                  '{:2.3g}'.format(np.nanmean(dataPercent[iday+arrayindex,:,:])))
            else:
                titlesdata.append('>' + str(int(tbound1[iday])) + ' ' +
                                  unit + ' events; mean = ' +
                                  '{:2.3g}'.format(np.nanmean(dataPercent[iday+arrayindex,:,:])))

    return(titlesdata,dataPercent,lats,lons)

# Get the data

titlescol1,datacol1,lats,lons = getplotdata(Data1,Version1,
                                  MinLatF,MaxLatF,
                                  anstartyr,anendyr)

titlescol2,datacol2,lats,lons = getplotdata(Data2,Version2,
                                  MinLatF,MaxLatF,
                                  anstartyr,anendyr)

# And now plot 
# Plot 1: Average density, and percentage easterly

figtitlein = (FigDir + 'Paper_' + plotvar + Data1 + '_' + Version1 + '_' + Data2 +
             '_' + Version2 + str(anstartyr) + '-' + str(anendyr) + '_' +
              str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit +
              fileOadd)

titlein = ('Events in ' + Data1 + " " + Version1 + '(L) and ' + Data2 + ' ' +
            Version2 + '(R), years ' + str(anstartyr) + '-' + str(anendyr))

if splittype == 'day':
    if Data1 in ["TRMM"]:
        clims1,clims2,clims3,clims4 = [0.0,98.0,0.0,0.0,0.0,0.0],[300,100.0,1.0,0.2,0.05,1.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,90.0,40.0,40.0,40.0,40.0]  
    elif Data1 in ["TRMMERAIgd"]:
        if minGB in [4,9]:
            clims1,clims2,clims3,clims4 = [0.0,70.0,0.0,0.0,0.0,0.0],[6,100.0,20.0,7.0,4.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,50.0,50.0,50.0,50.0]

        else:
            clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[10,100.0,15.0,5.0,2.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,40.0,50.0,50.0,50.0]
    elif Data1 in ["ERAI"]:
        if minGB in [4,9]:
            clims1,clims2,clims3,clims4 = [0.0,70.0,0.0,0.0,0.0,0.0],[6,100.0,20.0,7.0,4.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,50.0,50.0,50.0,50.0]
        else:
            clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[10,100.0,15.0,5.0,2.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,40.0,50.0,50.0,50.0]

#       clims1,clims2,clims3,clims4 = [0.0,95.0,0.0,0.0,0.0,0.0],[30,100.0,2.0,1.5,0.5,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]
    elif Data1 in ["ERA20C"]:
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[15,100.0,4.0,3.0,2.0,4.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]        
    elif Data1 in ["CESM"]:
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0],[50,100.0,5.0,3.0,0.5],[0.0,0.0,0.0,0.0,0.0],[100,90.0,40.0,40.0,40.0]  

elif splittype == 'maxspeed':
    if Data1 == "TRMM":
        clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,95.0,0.0,0.0,0.0],[0.5,0.75,5.0,100.0,5.0,0.75,0.5],[0.0,0.0,0.0,20.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]
    elif Data1 in ["TRMMERAIgd"]:
        clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,90.0,0.0,0.0,0.0],[0.5,1.0,5.0,100.0,5.0,1.0,0.5],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]
    elif Data1 in ["ERAI"]:
        clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,90.0,0.0,0.0,0.0],[0.5,1.0,5.0,100.0,5.0,1.0,0.5],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]
#               clims1,clims2,clims3,clims4 = [0.0,95.0,0.0,0.0,0.0,0.0],[30,100.0,2.0,1.5,0.5,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]

    elif Data1 == "CESM":
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0],[50,100.0,5.0,3.0,0.5],[0.0,0.0,0.0,0.0,0.0],[20,90.0,40.0,40.0,40.0]

    elif Data1 in ["ERA20C"]:
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[15,100.0,4.0,3.0,2.0,4.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]


if plotvar == 'TDensity':
    datacol1 = conserveRegrid(datacol1,
                                  'lat','lon',
                                  lats,lons,
                                  sumlats,sumlons)
    datacol2 = conserveRegrid(datacol2,
                                  'lat','lon',
                                  lats,lons,
                                  sumlats,sumlons)

    plotmap(datacol1,datacol2,
            clims1,clims2,clims1,clims2,
            titlescol1,titlescol2,
            titlein,figtitlein,
            lons,lats,
            MinLonF,MaxLonF,MinLatF,MaxLatF,
            FillValue)

elif plotvar == 'TPrecip':
    plotmap(datacol1,datacol2,
            clims3,clims4,clims3,clims4,
            titlescol1,titlescol2,
            titlein,figtitlein,
            lons,lats,
            MinLonF,MaxLonF,MinLatF,MaxLatF,
            FillValue)

Ngl.end()




