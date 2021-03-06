"""
Code to plot figures of precipitation events density and total rainfall using command line arguments

Author: Rachel White, rachel.white@cantab.net

Created: July 2016

Example use:
python paperplot_compare_arg.py --Data1 TRMMERAIgd --Version1 Standard \
--anstartyr 1998 --anendyr 2014 --Data2 ERAI --Version2 Standard \
--tbound1 0 1 2 5 --tbound2 1 2 5 100 --unit day --splittype day --minlat -40 \
 --maxlat 40 --sumlons 3 --sumlats 3 --minGB 9 --plotvar1 TPrecip --plotvar2 \
 TPrecip
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
            MinLatF,MaxLatF,MinLonF,MaxLonF,
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

    FileIn = getdenfilename(mapping, Datain, Versionin, Fstartyrin, Fendyrin,iday,
                 splittype, unit, 4, minGB, tbound1, tbound2,
                 'Mon',-1,-1)


    # Use MinLonF and MaxLonF to capture all possible longitudes
    # And then work out which to use
    VarIn = geteventmapdata(iday,'Ann','TDensity',FileIn,
                    1998, 2014,
                    MinLatF-5,MaxLatF+5,-180,360,[])

    # Deal with different possible longitudes: 0 to 360 or -180 to 180
    if VarIn.lon[0] < 0 and VarIn.lon[-1] <= 180:
        if MaxLonF > 180:
            MaxLonF = MaxLonF - 180
            MinLonF = MinLonF - 180

    VarIn = geteventmapdata(iday,'Ann','TDensity',FileIn,
                    1998, 2014,
                    MinLatF-5,MaxLatF+5,-180,360,[])

    # Always read in file to go from South to North
    # Add 5 lats in either direction for regridded averages
    if VarIn.lat[0] > VarIn.lat[1]:
        switchlats = True
        lats = VarIn.lat[::-1].values
    else:
        switchlats = False
        lats = VarIn.lat.values

    lons = VarIn['lon'].values

    yearstr = 0
    try:
        years = VarIn['years'].sel(years=slice(anstartyr,anendyr)).values
    except KeyError:    # years is called year
        try:
            years = VarIn['year'].sel(year=slice(anstartyr,anendyr)).values
        except KeyError:    # time is 'time'
            years = VarIn['time'].sel(time=slice(
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

        FileIn = getdenfilename(mapping, Datain, Versionin, 
                                Fstartyrin, Fendyrin,iday,
                                splittype, unit, 4, minGB, 
                                tbound1, tbound2,
                                'Mon',-1,-1)

        # Adding +/- 5 lats for regridding purposes
        VarIn = geteventmapdata(iday,'Ann',plotvarin,FileIn,
                    anstartyr, anendyr,
                    MinLatF-5,MaxLatF+5,MinLonF,MaxLonF,
                    [])

        if switchlats:
            data[iday,:,:,:] = (VarIn[:,::-1,:].values)
        else:
            data[iday,:,:,:] = (VarIn.values)

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
    elif plotvarin == 'LocalDensity':
        dataAllPercent = dataAllAvg
        dataSum = (np.nansum(dataAllAvg.sel(
                                    lat=slice(MinLatF,MaxLatF))))
    if args.totals == 1:
        arraynbounds = nbounds + 1
    else:
        arraynbounds = nbounds

    dataPercent = np.zeros([arraynbounds,nlats,nlons]) 
    titlesdata = []
    paneltitles = []
    labeltitles = []

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
            titlesdata.append('')
            labeltitles.append('events/yr/deg~S1~2')
            #titlesdata.append('Density, events/' + 
            #                  'yr/deg~S1~2; mean = ' +
            #                  strtotal)
            titlemain = 'density'
        elif plotvarin == 'TPrecip':
            labeltitles.append('%')

            #titlesdata.append('Precipitation; ' + 
            #                  ' mean = ' +
            #                  strtotal + '%')
            titlesdata.append('')
            titlemain = 'precip'
        elif plotvarin == 'LocalDensity':
            labeltitles.append('events/yr')
            #titlesdata.append('Density, events/yr')
            titlesdata.append('')
        paneltitles.append('Total') 

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
            strmean = '{:2.2g}'.format(
                            np.nanmean(dataPercent[iday+arrayindex,:,:]))

        if (iday == 0):
            paneltitles.append('<' + str(int(tbound2[iday])) + ' ' + unit) 
            titlesdata.append('')
            labeltitles.append('%')
        elif (iday == nbounds - 1):
            paneltitles.append('>' + str(int(tbound1[iday])) + ' ' + unit) 
            titlesdata.append('')
            labeltitles.append('%')
        else:
            paneltitles.append(str(int(tbound1[iday])) + 
                              '-' + str(int(tbound2[iday])) + ' ' + unit)
            titlesdata.append('')
                              #str(int(tbound1[iday])) + 
                              #'-' + str(int(tbound2[iday])) + ' ' + unit + 
                              #'; mean = ' +
                              #strmean + '%')
            labeltitles.append('%')
    newDA = xray.DataArray(dataPercent,
                           coords=[('lifetimemin',lifetimemin),
                                   ('lat',lats),
                                   ('lon',lons)])
    return(titlesdata,paneltitles,labeltitles,newDA)

# Get the data

titlescol1,paneltitles,labeltitles1,datacol1 = getplotdata(plotvar1,
                                  Data1,Version1,
                                  MinLatF,MaxLatF,
                                  MinLonF,MaxLonF,
                                  anstartyr,anendyr)

titlescol2,paneltitles,labeltitles2,datacol2 = getplotdata(plotvar2,
                                  Data2,Version2,
                                  MinLatF,MaxLatF,
                                  MinLonF,MaxLonF,
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

if plotvar1 in ['LocalDensity','TDensity','TPrecip'] and sumlats > 0:
    datacol1 = conserveRegrid(datacol1,
                                  'lat','lon',
                                  sumlats,sumlons)

if plotvar2 in ['LocalDensity','TDensity','TPrecip'] and sumlats > 0:
    datacol2 = conserveRegrid(datacol2,
                                  'lat','lon',
                                  sumlats,sumlons)
print paneltitles
if args.totals == 1:
    print('plotting totals')
    plotmap(datacol1,datacol2,
            cb1,cb2,cb3,cb4,
            titlescol1,titlescol2,
            titlein,figtitlein,
            MinLonF,MaxLonF,MinLatF,MaxLatF,
            FillValue,paneltitles,
            labeltitles1,labeltitles2)
else:
    # if not plotting totals, take 1st cb values out of array
    plotmap(datacol1,datacol2,
            cb1[1:-1],cb2[1:-1],cb3[1:-1],cb4[1:-1],
            titlescol1,titlescol2,
            titlein,figtitlein,
            MinLonF,MaxLonF,MinLatF,MaxLatF,
            FillValue,paneltitles,
            labeltitles1,labeltitles2)




