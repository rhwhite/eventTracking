"""
Code to sum characteristics of events split by either time or speed and write
out the results to a table.

Author: Rachel H White rhwhite@uw.edu
Created: Oct 2016

Example use:
python CalculateDomainAverages.py --Data TRMM --Version Standard \
--anstartyr 1998 --anendyr 2014 \
--tbound1 0 1 2 5 --tbound2 1 2 5 100 --splittype day \
--unit day --minlat -40 --maxlat 40
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

parser = argparse.ArgumentParser(description="Calculate domain averages")
parser.add_argument('--minlat',type=int,nargs='?',default=-45,help='min lat')
parser.add_argument('--maxlat',type=int,nargs='?',default=45,help='max lat')
parser.add_argument('--minlon',type=int,nargs='?',default=0,help='min lon')
parser.add_argument('--maxlon',type=int,nargs='?',default=360,help='max lon')
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,
                    help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs='?',default=4,
                    help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound1',type=float,nargs='+',
                    help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=float,nargs="+",
                    help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')
parser.add_argument('--Data',type=str,nargs=1,
                    help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,
                    help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--anstartyr',type=int,nargs=1,
                    help='start year for analysis')
parser.add_argument('--anendyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--test',type=int,nargs='?',default=0,help='1 for test')
parser.add_argument('--extra',type=int,nargs='?',default=0,help='1 for test')
parser.add_argument('--tperday',type=int,nargs='?',default=8,help=
        'timesteps per day, default is 3hourly = 8')
parser.add_argument('--minGB',type=int,nargs=1,default=0,help='min gridboxes')


args = parser.parse_args()

print(str(args))
print(str(dt.date.today()))

# put inputs into the type and variable names we want
splittype = args.splittype[0]
speedtspan = args.speedtspan
# multiply tbound in by 24 to get hours rather than days
tbound1 = np.multiply(args.tbound1,24.0)
tbound2 = np.multiply(args.tbound2,24.0)
unit = args.unit[0]
Data = args.Data[0]
Version = args.Version[0]
anstartyr = args.anstartyr[0]
anendyr = args.anendyr[0]
minlon = args.minlon
maxlon = args.maxlon
minlat = args.minlat
maxlat = args.maxlat
test = args.test
extra = args.extra
minGB = args.minGB[0]

diradd = getdirectory(splittype)
nbounds = len(tbound1)

print(tbound1)

print minGB
print str(minGB)

R = 6371000     # radius of Earth in m

nyears = anendyr - anstartyr + 1

minevent = 100000

if Data == "TRMM":
    FileInLats = ('/home/disk/eos4/rachel/Obs/TRMM/'
                    'SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc')
    Fstartyr = 1998
    Fendyr = 2014

elif Data == "TRMMERAIgd":
    FileInLats = ('/home/disk/eos4/rachel/Obs/TRMM/'
                    'regrid2ERAI_TRMM_3B42_1998-2014.nc')
    Fstartyr = 1998
    Fendyr = 2014

elif Data == "ERAI":
    FileInLats = ('/home/disk/eos4/rachel/Obs/ERAI/3hrly/Precip_3hrly/'
                'SeasAnn_ERAI_Totalprecip_1980-2014_preprocess.nc')
    Fstartyr = 1980
    Fendyr = 2014

elif Data == "ERA20C":
    FileInLats = '/home/disk/eos4/rachel/Obs/ERA_20C/ERA_20C_LatLon.nc'

elif Data == "CESM":
    FileInLats = ('/home/disk/eos4/rachel/EventTracking/Inputs/CESM/'
                    'f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
                    'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc')
else:
    print("unexpected data type")
    exit()

DirI = ('/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/'
            + Version + str(Fstartyr) + '/proc/')

DirO = DirI + diradd + '/'

# Work out how many to skip at beginning and end of file
filedate = dt.date(Fstartyr,1,1)
startdate = dt.date(anstartyr,1,1)
enddate = dt.date(anendyr+1,1,1)

if Fstartyr == anstartyr:
    starttstart = 0
elif Fstartyr > anstartyr:
    print("Fstartyr is ", Fstartyr, " and anstartyr is ",anstartyr)
    exit("cannot start analysing data before there is data!")
elif Fstartyr < anstartyr:
    diffdays = (startdate-filedate).days
    starttstart = diffdays * args.tperday

endtstart = starttstart + (enddate - startdate).days * args.tperday

# In[4]:

#Get lons and lats
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

# initialize data
averageydist = np.zeros([nbounds],float)
averagexdist = np.zeros([nbounds],float)
averageprecipperhr = np.zeros([nbounds],float)
averageprecipperareahr = np.zeros([nbounds],float)
averagetime = np.zeros([nbounds],float)
averagem2 = np.zeros([nbounds],float)
averagegridboxes = np.zeros([nbounds],float)
precipvolume = np.zeros([nbounds],float)

count = np.zeros([nbounds],int)


# In[5]:

# open main dataset and read in data
FileI1 = ('All_Precip_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Data + '_' +
            Version + '.nc')

datain = xrayOpen(DirI + FileI1,decodetimes=False)

nevents = len(datain.events)

ycenterstart = datain.ycenterstart[0:nevents].values
xcenterstart = datain.xcenterstart[0:nevents].values

ycenterend = datain.ycenterend[0:nevents].values
xcenterend = datain.xcenterend[0:nevents].values


ycentermean = datain.ycentermean[0:nevents].values
xcentermean = datain.xcentermean[0:nevents].values


tstart = datain.tstart[0:nevents].values
timespan = datain.timespan[0:nevents].values
totalprecip = datain.totalprecip[0:nevents].values
totalprecipSA = datain.totalprecipSA[0:nevents].values
gridboxspan = datain.gridboxspan[0:nevents].values
gridboxspanSA = datain.gridboxspanSA[0:nevents].values
unigridboxspan = datain.uniquegridboxspan[0:nevents].values


# In[6]:

# Set fileminlat and filemaxlat if we ran FiT on a subsection of the data
fileminlat = -90
filemaxlat = 90

# In[ ]:

def isinregion(ilat,ilon):
    if ilat < minlat or ilat > maxlat:
        return(False)
    else:
        if checklons:
            if ilon < minlon or ilon > maxlon:
                return False
        else:
            return(True)


# In[ ]:
checklons = True
if minlon == 0 and maxlon == 360:
    checklons = False
elif minlon ==-180 and maxlon == 180:
    checklons = False

if test == 1:
    nevents = 10000
    filenameadd = "test_"
else:
    nevents = len(datain.events)
    filenameadd = ""

for ievent in range(0,nevents):
    if (ievent % 100000 == 0):
        print "ievent: " + str(ievent)
    # check is in timeframe of analysis
    # assuming that tstart = 0 is equal to Fstartyr

    if unigridboxspan[ievent] > minGB:    # if unique
                                                  # gridboxes exceed specified threshold

        if tstart[ievent] >= starttstart and tstart[ievent] <= endtstart:

            # check if in region
            if isinregion(lats[ycenterstart[ievent]],
                          lons[xcenterstart[ievent]]):

                if isinregion(lats[ycenterend[ievent]],
                              lons[xcenterend[ievent]]):

                    for ibound in range(0,nbounds):
                        if timespan[ievent] < tbound2[ibound]:
                            if extra == 1:
                                averageydist[ibound] += (lats[ycenterend[ievent]] -
                                                        lats[ycenterstart[ievent]])
                                averagexdist[ibound] += (lons[xcenterend[ievent]] -
                                                        lons[xcenterstart[ievent]])
                            # if negative then that's fine for NH, positive is fine
                            # for southern hemisphere, so get the average event
                            # distance travelled

                            averagetime[ibound] += timespan[ievent]
                            averageprecipperhr[ibound] += (
                                                            totalprecip[ievent]/timespan[ievent])
                            # Include factor of 3 to convert to hours, not timesteps
                            averageprecipperareahr[ibound] += (
                                                                totalprecip[ievent]/(3.0 *
                                                                gridboxspan[ievent]))
                            # Include factor of 3 because timespan is in hours, not timesteps
                            averagem2[ibound] += ( gridboxspanSA[ievent] * 3.0 /
                                                          timespan[ievent])
                            # Include factor of 3 because timespan is in hours, not
                            # timesteps
                            averagegridboxes[ibound] += ( gridboxspan[ievent] * 3.0 /
                                                            timespan[ievent])
                            precipvolume[ibound] += totalprecipSA[ievent]
                            count[ibound] += 1
                            break


# In[ ]:
if extra == 1:
    averageydist = averageydist / count
    averagexdist = averagexdist / count

averagetime = averagetime / count
averageprecipperhr = averageprecipperhr / count
averageprecipperareahr = averageprecipperareahr / count
averagem2 = averagem2/count
averagegridboxes = averagegridboxes/count

# convert from m2 to km2
averagekm2 = averagem2 / (1000.0 * 1000.0)

# In[ ]:

# Write out to a text file

filename = (filenameadd + 'Averages_' + '{:d}'.format(minlat) + 'N-' +
                    '{:d}'.format(maxlat) + 'N_min' + str(minGB) + 'GB.txt')

with open(DirI + filename, 'w') as text_file:
    text_file.write('Domain averages for ' + '{:d}'.format(minlat) + 'N-' +
            '{:d}'.format(maxlat) + 'N and ' + '{:d}'.format(minlon) + 'E-'
            + '{:d}'.format(maxlon) + 'E \n')
    text_file.write('Created by CalculateDomainTotalAverages.py on ' +
                            str(dt.date.today()) + '\n')
    text_file.write('Input Arguments: \n') 
    text_file.write(str(args) + '\n')

    if extra == 1:
        text_file.write('timespan (hours), \t count (events/yr), \t average '
                    'latitude distance (degs), \t average longitude distance ' 
                    '(degrees) \t averagepreciphr (mm/hr), \t averagepreciphr ' 
                    '(mm/gridbox/hr) \t total precip (m3 /yr) \n')
        for ibound in range(0,nbounds):
            text_file.write('{:.1f}'.format(tbound1[ibound]) + '-' +
                        '{:.1f}'.format(tbound2[ibound]) + 'hours,   ' +
                        '{:.2e}'.format(count[ibound]/nyears) +
                        '; ' +
                        '{:.2f}'.format(averageydist[ibound]) +
                        '; ' +
                        '{:.2f}'.format(averagexdist[ibound]) +
                        '; ' +
                        '{:.2f}'.format(averagetime[ibound]) +
                        ';  ' +
                        '{:.2e}'.format(averagekm2[ibound]) +
                        '; ' +
                        '{:.2e}'.format(averageprecipperhr[ibound]) +
                        ';  ' +
                        '{:.2e}'.format(averageprecipperareahr[ibound])
                        + ';  ' +
                        '{:.2e}'.format(precipvolume[ibound]/nyears)
                        + ' \n')
    else:
        text_file.write('timespan (hours), \t count (events/yr), \t average '
                    'time, \t average footprint (km2) '
                    '\t average footprint (gridboxes) '
                    '\t averagepreciphr (gridboxes mm/hr), \taveragepreciphr '
                    '(mm/hr) \t total precip (m3 /yr) \n')
        for ibound in range(0,nbounds):
            text_file.write('{:.1f}'.format(tbound1[ibound]) + '-' +
                        '{:.1f}'.format(tbound2[ibound]) + 'hours; \t' +
                        '{:.2e}'.format(count[ibound]/nyears) +
                        '; \t' +
                        '{:.2f}'.format(averagetime[ibound]) +
                        '; \t' +
                        '{:.2e}'.format(averagekm2[ibound]) +
                        '; \t' +
                        '{:.2e}'.format(averagegridboxes[ibound]) +
                        '; \t' +
                        '{:.2e}'.format(averageprecipperhr[ibound]) +
                        '; \t' +
                        '{:.2e}'.format(averageprecipperareahr[ibound])
                        + '; \t ' +
                        '{:.2e}'.format(precipvolume[ibound]/nyears)
                        + ' \n')



# In[ ]:

datain.close()


