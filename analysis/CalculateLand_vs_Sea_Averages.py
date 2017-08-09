"""
Code to sum characteristics of events split by either time or speed and write
out the results to a table.

Author: Rachel H White rachel.white@cantab.net
Created: Oct 2016

Example use:
python CalculateLand_vs_SeaAverages.py --Data TRMM --Version Standard \
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

parser = argparse.ArgumentParser(description="Calculate land/sea averages")
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


args = parser.parse_args()

print args

# put inputs into the type and variable names we want
splittype = args.splittype[0]
speedtspan = args.speedtspan
# multiply tbound in by 24 to get hours rather than days
tbound1 = np.multiply(args.tbound1,24.0)
tbound2 = np.multiply(args.tbound2,24.0)
unit = args.unit[0]
Data = args.Data[0]
Version = args.Version[0]
startyr = args.anstartyr[0]
endyr = args.anendyr[0]
minlon = args.minlon
maxlon = args.maxlon
minlat = args.minlat
maxlat = args.maxlat
test = args.test

diradd = getdirectory(splittype)
nbounds = len(tbound1)

print(tbound1)

R = 6371000     # radius of Earth in m

nyears = endyr - startyr + 1

minevent = 100000

DirI = ('/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' +
        Version + str(startyr) + '/proc/')


if Data == "TRMM":
    if Version == '6th_from6' or Version == '5th_from48':
        DirI = ('/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' +
                Version + '/proc/')
    FileInLats = ('/home/disk/eos4/rachel/Obs/TRMM/'
                    'SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc')
elif Data == "TRMMERAIgd":
    FileInLats = ('/home/disk/eos4/rachel/Obs/TRMM/'
                    'regrid2ERAI_TRMM_3B42_1998-2014.nc')

elif Data == "ERAI":
    FileInLats = ('/home/disk/eos4/rachel/Obs/ERAI/3hrly/Precip_3hrly/'
                'SeasAnn_ERAI_Totalprecip_' +
                str(startyr) + '-' + str(endyr) + '_preprocess.nc')

elif Data == "ERA20C":
    FileInLats = '/home/disk/eos4/rachel/Obs/ERA_20C/ERA_20C_LatLon.nc'

elif Data == "CESM":
    DirI = ('/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' +
            Version + str(startyr) + '/proc/')
    FileInLats = ('/home/disk/eos4/rachel/EventTracking/Inputs/CESM/'
                    'f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
                    'f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc')
else:
    print("unexpected data type")
    exit()

DirO = DirI + diradd + '/'


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
averageydist = np.zeros([nbounds,2],float)
averagexdist = np.zeros([nbounds,2],float)
averageprecipperhr = np.zeros([nbounds,2],float)
averageprecipperareahr = np.zeros([nbounds,2],float)
averagetime = np.zeros([nbounds,2],float)
averagegridboxes = np.zeros([nbounds,2],float)

precipvolume = np.zeros([nbounds,2],float)

count = np.zeros([nbounds,2],int)


# In[5]:

# open main dataset and read in data
FileI1 = ('All_Precip_' + str(startyr) + '-' + str(endyr) + '_' + Data + '_' +
            Version + '.nc')

datain = xrayOpen(DirI + FileI1,decodetimes=False)

nevents = len(datain.events)

timespan = datain.timespan[0:nevents].values
ycenterstart = datain.ycenterstart[0:nevents].values
ycenterend = datain.ycenterend[0:nevents].values

ycentermean = datain.ycentermean[0:nevents].values
xcentermean = datain.xcentermean[0:nevents].values

timespan = datain.timespan[0:nevents].values
totalprecip = datain.totalprecip[0:nevents].values
totalprecipSA = datain.totalprecipSA[0:nevents].values
gridboxspan = datain.gridboxspan[0:nevents].values

startlats = lats[datain.ycenterstart[0:nevents].astype(int)]
endlats = lats[datain.ycenterend[0:nevents].astype(int)]


# In[6]:

# Set fileminlat and filemaxlat if we ran FiT on a subsection of the data
fileminlat = -90
filemaxlat = 90

# Need to get a land-sea mask at correct resolution
LandSeaMask = '/home/disk/eos4/rachel/Obs/TRMM/TMPA_land_sea_mask.nc'

LandSeaFile = xray.open_dataset(LandSeaMask)
LandSea = LandSeaFile['landseamask'].sel(lat=slice(fileminlat,filemaxlat))


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
    # check if in region
    if isinregion(lats[ycenterstart[ievent]],
                  lons[xcenterstart[ievent]]):

        if isinregion(lats[ycenterend[ievent]],
                      lons[xcenterend[ievent]]):

            if LandSea[ycentermean[ievent].astype(int),xcentermean[ievent].astype(int)] > 50:
                LSindex = 0
            else:
                LSindex = 1
            for ibound in range(0,nbounds):
                if timespan[ievent] < tbound2[ibound]:
                    averageydist[ibound,LSindex] += (lats[ycenterend[ievent]] -
                                            lats[ycenterstart[ievent]])
                    averagexdist[ibound,LSindex] += (lons[xcenterend[ievent]] -
                                            lons[xcenterstart[ievent]])
                    # if negative then that's fine for NH, positive is fine
                    # for southern hemisphere, so get the average event
                    # distance travelled

                    averagetime[ibound,LSindex] += timespan[ievent]
                    averageprecipperhr[ibound,LSindex] += (
                                                    totalprecip[ievent]/timespan[ievent])
                    # Include factor of 3 to convert to hours, not timesteps
                    averageprecipperareahr[ibound,LSindex] += (
                                                        totalprecip[ievent]/(3.0 *
                                                        gridboxspan[ievent]))
                    averagegridboxes[ibound,LSindex] += gridboxspan[ievent]
                    precipvolume[ibound,LSindex] += totalprecipSA[ievent]
                    count[ibound,LSindex] += 1
                    break


# In[ ]:

averageydist = averageydist / count
averagexdist = averagexdist / count

averagetime = averagetime / count
averageprecipperhr = averageprecipperhr / count
averageprecipperareahr = averageprecipperareahr / count
averagegridboxes = averagegridboxes/count


# In[ ]:

# Write out to a text file

for LSindex in range(0,2):
    if LSindex == 0:
        filename = (filenameadd + 'Averages_Seas_' + '{:d}'.format(minlat) + 'N-' +
                    '{:d}'.format(maxlat) + 'N.txt')
    elif LSindex == 1:
        filename = (filenameadd + 'Averages_Land_' + '{:d}'.format(minlat) + 'N-' +
                    '{:d}'.format(maxlat) + 'N.txt')
    with open(DirI + filename, 'w') as text_file:
        text_file.write('Domain averages for ' + '{:d}'.format(minlat) + 'N-' +
                '{:d}'.format(maxlat) + 'N and ' + '{:d}'.format(minlon) + 'E-'
                + '{:d}'.format(maxlon) + 'E \n')
        text_file.write('timespan (hours), \t count (events/yr), \t average '
                        'latitude distance (degs), \t average longitude distance ' 
                        '(degrees) \t averagepreciphr (mm/hr), \taveragepreciphr ' 
                        '(mm/gridbox/hr) \t total precip (m3 /yr) \n')
        for ibound in range(0,nbounds):
            text_file.write('{:.1f}'.format(tbound1[ibound]) + '-' +
                            '{:.1f}'.format(tbound2[ibound]) + 'hours,   ' +
                            '{:.2e}'.format(count[ibound,LSindex]/nyears) +
                            '; ' +
                            '{:.2f}'.format(averageydist[ibound,LSindex]) +
                            '; ' +
                            '{:.2f}'.format(averagexdist[ibound,LSindex]) +
                            '; ' +
                            '{:.2f}'.format(averagetime[ibound,LSindex]) +
                            ';  ' +
                            '{:.2e}'.format(averagegridboxes[ibound,LSindex]) +
                            '; ' +
                            '{:.2e}'.format(averageprecipperhr[ibound,LSindex]) +
                            ';  ' +
                            '{:.2e}'.format(averageprecipperareahr[ibound,LSindex])
                            + '; /t ' +
                            '{:.2e}'.format(precipvolume[ibound,LSindex]/nyears)
                            + ' \n')




# In[ ]:

datain.close()


