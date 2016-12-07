# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os, errno
import numpy as np
import numpy.ma as ma
import netCDF4
import datetime
from netCDF4 import Dataset
import datetime as dt
import re
import sys
import Ngl
import xray

histnbins = 10 # number of bins

startyr = 1998
endyr = 2014 

nyears = endyr - startyr
nmonths = (nyears + 1) * 12

minevent = 100000

DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '5th_nanto25.nc'

FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
#figtitle1 = 'logs_precip_tspan_sspan_' + str(histnbins) + 'bins_diffs_' + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)  + '_' + str(plotdensity) + '.nc' 
#FigDir + figtitle1

DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
FileO = 'Monthly_Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '5th_nanto25.nc'

filetimespan = "3hrly"

datain = xray.open_dataset(DirI + FileI1)
datain2 = xray.open_dataset(DirI + FileI2)

print(datain.coords)

tspan=datain['timespan'].values
tstart = datain2['tstart'].values
sspan = datain2['uniquegridboxspan'].values
totalprecip = datain2['totalprecip'].values
xcentermean = datain['xcentermean'].values
ycentermean = datain['ycentermean'].values
#xmin1 = datain['xmin']
#xmin2 = datain2['xmin']

#if (xmin1 != xmin2).any():
#	sys.exit('xmin from two files not equal')

nevents = tspan.shape[0]
print nevents
print int(nevents)
attrs = {'units':'hours since 1998-01-01'}

print tstart
ds = xray.Dataset({'time': ('time',tstart,attrs)})

dates = xray.decode_cf(ds)

years = dates['time.year']
months = dates['time.month']


try:
    os.remove(DirO + FileO)
except OSError:
    pass

dimsize1 = nevents/10
dimsize2 = nmonths

print dimsize2, dimsize1

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('events', dimsize1)
ncfile.createDimension('months', dimsize2)

Msspan = ncfile.createVariable('uniquegridboxspan','f4',('months','events'),fill_value=-9999)
Mprecip = ncfile.createVariable('totalprecip','f4',('months','events'),fill_value=-9999)
Mtimespan = ncfile.createVariable('timespan','f4',('months','events'),fill_value=-9999)
#Mtimestart = ncfile.createVariable('tstart','f4',('months','events'),fill_value=-9999)
#Mtimemean = ncfile.createVariable('tmean','f4',('months','events'),fill_value=-9999)
Mxcentmean = ncfile.createVariable('ycentermean','f4',('months','events'),fill_value=-9999)
Mycentmean = ncfile.createVariable('xcentermean','f4',('months','events'),fill_value=-9999)
Mdatamonth = ncfile.createVariable('month','f4',('months'),fill_value=-9999)
#Mdatayear = ncfile.createVariable('year','f4',('months'),fill_value=-9999)
Meventsnum = ncfile.createVariable('eventsnum','f4',('months'),fill_value=-9999)

idxmonth = 0
idxevent = 0
curmonth = 1

Mdatamonth[idxmonth] = months[0]
print months[0]
print months[nevents-1]

exit

for ievent in range(0,nevents):
	if months[ievent] != curmonth:
		print 'new month', months[ievent]
		curmonth = months[ievent] 
		Meventsnum[idxmonth] = idxevent + 1
		idxmonth = idxmonth + 1
		idxevent = 0
		print months[ievent]
		print years[ievent]
                print "ievent: " + str(ievent)

	Msspan[idxmonth,idxevent] = sspan[ievent]		
        Mprecip[idxmonth,idxevent] = totalprecip[ievent]
	Mtimespan[idxmonth,idxevent] = tspan[ievent]	
#        Mtimemean[idxmonth,idxevent] = tmean[ievent]
#        Mtimestart[idxmonth,idxevent] = tstart[ievent]
	Mxcentmean[idxmonth,idxevent] = xcentermean[ievent]
        Mycentmean[idxmonth,idxevent] = ycentermean[ievent]

	Mdatamonth[idxmonth] = months[ievent]
#        Mdatayear[idxmonth] = years[ievent]

	idxevent = idxevent + 1

ncfile.close()





