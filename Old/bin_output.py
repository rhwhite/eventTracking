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
import re
import sys

import xray

histnbins = 40 # number of bins

startyr = 1998
endyr = 2014 

plotdensity = False

mintsteps = 23360	# 
maxtsteps = 46720 # number of timesteps (3hrly) to include; 2920 is one year # 43800 is 15 years, 46720 is 16 years

minevent = 100000

DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Precip/'
FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + 'tspan.nc'
FileI2 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '.nc'

DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Precip/'
FileO = 'hist_precip_tspan_sspan_bins_' + str(mintsteps) + '-' + str(maxtsteps) + '_' + str(plotdensity) + '.nc'


filetimespan = "3hrly"


datain = xray.open_dataset(DirI + FileI1)
datain2 = xray.open_dataset(DirI + FileI2)

print(datain.coords)

tspan=datain['timespan']
tmean = datain['tmean']
tstart = datain2['tstart'].values
sspan = datain2['span']
totalprecip = datain2['totalprecip']
xmin1 = datain['xmin']
xmin2 = datain2['xmin']

if (xmin1 != xmin2).any():
	sys.exit('xmin from two files not equal')

nevents = tspan.shape[0]

curtime = tmean[mintsteps]
mintimeidx = 0

while curtime < mintsteps:
	mintimeidx = mintimeidx + 1
	curtime = tstart[mintimeidx]


print mintimeidx
maxtimeidx = mintimeidx
curtime = tmean[maxtimeidx]
print curtime

while curtime < maxtsteps:
	maxtimeidx = maxtimeidx + 1
	curtime = tstart[maxtimeidx]

print nevents
print 'maximum event is: ', maxtimeidx

#Bin data 

try:
    os.remove(DirO + FileO)
except OSError:
    pass

outfile = Dataset(DirO + FileO, 'w')
outfile.createDimension('bins', histnbins)

tspanO = outfile.createVariable('tspanhist','f4',('bins'),fill_value=-9999)
sspanO = outfile.createVariable('sspanhist','f4',('bins'),fill_value=-9999)
precipO = outfile.createVariable('preciphist','f4',('bins'),fill_value=-9999)
tspan_binsO = outfile.createVariable('tspanhist_bins','f4',('bins'),fill_value=-9999)
sspan_binsO = outfile.createVariable('sspanhist_bins','f4',('bins'),fill_value=-9999)
precip_binsO = outfile.createVariable('preciphist_bins','f4',('bins'),fill_value=-9999)


tspantemp,tspan_bin = np.histogram(tspan[mintimeidx:maxtimeidx],bins=histnbins,density=plotdensity)
sspantemp,sspan_bin = np.histogram(sspan[mintimeidx:maxtimeidx],bins=histnbins,density=plotdensity)
preciptemp,precip_bin = np.histogram(totalprecip[mintimeidx:maxtimeidx],bins=histnbins,density=plotdensity)

print tspantemp
tspanO[:] = tspantemp
sspanO[:] = sspantemp
precipO[:] = preciptemp

for ibin in range(0,histnbins):
        tspan_binsO[ibin] = 0.5 *(tspan_bin[ibin] + tspan_bin[ibin+1])
        sspan_binsO[ibin] = 0.5 *(sspan_bin[ibin] + sspan_bin[ibin+1])
        precip_binsO[ibin] = 0.5 *(precip_bin[ibin] + precip_bin[ibin+1])

datain.close()
datain2.close()
outfile.close()
