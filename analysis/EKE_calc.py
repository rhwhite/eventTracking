# coding: utf-8
# Script to calculate EKE from ERA-interim files (currently 6-hourly)
# Taken from ipython notebook
#
# Written by Rachel White (rhwhite@uw.edu) on 15-12-2016

import os, errno
import netCDF4
import numpy as np
import datetime as dt
import pandas as pd
import xray as xr
#import Ngl
#import math
from scipy import stats
from rhwhitepackages.readwrite import shiftlons
from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.stats import regressmaps
from rhwhitepackages.readwrite import getdenfilename

# plotting
import xray.plot as xplt

uvindir = '/home/disk/eos4/rachel/Obs/ERAI/uv'
startyr = 1998
endyr = 2015
for iyear in range(startyr,endyr):
    uvfile = xr.open_mfdataset(uvindir + '/interim_daily_' + str(iyear) +
'*.grb')
    ulev, vlev = uvfile['u'], uvfile['v']
    udash = ulev - ulev.mean(dim='longitude')
    vdash = vlev - vlev.mean(dim='longitude')
    EKEall = 0.5 * ((ulev * ulev) + (vlev * vlev))
    EKEyears = EKEall.groupby('time.month').sum(dim='time')
    EKEyears = EKEyears.rename({'month':'time'})
    EKEyears = EKEyears.rename({'latitude':'lat'})
    EKEyears = EKEyears.rename({'longitude':'lon'})
    EKEds = xr.Dataset({'EKE':EKEyears})
    EKEds.to_netcdf(uvindir + '/EKE_' + str(iyear) + '.nc',mode='w')


uvindir = '/home/disk/eos4/rachel/Obs/ERAI/uv'
startyr = 1998
endyr = 2015
for iyear in range(startyr,endyr):
    for imonth in range(0,12):
        uvfile = xr.open_mfdataset(uvindir + '/interim_daily_' + str(iyear) +
'{:02d}'.format(imonth+1) + '.grb')
        ulev, vlev = uvfile['u'], uvfile['v']
        udash = ulev - ulev.mean(dim='time')
        vdash = vlev - vlev.mean(dim='time')
        EKEall = 0.5 * ((udash * udash) + (vdash * vdash))
        EKEmonth = EKEall.mean(dim='time')
        EKEmonth = EKEmonth.rename({'latitude':'lat'})
        EKEmonth = EKEmonth.rename({'longitude':'lon'})
        EKEtime = ulev.time
        s = pd.Series(EKEtime)
        smin = pd.to_datetime(s.min(),unit='ms')
        averagetime = (s-s.min()).astype('m8[ms]').mean()
        averagedate = smin + dt.timedelta(milliseconds=averagetime)
        EKEds = xr.Dataset({'EKE':EKEmonth},{'time':averagedate})

        EKEds.to_netcdf(uvindir + '/EKE_' + str(iyear) +
'{:02d}'.format(imonth+1) + '.nc',mode='w')


