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

uvindir = '/home/disk/eos4/rachel/Obs/ERAI/uv'
uvfile = xr.open_mfdataset(uvindir + '/*.grb')

ulev, vlev = uvfile['u'], uvfile['v']
EKE = 0.5 * ((ulev * ulev) + (vlev * vlev))

EKEda = xr.DataArray(EKE,coords=[('time',ulev.time),('lat',ulev.latitude),('lon',ulev.longitude),('lev',ulev.level)])

EKEds = xr.Dataset({variable:EKE})

#EKEds.to_netcdf(uvindir + '/EKE_1998-2015.nc',mode='w',format='NETCDF4')

