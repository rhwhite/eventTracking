# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:45:25 2015

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
from scipy import stats

Data = "TRMM"
Version = "ERAIgd"

print "***** caution, this script assumes that all years are full years, and calculates mean accordingly***"

Seas = ['DJF','MAM','JJA','SON','Ann']
nseas = 4

FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2014
	DirIn = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'

        if Version == "ERAIgd":
                FileIn = 'regrid2ERAI_TRMM_3B42_1998-2014.nc'
        else:
                FileIn = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'

elif Data == "ERAI":
	DirIn = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
	FileIn = 'ERAI_Totalprecip_1980-2015_preprocess.nc'
elif Data == "CMAP":
	DirIn = "/home/disk/eos4/rachel/Obs/CMAP/"
	FileIn = "precip.mon.mean.nc"
elif Data == "ERA20C":
        startyr = 1980 # Don't change - tied to file names!
	endyr = 2011
	DirIn = '/home/disk/eos4/rachel/Obs/ERA_20C/'
        FileIn = 'ERA_20C_Totalprecip_' + str(startyr) + '-' + str(endyr) + '_preprocess.nc'

print DirIn + FileIn
Outputfile = DirIn + "SeasAnn_" + FileIn
print Outputfile

#Get lons and lats
FileIn = xray.open_dataset(DirIn + FileIn)

if Data == "CMAP":
	lats = FileIn['lat']
	lons = FileIn['lon']
	Precip = FileIn['precip']
elif Data == "ERA20C":
        lats = FileIn['latitude']
        lons = FileIn['longitude']
        Precip = FileIn['tpnew']
elif Data == "TRMM":
        lats = FileIn['latitude']
        lons = FileIn['longitude']
	Precip = FileIn['pcp']
else:
	lats = FileIn['latitude']
	lons = FileIn['longitude']
	
nlats = len(lats)
nlons = len(lons)

dates = FileIn['time']

years = dates['time.year']
seasons = dates['time.season']
months = dates['time.month']
nyears = np.max(years) - np.amin(years) + 1

yearoutput = range(np.min(years),np.max(years)+1)
print yearoutput
ntimes = len(years)
print ntimes
exit()
seasprecip = np.zeros([nyears,4,nlats,nlons],np.float)
annprecip = np.zeros([nyears,nlats,nlons],np.float)

for icount in range(0,nyears):	
	print icount
	iyear = years[icount] - years[0]
	if months[icount] == 1:
		annprecip[iyear,:,:] = (Precip[icount + 0,:,:] * 31.0 + Precip[icount+1,:,:] * 28.0 + Precip[icount+2,:,:] * 31.0 + Precip[icount+3,:,:] * 30.0 + Precip[icount+4,:,:] * 31.0 + Precip[icount+5,:,:] * 30.0 + Precip[icount+6,:,:] * 31.0 + Precip[icount+7,:,:] * 31.0 + Precip[icount+8,:,:] * 30.0 + Precip[icount+9,:,:] * 31.0 + Precip[icount+10,:,:] * 30.0 + Precip[icount+11,:,:] * 31.0) / 365.0

	"""	
	elif months[icount] == 3:
		seasprecip[iyear,1,:,:] = (Precip[icount,:,:] * 31.0 + Precip[icount + 1,:,:] * 30.0 + Precip[icount+2,:,:] * 31.0)/92.0
	elif months[icount] == 6:
                seasprecip[iyear,2,:,:] = (Precip[icount,:,:] * 30.0 + Precip[icount + 1,:,:] * 31.0 + Precip[icount+2,:,:] * 31.0)/92.0
	elif months[icount] == 9:
                seasprecip[iyear,3,:,:] = (Precip[icount,:,:] * 30.0 + Precip[icount + 1,:,:] * 31.0 + Precip[icount+2,:,:] * 30.0)/91.0
        elif months[icount] == 12:
                seasprecip[iyear,0,:,:] = (Precip[icount,:,:] * 31.0 + Precip[icount + 1,:,:] * 31.0 + Precip[icount+2,:,:] * 28.0)/90.0
	"""
ncfile = Dataset(Outputfile,'w')
ncfile.createDimension('lon', nlons)
ncfile.createDimension('lat', nlats)
ncfile.createDimension('year',nyears)
#ncfile.createDimension('season',4)

#PrecipSeas = ncfile.createVariable('PrecipSeas','f4',('year','season','lat','lon'),fill_value=-9999)
PrecipAnn = ncfile.createVariable('PrecipAnn','f4',('year','lat','lon'),fill_value=-9999)
PrecipYears = ncfile.createVariable('Year','i4',('year'))
PrecipLats = ncfile.createVariable('Latitude','f4',('lat'))
PrecipLons = ncfile.createVariable('Longitude','f4',('lon'))

#PrecipSeas[:,:,:,:] = seasprecip
PrecipAnn[:,:,:] = annprecip
PrecipYears[:] = yearoutput
PrecipLats[:] = lats.values
PrecipLons[:] = lons.values
ncfile.close()

Ngl.end()




