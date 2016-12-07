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
from rhwhitepackages.calc_seas_ann import calcann
from rhwhitepackages.calc_seas_ann import calcseas


Data = "TRMM"
Version = "ERAIgd"

print Data

startyr = 1998 
endyr = 2014

if Data == "TRMM":
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
elif Data == "OLR":
	DirIn = "/home/disk/eos4/rachel/Obs/OLR/"
	FileIn = "olr.mon.mean.nc"

print DirIn + FileIn
OutputfileA = DirIn + "Ann_" + FileIn
OutputfileS = DirIn + "Seas_" + FileIn
print OutputfileA

#Get lons and lats
FileIn = xray.open_dataset(DirIn + FileIn)

if Data == "CMAP":
	lats = FileIn['lat']
	lons = FileIn['lon']
	vardata = FileIn['precip']
elif Data == "OLR":
	lats = FileIn['lat']
	lons = FileIn['lon']
	vardata = FileIn['olr']
elif Data == "ERAI":
	lats = FileIn['latitude']
	lons = FileIn['longitude']
	vardata = FileIn['tpnew']
elif Data == "TRMM":
        lats = FileIn['latitude']
        lons = FileIn['longitude']
	vardata = FileIn['pcp']
else:
	print "unknown data type ", Data
	exit()

nlats = len(lats)
nlons = len(lons)

dates = FileIn['time']

years = dates['time.year']
seasons = dates['time.season']
months = dates['time.month']
nyears = np.max(years) - np.amin(years)

yearoutput = range(np.min(years),np.max(years))
ntimes = len(years)

startyear = int(years[0])
endyear = int(years[ntimes-1])

# Create output file
#Write out to new file
ncfileA = Dataset(OutputfileA,'w')
ncfileA.createDimension('lon', nlons)
ncfileA.createDimension('lat', nlats)
ncfileA.createDimension('year',None)
ncfileA.createDimension('season',4)

ncfileS = Dataset(OutputfileS,'w')
ncfileS.createDimension('lon', nlons)
ncfileS.createDimension('lat', nlats)
ncfileS.createDimension('year',None)
ncfileS.createDimension('season',4)


if Data == "OLR":
        DataSeasS = ncfileS.createVariable('OLRSeas','f4',('year','season','lat','lon'),fill_value=-9999)
        DataClimSeasS = ncfileS.createVariable('OLRClimSeas','f4',('season','lat','lon'),fill_value=-9999)

        DataAnnA = ncfileA.createVariable('OLRAnn','f4',('year','lat','lon'),fill_value=-9999)
        DataClimAnnA = ncfileA.createVariable('OLRClimAnn','f4',('lat','lon'),fill_value=-9999)
else:
        DataSeasS = ncfileS.createVariable('PrecipSeas','f4',('year','season','lat','lon'),fill_value=-9999)
        DataClimSeasS = ncfileS.createVariable('PrecipClimSeas','f4',('season','lat','lon'),fill_value=-9999)

        DataAnnA = ncfileA.createVariable('PrecipAnn','f4',('year','lat','lon'),fill_value=-9999)
	DataClimAnnA = ncfileA.createVariable('PrecipClimAnn','f4',('lat','lon'),fill_value=-9999)

DataYearsA = ncfileA.createVariable('Year','i4',('year'))
DataLatsA = ncfileA.createVariable('Latitude','f4',('lat'))
DataLonsA = ncfileA.createVariable('Longitude','f4',('lon'))

DataYearsS = ncfileS.createVariable('Year','i4',('year'))
DataLatsS = ncfileS.createVariable('Latitude','f4',('lat'))
DataLonsS = ncfileS.createVariable('Longitude','f4',('lon'))
DataSeasonsS = ncfileS.createVariable('Season','S3',('season'))

DataLatsA[:] = lats.values
DataLonsA[:] = lons.values

DataSeasonsS[:] = Seas[:]
DataLatsS[:] = lats.values
DataLonsS[:] = lons.values


# Calculate values

if ntimes > 10000:
	# Doing one season and one year at a time for lots of data
#
	print 'need to know in advance how big the arrays are -> need to know whether first and last years count!'
	seasdata = np.zeros([nyears,4,nlats,nlons],np.float)
	anndata = np.zeros([nyears,nlats,nlons],np.float)

	for iyear in range(0,nyears):
		print iyear
		selyear = startyear + iyear
		selmonths = months.sel(time=str(selyear))
		print selmonths[-1]
		if selmonths[-1] != 12:
			print "not a full year"
			break
		npoints = vardata.sel(time=str(selyear)).shape[0]
		print npoints
		print years.sel(time=str(selyear))
		print vardata.sel(time=str(selyear)).shape
		anndata[iyear,...] = calcann(0,npoints,vardata.sel(time=str(selyear)),years.sel(time= str(selyear)))
	
        DataAnn[:,:,:] = anndata
	DataClimAnn[:,:] = np.nanmean(anndata[0:nyears,:,:],axis = 0)
	
	seasdata[...] = np.nan
	for iyear in range(0,nyears-1):
		print iyear
		for iseas in range(0,4):
			selyear = startyear + iyear
			selyear2 = selyear
			Sseason = seasonstarts[iseas]
			Eseason = Sseason + 2
			if Eseason > 12:
				Eseason = Eseason-12
				selyear2 = selyear + 1

			npoints = vardata.sel(time=slice(str(selyear)+'-'+str(Sseason),str(selyear2)+'-' + str(Eseason))).shape[0]

			seasdata[iyear,iseas,...] = calcseas(0,npoints,vardata.sel(time=slice(str(selyear)+'-'+ str(Sseason),str(selyear2)+'-'+str(Eseason))),years.sel(time=slice(str(selyear)+'-'+str(Sseason),str(selyear2)+'-'+str(Eseason))),months.sel(time = slice(str(selyear)+'-'+str(Sseason),str(selyear2)+'-'+str(Eseason))).values, seasonstarts,seasmonths)
	
        DataSeas[:,:,:,:] = seasdata
	DataClimSeas[:,:,:] = np.nanmean(seasdata[0:nyears-1,:,:,:],axis = 0)

        DataYearsA[:] = yearoutputA
        DataYearsS[:] = yearoutputS


else:
	# Do whole lot at one time:
	anndata,yearoutputA = calcann(0,ntimes,vardata.values,years.values)
	

	seasdata,yearoutputS = calcseas(0,ntimes,vardata,years.values,months.values,seasonstarts,seasmonths)

        DataAnnA[:,:,:] = anndata
        DataClimAnnA[:,:] = np.nanmean(anndata[0:nyears,:,:],axis = 0)
	DataYearsA[:] = yearoutputA

        DataSeasS[:,:,:,:] = seasdata
        DataClimSeasS[:,:,:] = np.nanmean(seasdata[0:nyears-1,:,:,:],axis = 0)
        DataYearsS[:] = yearoutputS

ncfileA.close()

ncfileS.close()

Ngl.end()




