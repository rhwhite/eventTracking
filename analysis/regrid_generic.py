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
import Ngl
import xray
import math
from rhwhitepackages.readwrite import XrayOpen

#Data = "CESM"
#DirIn = "/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/"
#FileIn = "ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc"
#latin = 'latitude'
#lonin = 'longitude'
#datain = ['
#DirIn = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
#FileIn = 'TMPA_land_sea_mask.nc'

#DirIn = '/home/disk/eos4/rachel/Obs/ERA_20C/'
#FileIn = 'ERA_20C_Ann_Totalprecip_1980-2011.nc'

Data = "TRMM"
DirIn = '/home/disk/eos4/rachel/Obs/TRMM/'
FileIn = 'TRMM_3B42_1998-2014_annmean.nc'
latin = 'latitude'
lonin = 'longitude'

#Data = "ERAI"
#DirIn = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
#FileIn = 'ncra_ERAI_Totalprecip_1980-2015_preprocess.nc'

#latin = 'latitude'
#lonin = 'longitude'
#datain = ['tpnew']

sumlats = 16
sumlons = 16

filetimespan = '3hrly'

DirO = DirIn 
FileO = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_' + FileIn
#pcp, time, ongitude, latitude
precipdata = XrayOpen(DirIn + FileIn)

lons = precipdata[lonin]
lats = precipdata[latin]

nlons = len(lons)
nlats = len(lats)

nlonsnew = nlons/sumlons
nlatsnew = nlats/sumlats

Latsnew = np.zeros(nlatsnew,np.float)
Lonsnew = np.zeros(nlonsnew,np.float)

nlonsnew = math.floor(nlons/sumlons)
nlatsnew = math.floor(nlats/sumlats)

nlats2 = np.int(nlatsnew * sumlats)
nlons2 = np.int(nlonsnew * sumlons)


for variable in precipdata.variables.keys():
	print variable
	indata = precipdata[variable]

	if variable in ['time','years','lat','lon','season','Longitude','Latitude','latitude','longitude','bnds','time_bnds']:
		print 'not doing anything for dimensions variables'
		continue
	try:
		units = indata.units
		print units
		print indata.longname
		if units == "mm/hr":
			pass
		elif units == "m/s":
			indata = indata * 1000.0 * 60.0 * 60.0	# convert from m/s to mm/hr
		else:
			if indata.longname == "precipitation (mm/hr)":
				units = "mm/hr"
			else:
				exit("unexpected unit: " + indata.unit)
	except AttributeError:
		if indata.long_name == "precipitation (mm/hr)":
			units = "mm/hr"
		else:
			exit("unexpected unit: " + indata.unit)

	indatacoords = indata.coords
	indatadims = indata.dims
	print indatadims
	newdimsize = []
	nnewdims = 0
	countdim=0
	newdims=[]
	for idim in indatadims:
		print idim
		if idim == latin:
			latdimnum = countdim
			newdims.append('lat')
		elif idim == lonin:
			londimnum = countdim
			newdims.append('lon')
		else:
			nnewdims = nnewdims + 1
			newdimsize.append(len(indata[idim]))
			newdims.append(idim)
		countdim += 1

	print latdimnum,londimnum
	newdimsize.append(nlatsnew)
	newdimsize.append(nlonsnew)

	New = np.zeros((newdimsize),np.float)
#'regrid' by averaging over large boxes (averaging lons and lats)

	inlat = 0
	for ilats in range(0,nlats2,sumlats):
		print ilats
		inlon = 0
		Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlats])
		for ilons in range(0,nlons2,sumlons):
			Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlons])
			New[...,inlat,inlon] = np.mean(indata[...,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(latdimnum,londimnum),dtype=np.float)

			inlon += 1
		inlat += 1
	
	foo = xray.DataArray(New,coords=[('time',indata.time),('lat',Latsnew),('lon',Lonsnew)],attrs=[('units',units)])

	newDataset = xray.Dataset({variable:foo})
	
	newDataset.to_netcdf(DirO + FileO,mode='w')

Ngl.end()


