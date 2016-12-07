#-*- coding: utf-8 -*-
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

Data = "CESM"
splittype = "days"      # "speed","maxspeed"
speedtspan = 0
addTAnn = 0
# In m/s
if splittype == "speed" or splittype == "maxspeed":
        tbound = np.array([-1000,-30,-10,-6,-3,3,6,10,30]) 
        tbound2 = np.array([-30,-10,-6,-3,3,6,10,30,1000])
        unit = "ms"
elif splittype == "days":
        tbound = np.array([0,1,2,5,1])
        tbound2 = np.array([1,2,5,100,5])
        unit = "day"
nbounds = len(tbound)

mapping = "centre"

mapprecip = 1
Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

print Version

sumlats = 8 #10 8 is approximately 2 x 2 degrees for 0.25 degree data
sumlons = 8 #20

nseas = 4

if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2015
	nyears = endyr-startyr # No extra + 1 as actually finishes in 2013
	if Version == 'Standard':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'

	elif Version == '5th_nanto25':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'

	elif Version == '5th_nantozero':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'

	elif Version == '7thresh':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
	elif Version == '6th_from6':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
	elif Version == '5th_from48':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'
	else:
        	sys.exit('unexpected Version')
elif Data == "ERAI":
	startyr = 1980
	endyr = 2014
	nyears = endyr-startyr 	# no extra + 1 as actually finishes in 2013
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	DirO = DirI
elif Data == "CESM":
	Dstartyr = 1990 # Don't change - tied to file names!
	Dendyr = 2014
	startyr = 1990
	endyr = 2014
	nyears = endyr - startyr # no extra + 1 as actually finished in 2013
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Dstartyr) + '-' + str(Dendyr) + '/Precip/'
        DirO = DirI

print nyears
sumvars = ["TDensity","EDensity","WDensity","SDensity"]

if mapprecip == 0:
    add = "noP"
else:
    add = ""

for iloop in range(0,nbounds):
        if speedtspan == 0:
                fileadd = ""
        else:
                fileadd = str(speedtspan) + "ts_"

        tboundtitle = str(int(tbound[iloop])) + '-' + str(int(tbound2[iloop]))
        if splittype == "maxspeed":
                fileadd = fileadd + "MaxSpeeds_" + add
        elif splittype == "speed":
                fileadd = fileadd + "" + add
        elif splittype == "days":
                fileadd = fileadd + "Sizes_" + add

	FileI = 'DenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
	FileO = 'DenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version +'_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'


	filetimespan = "3hrly"

	print DirI + FileI
	file = xray.open_dataset(DirI + FileI)
	lons = file['Longitude'].values
	lats = file['Latitude'].values
	years = file['years'].values
	Seasons = file['Seasons'].values

	print Seasons

	nlons = len(lons)
	nlats = len(lats)
	print "nlons: ",nlons
	print "nlats: ",nlats

	print nlats, nlons
	nlonsnew = math.floor(nlons/sumlons)
	nlatsnew = math.floor(nlats/sumlats)
	print nlonsnew,nlatsnew

	nlats2 = np.int(nlatsnew * sumlats)
	nlons2 = np.int(nlonsnew * sumlons)

	print nlons2
	print nlats2

	Latsnew = np.zeros(nlatsnew,np.float)
	Lonsnew = np.zeros(nlonsnew,np.float)

	# Create Outfile structure
	ncfile = Dataset(DirO + FileO, 'w')
	ncfile.createDimension('years', nyears) 
	ncfile.createDimension('seas', nseas)
	ncfile.createDimension('size', 4)
	ncfile.createDimension('lon', nlonsnew)
	ncfile.createDimension('lat', nlatsnew)
	OLongitude = ncfile.createVariable('lon','f4',('lon'),fill_value=-9999)
	OLatitude = ncfile.createVariable('lat','f4',('lat'),fill_value=-9999)
	OYears = ncfile.createVariable('years','f4',('years'),fill_value=-9999)
	OSeas = ncfile.createVariable('seas',str,('seas'))

	OYears[:] = years
	OSeas[:] = Seasons[:]


        if addTAnn == 1:
		variable = "TDensity"
		print variable
                VarMap = file["WDensity"] + file["EDensity"] + file["SDensity"]
                # Catch all silly values!
                VarMap[np.where(VarMap > 10E30)] = np.nan
                VarMap[np.where(VarMap == -9999)] = np.nan
                VarMap[np.where(VarMap < -10E30)] = np.nan


                VarMapnew = np.zeros((nyears,nseas,nlatsnew,nlonsnew),np.float64)
                VarMapnew[...] = np.nan

		inlat=0
                for ilats in range(0,nlats2,sumlats):
                        inlon = 0
                        Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlats])
                #       print np.mean(lats[ilats:ilats+sumlats])
                        for ilons in range(0,nlons2,sumlons):
                                Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlons])
                                if variable in sumvars:
                                        VarMapnew[:,:,inlat,inlon] = np.nansum(VarMap[:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(2,3),dtype=np.float64)
                                else:
                                        VarMapnew[:,:,inlat,inlon] = np.nanmean(VarMap[:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(2,3),dtype=np.float64)

                                inlon += 1
                        inlat += 1
                # Create seasonal and or means!

                VarMapAnn = np.nansum(VarMapnew,axis=1)

                OVarMap = ncfile.createVariable("TDensityAnn",'f4',('years','lat','lon'),fill_value=-9999)
                OVarMapSeas = ncfile.createVariable("TDensitySeas",'f4',('years','seas','lat','lon'),fill_value=-9999)

                OVarMap[:,:,:] = VarMapAnn
                OVarMapSeas[:,:,:,:] = VarMapnew

	for variable in file.variables.keys(): 
		#["ETotalPrecip","WTotalPrecip","EDensity","WDensity","ESpeed","WSpeed","EDistance","WDistance","ESize","WSize","EUniSize","WUniSize","ETSpan","WTSpan","EPerc"]:
		print variable
		if variable in ['years','Seasons','lat','lon','season','Longitude','Latitude']:
			print 'not doing anything for dimensions variables'
		else:
			VarMap = file[variable].values

			# Catch all silly values!
			VarMap[np.where(VarMap > 10E30)] = np.nan
			VarMap[np.where(VarMap == -9999)] = np.nan
			VarMap[np.where(VarMap < -10E30)] = np.nan
			

			VarMapnew = np.zeros((nyears,nseas,nlatsnew,nlonsnew),np.float64)
			VarMapnew[...] = np.nan

			#'regrid' by summing over large boxes (averaging lons and lats)
			inlat = 0
			for ilats in range(0,nlats2,sumlats):
				inlon = 0
				Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlats])
			#	print np.mean(lats[ilats:ilats+sumlats])
				for ilons in range(0,nlons2,sumlons):
					Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlons])
					if variable in sumvars:
						VarMapnew[:,:,inlat,inlon] = np.nansum(VarMap[:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(2,3),dtype=np.float64)	
					else:
                                                VarMapnew[:,:,inlat,inlon] = np.nanmean(VarMap[:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(2,3),dtype=np.float64)         
	
					inlon += 1
				inlat += 1
			# Create seasonal and or means!

			VarMapAnn = np.nansum(VarMapnew,axis=1)
			
			OVarMap = ncfile.createVariable(variable + "Ann",'f4',('years','lat','lon'),fill_value=-9999)
			OVarMapSeas = ncfile.createVariable(variable + "Seas",'f4',('years','seas','lat','lon'),fill_value=-9999)
			
			if Latsnew[0] > Latsnew[-1]:
				OVarMap[:,:,:] = VarMapAnn[:,::-1,:]
				OVarMapSeas[:,:,:,:] = VarMapnew[:,::-1,:]
			else:
				OVarMap[:,:,:] = VarMapAnn
				OVarMapSeas[:,:,:,:] = VarMapnew

	# After last variable, write lats and lons, which will be identical for all variables and just need to be written once
	if Latsnew[0] > Latsnew[-1]:
		OLatitude[:] = Latsnew[::-1]
	else:
		OLatitude[:] = Latsnew
	OLongitude[:] = Lonsnew

	ncfile.close()

Ngl.end()



