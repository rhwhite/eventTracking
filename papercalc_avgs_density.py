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
import math
from scipy import stats
from rhwhitepackages.readwrite import shiftlons
from rhwhitepackages.readwrite import XrayOpen

day1 = [0,1,2,5]
day2 = [1,2,5,100]

ndayplots = len(day1)

FillValue = -9999
Data = "TRMM"
mapping = 'centre'
Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

MinLonF = 0
MaxLonF = 360
MinLatF = -45
MaxLatF = 45

sumlats = 16
sumlons = 16

Seas = ['MAM','JJA','SON','DJF']
nseas = 4

# Time period for analysis
astartyr = 1998
aendyr = 2015
nyears = aendyr - astartyr + 1
print nyears
FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2015
	if Version == 'Standard':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/Old/'
		print "!!!!!!!!!!!!!!!!!!!!!!!!"
		print "YOU SHOULD BE AWARE THAT THIS IS USING OLD DATA!! REMEMBER TO UPDATE!!!"
		print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	elif Version == '5th_nanto25':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
	elif Version == '5th_nantozero':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
	elif Version == '7thresh':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
	elif Version == '6th_from6':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
	elif Version == '5th_from48':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'

	PrecipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/3hrly/"
	PrecipClimFile = "Regrid_" + str(sumlats) + "_" + str(sumlons) + "_SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc"

	FileInPrecip = XrayOpen(PrecipClimDir + PrecipClimFile)
	PrecipIn = FileInPrecip["PrecipAnnClim"].sel(lat=slice(MinLatF,MaxLatF))

	LandSeaMask = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/Regrid_' + str(sumlats) + '_' + str(sumlons) + '_TMPA_land_sea_mask.nc'

elif Data == "ERAI":
        if sumlats == 16:
		sumlats = 5
	if sumlons == 16:
		sumlons = 5
	startyr = 1980 # Don't change - tied to file names!
        endyr = 2014

        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	
	PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
	PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ncra_ERAI_Totalprecip_1980-2015_preprocess.nc' #'SeasAnn_ERAI_Totalprecip_1980-2015_preprocess.nc'

	FileInPrecip = XrayOpen(PrecipClimDir + PrecipClimFile)
	latin = FileInPrecip['lat']
	
	if latin[0] > latin[1]:
		PrecipIn = FileInPrecip['tpnew'][:,::-1,:].sel(lat=slice(MinLatF,MaxLatF)) #["PrecipClimAnn"]
	else:
                PrecipIn = FileInPrecip['tpnew'].sel(lat=slice(MinLatF,MaxLatF)) #["PrecipClimAnn"]

	print PrecipIn.shape

elif Data == "CESM":
	startyr = 1990 # Don't change - tied to file names!
	endyr = 2014
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
	
	PrecipClimDir = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
	PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'
	FileInPrecip = XrayOpen(PrecipClimDir + PrecipClimFile)
        latin = FileInPrecip['lat']

        if latin[0] > latin[1]:
                PrecipIn = FileInPrecip['PRECT'][:,::-1,:].sel(lat=slice(MinLatF,MaxLatF))
        else:
                PrecipIn = FileInPrecip['PRECT'].sel(lat=slice(MinLatF,MaxLatF))

	print PrecipIn.shape

if PrecipIn.units == "mm/hr":
	PrecipIn = PrecipIn * 24.0 #convert to mm/day
elif PrecipIn.units == "mm/day":
	pass	
else:
	error("unexpected unit in Precip file")

def regressmaps(m,c,r,p,stderr,A,linreg,nlats,nlons):
# Calculate linear trend
	for ilat in range(0,nlats):
		for ilon in range(0,nlons):
			m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(A,linreg[:,ilat,ilon])

def plotmap(plotvars1,plotvars2,plotmin1,plotmax1,plotmin2,plotmax2,vartitle1,vartitle2,title,figtitle,lons,lats,minlon,maxlon,minlat,maxlat):

	nplots = plotvars1.shape[0]
	print nplots
	wkres = Ngl.Resources()
	wkres.wkColorMap = "precip_diff_12lev"
	wks_type = "eps"
	wks = Ngl.open_wks(wks_type,figtitle,wkres)

	res = Ngl.Resources()
	res.cnInfoLabelOn         = False    # Turn off informational
						      # label.
	res.pmLabelBarDisplayMode = "Always" # Turn on label bar.
	res.cnLinesOn             = False    # Turn off contour lines.
	res.nglDraw  = False
	res.nglFrame = False

	res.sfMissingValueV = FillValue

	res.cnFillOn = True
	res.cnMissingValFillColor = "white"
	res.cnLineLabelsOn       = False
	res.pmLabelBarDisplayMode = "Always"
	res.cnLinesOn =  False

	# if lons start negative, shift everything over so there isn't a line down the middle of the Pacific
	if lons[0] < 0:
		nlonhalf = nlons/2
		lonsnew = np.zeros(lons.shape,np.float)
		lonsnew[0:nlonhalf] = lons[nlonhalf:nlons]
		lonsnew[nlonhalf:nlons] = lons[0:nlonhalf] + 360.0
		lons = lonsnew

		for iplot in range(0,nplots):
			plotvars1[iplot] = shiftlons(plotvars1[iplot],lons)
                        plotvars2[iplot] = shiftlons(plotvars2[iplot],lons)
	else:
		lonsnew = lons

	res.sfXCStartV = float(lonsnew[0])
	res.sfXCEndV = float(lonsnew[len(lons)-1])
	res.sfYCStartV = float(lats[0])
	res.sfYCEndV = float(lats[len(lats)-1])

	res.mpProjection = "CylindricalEquidistant" # Change the map projection.
	res.mpCenterLonF = 180.           # Rotate the projection.
	res.mpFillOn     = True           # Turn on map fill.

	res.lbOrientation   = "Vertical"
	res.mpLimitMode = "LatLon"    # Limit the map view.
	res.mpMinLonF = MinLonF
	res.mpMaxLonF = MaxLonF
	res.mpMinLatF = MinLatF
	res.mpMaxLatF = MaxLatF
	res.mpOutlineBoundarySets = "AllBoundaries"

	res.lbLabelFontHeightF = 0.0125
	res.lbTitleFontHeightF = 0.0125
	
	res.tiMainFontHeightF = 0.015

	res.cnLevelSelectionMode = "ManualLevels" # Define your own

	toplot = []


	for iplot in range(0,nplots):
		tempplot = plotvars1[iplot]
		tempplot[np.where(np.isnan(tempplot))] = FillValue
		res.cnMinLevelValF       = plotmin1[iplot]          # contour levels.
		res.cnMaxLevelValF       = plotmax1[iplot]
		res.cnLevelSpacingF      = ((plotmax1[iplot]-plotmin1[iplot])/10.0)
		if iplot == 0:
			res.tiMainString = vartitle1[iplot]
		else:
			res.tiMainString = vartitle1[iplot] + "; mean = {:2.3g}".format(np.nanmean(tempplot)) + "%"
		toplot.append(Ngl.contour_map(wks,tempplot,res))

                tempplot = plotvars2[iplot]
                tempplot[np.where(np.isnan(tempplot))] = FillValue
                res.cnMinLevelValF       = plotmin2[iplot]          # contour levels.
                res.cnMaxLevelValF       = plotmax2[iplot]
                res.cnLevelSpacingF      = ((plotmax2[iplot]-plotmin2[iplot])/10.0)
                if iplot == 0:
			res.tiMainString = vartitle2[iplot] 
		else:
			res.tiMainString = vartitle2[iplot] + "; mean = {:2.3g}".format(np.nanmean(tempplot)) + "%"
                toplot.append(Ngl.contour_map(wks,tempplot,res))

	
	textres = Ngl.Resources()
	textres.txFontHeightF = 0.015
	Ngl.text_ndc(wks,title,0.5,0.87,textres)


	panelres = Ngl.Resources()
	panelres.nglPanelLabelBar = True
	#panelres.nglPanelYWhiteSpacePercent = 5.
	#panelres.nglPanelXWhiteSpacePercent = 5.

	panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
	panelres.nglPanelTop                      = 0.95
	panelres.nglPanelBottom                      = 0.01

	#panelres.nglPanelFigureStrings            = ["a","b","c","d","e","f"]
	#panelres.nglPanelFigureStringsJust        = "BottomLeft"

	panelres.nglPaperOrientation = "Auto"

	plot = Ngl.panel(wks,toplot,[nplots,2],panelres)

# Get lats and lons
iday = 0
FileIn = 'DenDirSpd_Map_Sizes_' + str(day1[iday]) + '-' + str(day2[iday]) + 'day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
print DirIn + FileIn
FileIn = XrayOpen(DirIn + FileIn)

lats = FileIn['lat'].sel(lat=slice(MinLatF,MaxLatF))
lons = FileIn['lon'].values
years = FileIn['years']

nlats = len(lats)
nlons = len(lons)
nyears = len(years)

difflat = lats[1]-lats[0]
difflon = lons[1]-lons[0]

Dendays = np.zeros([ndayplots+1,nyears,nlats,nlons])	# +1 for total sum as first plot
Precipdays = np.zeros([ndayplots+1,nyears,nlats,nlons])	# + 1 for total precip as first plot

for iday in range(0,ndayplots):
        FileIn = 'DenDirSpd_Map_Sizes_' + str(day1[iday]) + '-' + str(day2[iday]) + 'day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

	#Get lons and lats
	FileIn = XrayOpen(DirIn + FileIn)
	Dendays[iday,:,:,:] = FileIn['TDensityAnn'].sel(lat=slice(MinLatF,MaxLatF)) #np.nansum([FileIn['WDensityAnn'],FileIn['EDensityAnn'],FileIn['SDensityAnn']],axis=0)
	Precipdays[iday,:,:,:] = FileIn['TPrecipAnn'].sel(lat=slice(MinLatF,MaxLatF))


Dendays[ndayplots,:,:] = np.nansum(Dendays[0:ndayplots-1,:,:],axis=0)
Precipdays[ndayplots,:,:] = np.nansum(Precipdays[0:ndayplots-1,:,:],axis=0)

#Create masks
#Create latitude array

latsarray = np.transpose(np.broadcast_to(lats,(nlons,nlats)))
minlat = np.amin(latsarray)
maxlat = np.amax(latsarray)
NHlats = np.ma.masked_greater(latsarray,0)
SHlats = np.ma.masked_less(latsarray,0)
EQlats = np.ma.masked_equal(latsarray,0)
# Need to get a land-sea mask at this resolution!! 
LandSeaFile = xray.open_dataset(LandSeaMask)
LandSea = LandSeaFile['landseamask'].sel(lat=slice(minlat,maxlat))


Lmask = np.ma.masked_less_equal(LandSea,50)
Smask = np.ma.masked_greater(LandSea,50)

Ones = np.zeros([nlats,nlons],int)
Ones[...] = 1
print "percentage of total this is land (not taking into account changes in surface area):"
print 100.0 * np.nansum(Ones * Lmask.mask)/np.nansum(Ones)
print "percentage of NH that is land: (not taking into account changes in surface area)"
print 100.0 * np.nansum(Ones * NHlats.mask * Lmask.mask)/np.nansum(Ones * NHlats.mask)
print "percentage of SH that is land:"
print 100.0 * np.nansum(Ones * SHlats.mask * Lmask.mask)/np.nansum(Ones * SHlats.mask)


DenAnnAvg = np.nanmean(Dendays,axis=1)
PrecipAnnAvg = np.nanmean(Precipdays,axis=1)

np.set_printoptions(precision=1)

noneqSum = np.nansum(PrecipAnnAvg*NHlats.mask,axis=(1,2)) + np.nansum(PrecipAnnAvg*SHlats.mask,axis=(1,2))
nhSum = np.nansum(PrecipAnnAvg*NHlats.mask,axis=(1,2))
shSum = np.nansum(PrecipAnnAvg*SHlats.mask,axis=(1,2))
totalSum = np.sum(PrecipAnnAvg,axis=(1,2)) 

noneqDSum = np.nansum(DenAnnAvg*NHlats.mask,axis=(1,2)) + np.nansum(DenAnnAvg*SHlats.mask,axis=(1,2))
nhDSum = np.nansum(DenAnnAvg*NHlats.mask,axis=(1,2))
shDSum = np.nansum(DenAnnAvg*SHlats.mask,axis=(1,2))
totalDSum = np.sum(DenAnnAvg,axis=(1,2))

# create total ignoring equator values:
print "Fractions relative to each category, ignoring equator ones"
print "Total Precip:", 100.0 * np.nansum(PrecipAnnAvg,axis=(1,2))/totalSum
print "NH Precip:   ", 100.0 * np.nansum(PrecipAnnAvg*NHlats.mask,axis=(1,2))/noneqSum
print "SH Precip:   ", 100.0 * np.nansum(PrecipAnnAvg*SHlats.mask,axis=(1,2))/noneqSum
print "Eq Precip:   ", 100.0 * np.nansum(PrecipAnnAvg*EQlats.mask,axis=(1,2))/totalSum

print "Land Precip: ", 100.0 * np.nansum(PrecipAnnAvg*Lmask.mask,axis=(1,2))/totalSum
print "Sea Precip:  ", 100.0 * np.nansum(PrecipAnnAvg*Smask.mask,axis=(1,2))/totalSum
print "and of the precip in the SH/NH:"
print "NH Land Precip: ", 100.0 * np.nansum(PrecipAnnAvg*Lmask.mask*NHlats.mask,axis=(1,2))/nhSum
print "SH Land Precip: ", 100.0 * np.nansum(PrecipAnnAvg*Lmask.mask*SHlats.mask,axis=(1,2))/shSum


np.set_printoptions(precision=1)
print "Total Density:", 100.0 * np.nansum(DenAnnAvg,axis=(1,2))/totalDSum
print "NH Density:  ", 100.0 * np.nansum(DenAnnAvg*NHlats.mask,axis=(1,2))/noneqDSum
print "SH Density:  ", 100.0 * np.nansum(DenAnnAvg*SHlats.mask,axis=(1,2))/noneqDSum
print "EQ Density:  ", 100.0 * np.nansum(DenAnnAvg*EQlats.mask,axis=(1,2))/totalDSum

print "Land Density:", 100.0 * np.nansum(DenAnnAvg*Lmask.mask,axis=(1,2))/totalDSum
print "Sea Density: ", 100.0 *np.nansum(DenAnnAvg*Smask.mask,axis=(1,2))/totalDSum

print "and of the events in the SH/NH:"
print "NH Land Density: ", 100.0 * np.nansum(DenAnnAvg*Lmask.mask*NHlats.mask,axis=(1,2))/nhDSum
print "SH Land Density: ", 100.0 * np.nansum(DenAnnAvg*Lmask.mask*SHlats.mask,axis=(1,2))/shDSum


noneqSumAll = np.nansum(PrecipAnnAvg[0:4,...]*NHlats.mask,axis=(0,1,2)) + np.nansum(PrecipAnnAvg[0:4,...]*SHlats.mask,axis=(0,1,2))
nhSumAll = np.nansum(PrecipAnnAvg[0:4,...]*NHlats.mask,axis=(0,1,2))
shSumAll = np.nansum(PrecipAnnAvg[0:4,...]*SHlats.mask,axis=(0,1,2))
totalSumAll = np.sum(PrecipAnnAvg[0:4,...],axis=(0,1,2))

noneqDSumAll = np.nansum(DenAnnAvg[0:4,...]*NHlats.mask,axis=(0,1,2)) + np.nansum(DenAnnAvg[0:4,...]*SHlats.mask,axis=(0,1,2))
nhDSumAll = np.nansum(DenAnnAvg[0:4,...]*NHlats.mask,axis=(0,1,2))
shDSumAll = np.nansum(DenAnnAvg[0:4,...]*SHlats.mask,axis=(0,1,2))
totalDSumAll = np.sum(DenAnnAvg[0:4,...],axis=(0,1,2))

print "Fractions relative to totals for all categories, excluding equator values"
print "Total Precip:", 100.0 * np.nansum(PrecipAnnAvg,axis=(1,2))/totalSumAll
print "NH Precip:   ", 100.0 * np.nansum(PrecipAnnAvg*NHlats.mask,axis=(1,2))/noneqSumAll
print "SH Precip:   ", 100.0 * np.nansum(PrecipAnnAvg*SHlats.mask,axis=(1,2))/noneqSumAll
print "Eq Precip:   ", 100.0 * np.nansum(PrecipAnnAvg*EQlats.mask,axis=(1,2))/totalSumAll

print "Land Precip: ", 100.0 * np.nansum(PrecipAnnAvg*Lmask.mask,axis=(1,2))/totalSumAll
print "Sea Precip:  ", 100.0 * np.nansum(PrecipAnnAvg*Smask.mask,axis=(1,2))/totalSumAll
print " "

np.set_printoptions(precision=1)
print "Total Density:", 100.0 * np.nansum(DenAnnAvg,axis=(1,2))/totalDSumAll
print "NH Density:  ", 100.0 * np.nansum(DenAnnAvg*NHlats.mask,axis=(1,2))/noneqDSumAll
print "SH Density:  ", 100.0 * np.nansum(DenAnnAvg*SHlats.mask,axis=(1,2))/noneqDSumAll
print "EQ Density:  ", 100.0 * np.nansum(DenAnnAvg*EQlats.mask,axis=(1,2))/totalDSumAll

print "Land Density:", 100.0 * np.nansum(DenAnnAvg*Lmask.mask,axis=(1,2))/totalDSumAll
print "Sea Density: ", 100.0 *np.nansum(DenAnnAvg*Smask.mask,axis=(1,2))/totalDSumAll




exit()

DenAll = np.nansum(DenAnnAvg,axis=0)
PrecipAll = np.nansum(PrecipAnnAvg,axis=0)

exit()


#PrecipAllPercent = PrecipAll/365.0	# Convert from mm/year to mm/day
PrecipAllPercent = 100.0 * np.divide(PrecipAll,PrecipIn*365.0)	# convert PrecipIn from mm/day to mm/yr before dividing

DenPercent = np.zeros(DenAnnAvg.shape)
PrecipPercent = np.zeros(PrecipAnnAvg.shape)

titlesDen = []
titlesPrec = []

DenPercent[0,:,:] = DenAll / (difflat * difflon)
PrecipPercent[0,:,:] = PrecipAllPercent	# Percent of total TRMM precip falling in these events

titlesDen.append("Total annual event density, events/(yr deg~S1~2 )")
titlesPrec.append("Percentage of total precipitation captured in events")

for iday in range(0,ndayplots):
	DenPercent[iday+1,:,:] = 100.0 * (np.where(DenAll > 0, np.divide(DenAnnAvg[iday,:,:],DenAll),0.0)) 
        PrecipPercent[iday+1,:,:] = 100.0 * (np.where(PrecipAll > 0, np.divide(PrecipAnnAvg[iday,:,:],PrecipAll),0.0))


# NH: find lat of 0 

print DenAnnAvg.shape





