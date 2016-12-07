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
astartyr = 1999
aendyr = 2015
nyears = aendyr - astartyr + 1
print nyears
FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2015
	if Version == 'Standard':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
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

lats = FileIn['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
lons = FileIn['lon'].values
years = FileIn['years']

print lats

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
	print DirIn + FileIn
	FileIn = xray.open_dataset(DirIn + FileIn)

	Dendays[iday,:,:,:] = FileIn['TDensityAnn'].sel(lat=slice(MinLatF,MaxLatF)) #np.nansum([FileIn['WDensityAnn'],FileIn['EDensityAnn'],FileIn['SDensityAnn']],axis=0)
	print day1,day2

	Precipdays[iday,:,:,:] = FileIn['TPrecipAnn'].sel(lat=slice(MinLatF,MaxLatF))


print Dendays.shape

DenAnnAvg = np.nanmean(Dendays,axis=1)
PrecipAnnAvg = np.nanmean(Precipdays,axis=1)

DenAll = np.nansum(DenAnnAvg,axis=0)
PrecipAll = np.nansum(PrecipAnnAvg,axis=0)

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

	if iday < ndayplots-1:
		titlesDen.append(str(day1[iday]) + ' to ' + str(day2[iday]) + ' day events')
	        titlesPrec.append(str(day1[iday]) + ' to ' + str(day2[iday]) + ' day events')
	else:
                titlesDen.append(">" + str(day1[iday]) + ' day events')
                titlesPrec.append(">" + str(day1[iday]) + ' day events')
	print np.nanmean(DenPercent[iday+1,:,:])
	print np.nanmean(PrecipPercent[iday+1,:,:])


# And now plot 

# Plot 1: Average density, and percentage easterly

figtitlein = FigDir + 'Paper1new_DenPrecipClim_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(day1[0]) + '_to_' + str(day2[ndayplots-1]) + 'Day'
titlein = Data #'Statistics of events in ' + Data + " " + Version

if Data == "TRMM":
	plotmap(DenPercent,PrecipPercent,[0.0,98.0,0.0,0.0,0.0],[300,100.0,1.0,0.2,0.05],[0.0,0.0,0.0,0.0,0.0],[100,90.0,40.0,40.0,40.0],titlesDen,titlesPrec,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)
elif Data == "ERAI":
	plotmap(DenPercent,PrecipPercent,[0.0,95.0,0.0,0.0,0.0],[30,100.0,2.0,1.5,0.5],[0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0],titlesDen,titlesPrec,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)
elif Data == "CESM":
        plotmap(DenPercent,PrecipPercent,[0.0,90.0,0.0,0.0,0.0],[50,100.0,5.0,3.0,0.5],[0.0,0.0,0.0,0.0,0.0],[20,90.0,40.0,40.0,40.0],titlesDen,titlesPrec,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

Ngl.end()




