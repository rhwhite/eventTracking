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

day1 = 2
day2 = 5

FillValue = -9999
Data = "TRMM"
mapping = 'centre'
Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'


if day1 == 1 and day2 == 2:
	plotmin1 = [0,-1,0,-5000,29,-0.1,0,-2.0]
	plotmax1 = [40,1,200000,5000,34,0.1,100,2.0]
	plotspace1 = [4,0.2,20000,1000,0.5,0.02,10,0.4]

elif day1 == 2 and day2 == 5:
        plotmin1 = [0,-0.2,0,-2500,50,-0.2,0,-5.0]
        plotmax1 = [8,0.2,200000,2500,72,0.2,200,5.0]
        plotspace1 = [0.8,0.04,20000,500,3,0.04,20,1.0]

MinLonF = 0
MaxLonF = 360
MinLatF = -45
MaxLatF = 45

sumlats = 16
sumlons = 16

# Time period for analysis
astartyr = 1998
aendyr = 2014
nyears = aendyr - astartyr + 1
print nyears
FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2014
	PrecipDir = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/' 
	Precipin = "Regrid_" + str(sumlats) + "_" + str(sumlons) + "_SeasAnn_TRMM_" + str(startyr) + "-" + str(endyr) + "_3B42_3hrly_nonan.nc"
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

	Allin = 'DenDirSpd_Map_monthly_regrid_TRMM_'+ str(startyr) + "-" + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_Standard.nc'

elif Data == "ERAI":
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
elif Data == "CESM":
	startyr = 1990 # Don't change - tied to file names!
	endyr = 2014
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'


FileIn = 'DenDirSpd_Map_' + str(day1) + '-' + str(day2) + 'Day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

FileInP = 'PrecipDen_Map_' + str(day1) + '-' + str(day2) + 'Day_' + mapping + '_' + Data + '_' + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

#Get lons and lats
print DirIn + FileIn
FileIn = xray.open_dataset(DirIn + FileIn)
FileInAll = xray.open_dataset(DirIn + Allin)
FileInP = xray.open_dataset(DirIn + FileInP)

FileInPrecip = xray.open_dataset(PrecipDir + Precipin)

lats = FileIn['Latitude']
lons = FileIn['Longitude']

latsAll = FileInAll['Latitude']
lonsAll = FileInAll['Longitude']

nlats = len(lats)
nlons = len(lons)

if nlats != len(latsAll) or nlons != len(lonsAll):
	print nlats, len(latsAll), nlons, len(lonsAll)
	sys.error("problem with inputs")

#Seas = FileIn['seas']

def regressmaps(m,c,r,p,stderr,A,linreg,nlats,nlons):
# Calculate linear trend
	for ilat in range(0,nlats):
		for ilon in range(0,nlons):
			linregin = np.ma.masked_array(linreg[:,ilat,ilon],mask=np.isnan(linreg[:,ilat,ilon])).compressed()
			Ain = np.ma.masked_array(A,mask=np.isnan(linreg[:,ilat,ilon])).compressed()
			if len(Ain) < len(A)*0.6:
				m[ilat,ilon] = np.nan
			else:
				m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(Ain,linregin)

def plotmap(plotvarsAnn,plotmins,plotmaxs,plotspaces,vartitles,title,figtitle,lons,lats,minlon,maxlon,minlat,maxlat):

	nplots = len(plotvarsAnn)
	
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

	res.sfXCStartV = float(lons[0])
	res.sfXCEndV = float(lons[len(lons)-1])
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

#	plotvarsAnn[np.where(np.isnan(plotvarsAnn))] = FillValue
	for iplot in range(0,nplots):
		temp =plotvarsAnn[iplot]
		temp[np.where(np.isnan(temp))] = FillValue	
		res.cnMinLevelValF       = plotmins[iplot]          # contour levels.
		res.cnMaxLevelValF       = plotmaxs[iplot]
		res.cnLevelSpacingF      = plotspaces[iplot]
		res.tiMainString = "Annual; " + vartitles[iplot]
		toplot.append(Ngl.contour_map(wks,plotvarsAnn[iplot],res))

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

	plot = Ngl.panel(wks,toplot,[nplots/2,2],panelres)



# And now plot 
c = np.zeros((4,nlats,nlons),np.float)
m = np.zeros((4,nlats,nlons),np.float)
r = np.zeros((4,nlats,nlons),np.float)
p = np.zeros((4,nlats,nlons),np.float)
stderr = np.zeros((4,nlats,nlons),np.float)

# Plot 1: Average density, and percentage easterly

figtitlein = FigDir + Data + '_' + Version + '_Annual_variables4_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

TDensity = FileIn['TDensityAnn']
linreg = TDensity[0:nyears,:,:]
A = np.array(range(0,nyears-1))
regressmaps(m[0,:,:],c[0,:,:],r[0,:,:],p[0,:,:],stderr[0,:,:],A,linreg,nlats,nlons)

TPrecip = FileInP['EPrecipAnn'] + FileInP['WPrecipAnn']
print np.nanmean(TPrecip)
linregP = TPrecip[0:nyears,:,:]
regnum = 1
regressmaps(m[regnum,:,:],c[regnum,:,:],r[regnum,:,:],p[regnum,:,:],stderr[regnum,:,:],A,linregP,nlats,nlons)

print np.sum(np.isfinite(m[0,:,:]))

TimespanAnn = np.nansum([FileIn['ETSpanAnn'].values, FileIn['WTSpanAnn'].values],axis=0)
TimespanAnn = np.where(TimespanAnn == 0,np.nan,TimespanAnn)
SpaceAnn = FileIn['ESizeAnn'].values + FileIn['WSizeAnn'].values 

linregT = np.divide(TimespanAnn[0:nyears,:,:],TDensity[0:nyears,:,:])
print np.sum(np.isfinite(linregT))
regnum = 2
regressmaps(m[regnum,:,:],c[regnum,:,:],r[regnum,:,:],p[regnum,:,:],stderr[regnum,:,:],A,linregT,nlats,nlons)
print np.sum(np.isfinite(m[2,:,:]))
print m.shape
TimePerEventAnn = np.nanmean(TimespanAnn/TDensity,axis=0)

SpacePerEvTS = np.nanmean(SpaceAnn/TimespanAnn,axis=0)
print np.nanmean(SpacePerEvTS)

linregS = np.divide(SpaceAnn[0:nyears,:,:],TimespanAnn[0:nyears,:,:])

regnum = 3
regressmaps(m[regnum,:,:],c[regnum,:,:],r[regnum,:,:],p[regnum,:,:],stderr[regnum,:,:],A,linregS,nlats,nlons)

print np.nanmean(m[3,:,:])



plotmap([np.mean(TDensity,axis=0),m[0,:,:],np.mean(TPrecip,axis=0),m[1,:,:],TimePerEventAnn* 3.0,m[2,:,:],SpacePerEvTS,m[3,:,:]],plotmin1,plotmax1,plotspace1,['Density of Events centered on gridbox','Regression of density','Average Precip from events, mm','Regression of precip','Average timespan per event, hours','Regression of timespan','Average gridboxes covered per event timestep','regression of gridboxes'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)


Ngl.end()




