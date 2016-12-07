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

day1 = 5
day2 = 10

FillValue = -9999
Data = "TRMM"
mapping = 'centre'
Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'


plotmin1 = 0
plotmax1 = 3
plotspace1 = 0.25

plotmin2 = -0.1
plotmax2 = 0.1
plotspace2 = 0.02

plotmin3 = 0
plotmax3 = 3
plotspace3 = 0.25

plotmin4 = -0.1
plotmax4 = 0.1
plotspace4 = 0.02


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

#Get lons and lats
print DirIn + FileIn
FileIn = xray.open_dataset(DirIn + FileIn)
FileInAll = xray.open_dataset(DirIn + Allin)

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
			m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(A,linreg[:,ilat,ilon])

def plotmap(plotvar1Ann,plotvar1Seas,plotvar2Ann,plotvar2Seas,plotmin,plotmax,plotspace,seasratio,vartitle,title,figtitle,lons,lats,minlon,maxlon,minlat,maxlat):
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

	plotvar1Ann[np.where(np.isnan(plotvar1Ann))] = FillValue
        plotvar2Ann[np.where(np.isnan(plotvar2Ann))] = FillValue
        plotvar1Seas[np.where(np.isnan(plotvar1Seas))] = FillValue
        plotvar2Seas[np.where(np.isnan(plotvar2Seas))] = FillValue

	res.cnMinLevelValF       = plotmin[0]          # contour levels.
	res.cnMaxLevelValF       = plotmax[0]
	res.cnLevelSpacingF      = plotspace[0]
	res.tiMainString = "Annual; " + vartitle[0]
	toplot.append(Ngl.contour_map(wks,plotvar1Ann,res))

        res.cnMinLevelValF       = plotmin[1]          # contour levels.
        res.cnMaxLevelValF       = plotmax[1]
        res.cnLevelSpacingF      = plotspace[1]
        res.tiMainString = "Annual; " + vartitle[1]
        toplot.append(Ngl.contour_map(wks,plotvar2Ann,res))

	for iseas in range(0,nseas):
		res.cnMinLevelValF       = plotmin[0]/seasratio[0]          # contour levels.
		res.cnMaxLevelValF       = plotmax[0]/seasratio[0]
		res.cnLevelSpacingF      = plotspace[0]/seasratio[0]
		res.tiMainString = Seas[iseas] + "; " + vartitle[0]
		toplot.append(Ngl.contour_map(wks,plotvar1Seas[iseas,:,:],res))
		
	        res.cnMinLevelValF       = plotmin[1]/seasratio[1]          # contour levels.
	        res.cnMaxLevelValF       = plotmax[1]/seasratio[1] 
	        res.cnLevelSpacingF      = plotspace[1]/seasratio[1] 
	        res.tiMainString = Seas[iseas] + "; " + vartitle[1]
	        toplot.append(Ngl.contour_map(wks,plotvar2Seas[iseas,:,:],res))

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

	plot = Ngl.panel(wks,toplot,[5,2],panelres)



# And now plot 

# Plot 1: Average density, and percentage easterly

figtitlein = FigDir + Data + '_' + Version + '_Density_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

EPercAnn = FileIn['EPercAnn'].values
EPercAnn[np.where(EPercAnn > 10E30)] = np.nan
EPercSeas = FileIn['EPercSeas'].values
EPercSeas[np.where(EPercSeas > 10E30)] = np.nan

plotmap(np.mean(FileIn['TDensityAnn'],axis=0),np.mean(FileIn['TDensitySeas'],axis=0),np.nanmean(EPercAnn,axis=0),np.nanmean(EPercSeas,axis=0),[0.0,0.0],[8.0,100.0],[0.8,10.0],[4,1],['Avg Number of Events centered on gridbox','Percentage of events easterly'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

# Plot 2: Westerly density, Stationary density
figtitlein = FigDir + Data + '_' + Version + '_West_Stat_Density_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version
plotmap(np.mean(FileIn['WDensityAnn'],axis=0),np.mean(FileIn['WDensitySeas'],axis=0),np.nanmean(FileIn['SDensityAnn'],axis=0),np.nanmean(FileIn['SDensitySeas'],axis=0),[0.0,0.0],[8.0,8.0],[0.8,0.8],[4,4],['Avg Number of Westerly Events centered on gridbox','Avg Number of Stationary events centered on gridbox'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)


# Plot 3: Precip vs Precip in 2-5 day events
figtitlein = FigDir + Data + '_' + Version + '_Precip_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

PrecipInraw1 = FileInPrecip['PrecipAnnClim']
PrecipInraw2 = FileInPrecip['PrecipSeasClim']

# Currently data has been summed over all events per season, but not divided by the number of 3hrly periods per season, hence factor of 91*8
# Annual is average of seasonal data, so factor is there same
PrecipIn1 = (np.nanmean(FileIn['SPrecipAnn'],axis=0) + np.nanmean(FileIn['WPrecipAnn'],axis=0) + np.nanmean(FileIn['EPrecipAnn'],axis=0))/(91*8)
PrecipIn2 = (np.nanmean(FileIn['SPrecipSeas'],axis=0) + np.nanmean(FileIn['WPrecipSeas'],axis=0) + np.nanmean(FileIn['EPrecipSeas'],axis=0))/(91*8)

plotmap(PrecipInraw1,PrecipInraw2,PrecipIn1,PrecipIn2,[0.0,0.0],[0.5,0.25],[0.05,0.025],[1,1],['Avg Precip mm','Avg Precip in ' + str(day1) + '-' + str(day2) + ' events, %'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

# Plot average spatial area and average timespan
figtitlein = FigDir + Data + '_' + Version + '_Timespan_Spatialspan_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'

titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

DensityAnn = FileIn['TDensityAnn']
DensitySeas = FileIn['TDensitySeas']

nyears = DensityAnn.shape[0]
print nyears

TimespanAnn = FileIn['TimespanAnn']
TimespanSeas = FileIn['TimespanSeas']

SpaceAnn = FileIn['GridboxspanAnn']
SpaceSeas = FileIn['GridboxspanSeas']

TimePerEventAnn = np.nansum(TimespanAnn/DensityAnn,axis=0)/nyears
TimePerEventSeas = np.nansum(TimespanSeas/DensityAnn,axis=0)/nyears



plotmap(,axis=0),np.mean(FileIn['WDensitySeas'],axis=0),np.nanmean(FileIn['SDensityAnn'],axis=0),np.nanmean(FileIn['SDensitySeas'],axis=0),[0.0,0.0],[8.0,8.0],[0.8,0.8],[4,4],['Avg Number of Westerly Events centered on gridbox','Avg Number of Stationary events centered on gridbox'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)



"""
# Plot 4 Precip and percentage of precip
figtitlein = FigDir + Data + '_' + Version + '_Precip_percent_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

plotmap(PrecipInraw1,PrecipInraw2,100.0*PrecipIn1/PrecipInraw1,100.0*PrecipIn2/PrecipInraw2,[0.0,0.0],[0.5,30.0],[0.05,3.0],[1,1],['Avg Precip mm','% of Precip in ' + str(day1) + '-' + str(day2) + ' events, %'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

# Plot 5 Precip and percentage of precip in easterly 2-5 day events
figtitlein = FigDir + Data + '_' + Version + '_EPrecip_percent_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day easterly events in ' + Data + " " + Version
# Currently data has been summed over all events per season, but not divided by the number of 3hrly periods per season, hence factor of 91*8
# Annual is average of seasonal data, so factor is there same
PrecipInE1 = np.nanmean(FileIn['EPrecipAnn'],axis=0)/(91*8)
PrecipInE2 = np.nanmean(FileIn['EPrecipSeas'],axis=0)/(91*8)

plotmap(PrecipInraw1,PrecipInraw2,100.0*PrecipInE1/PrecipInraw1,100.0*PrecipInE2/PrecipInraw2,[0.0,0.0],[0.5,30.0],[0.05,3.0],[1,1],['Avg Precip mm','% of Precip in ' + str(day1) + '-' + str(day2) + ' easterly events, %'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

# Plot 5 Precip and percentage of precip in easterly 2-5 day events
figtitlein = FigDir + Data + '_' + Version + '_WPrecip_percent_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day westerly events in ' + Data + " " + Version
# Currently data has been summed over all events per season, but not divided by the number of 3hrly periods per season, hence factor of 91*8
# Annual is average of seasonal data, so factor is there same
PrecipInW1 = np.nanmean(FileIn['WPrecipAnn'],axis=0)/(91*8)
PrecipInW2 = np.nanmean(FileIn['WPrecipSeas'],axis=0)/(91*8)

plotmap(PrecipInraw1,PrecipInraw2,100.0*PrecipInW1/PrecipInraw1,100.0*PrecipInW2/PrecipInraw2,[0.0,0.0],[0.5,30.0],[0.05,3.0],[1,1],['Avg Precip mm','% of Precip in ' + str(day1) + '-' + str(day2) + ' westerly events, %'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

# Now plot % of all events that are 2-5 days
figtitlein = FigDir + Data + '_' + Version + '_RelativeDensity_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version


# Sum over sizes and then mean over years
# In seasonal, sizes is axis 2

test = np.sum(FileInAll['DensityMapAnn'],axis=(1))
print np.nansum(test,axis=(1,2))

DensityAllA = np.nanmean(np.sum(FileInAll['DensityMapAnn'],axis=(1)),axis=0)
DensityAllS = np.nanmean(np.sum(FileInAll['DensityMapSeas'],axis=(2)),axis=0)

globmean = np.nanmean(DensityAllA)
print globmean

# For now, switch months
print "SWITCHING SEASONS BECAUSE NEW SEASONS START AT MAM"
DensityAllStemp = DensityAllS[...]
DensityAllS[0,:,:] = DensityAllStemp[1,:,:]
DensityAllS[1,:,:] = DensityAllStemp[2,:,:]
DensityAllS[2,:,:] = DensityAllStemp[3,:,:]
DensityAllS[3,:,:] = DensityAllStemp[0,:,:]

DensityAllA

PercDenA = 100.0*np.mean(FileIn['TDensityAnn'],axis=0)/DensityAllA
PercDenS = 100.0*np.mean(FileIn['TDensitySeas'],axis=0)/DensityAllS

print "making nan"
for ilat in range(0,nlats):
	for ilon in range(0,nlons):
		if DensityAllA[ilat,ilon] < globmean * 0.05:
			PercDenA[ilat,ilon] = FillValue

		for iseas in range(0,nseas):
			if DensityAllS[iseas,ilat,ilon] < globmean * 0.05:
				PercDenS[iseas,ilat,ilon] = FillValue

print "plotting"
plotmap(DensityAllA,DensityAllS,PercDenA,PercDenS,[0,0],[5000,0.2],[500,0.02],[4,1],['Average number of events (all sizes) per 4x4degree box','% of Events ' + str(day1) + '-' + str(day2) + ' day events, masked %'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

"""
"""
# Plot 4: Regression of Precip, regression of 2-5 day events
c = np.zeros((4,nseas+1,nlats,nlons),np.float)
m = np.zeros((4,nseas+1,nlats,nlons),np.float)
r = np.zeros((4,nseas+1,nlats,nlons),np.float)
p = np.zeros((4,nseas+1,nlats,nlons),np.float)
stderr = np.zeros((4,nseas+1,nlats,nlons),np.float)

A = np.array(range(0,nyears-2))

for iseas in range(0,nseas):
	linreg = FileIn["TDensitySeas"][0:nyears-2,iseas,:,:]
	regressmaps(m[0,iseas,:,:],c[0,iseas,:,:],r[0,iseas,:,:],p[0,iseas,:,:],stderr[0,iseas,:,:],A,linreg,nlats,nlons)

        #linreg = FileIn["EDensitySeas"][:,iseas,:,:]
        #regressmaps(m[1,iseas,:,:],c[1,iseas,:,:],r[1,iseas,:,:],p[1,iseas,:,:],stderr[1,iseas,:,:],A,linreg,nlats,nlons)

        linreg = FileInPrecip["PrecipSeas"][0:nyears-2,iseas,:,:]
        regressmaps(m[2,iseas,:,:],c[2,iseas,:,:],r[2,iseas,:,:],p[2,iseas,:,:],stderr[2,iseas,:,:],A,linreg,nlats,nlons)

#Annual 
linreg = FileIn["TDensityAnn"][0:nyears-2,:,:]
regressmaps(m[0,nseas,:,:],c[0,nseas,:,:],r[0,nseas,:,:],p[0,nseas,:,:],stderr[0,nseas,:,:],A,linreg,nlats,nlons)
#linreg = FileIn["EDensityAnn"][:,:,:]
#regressmaps(m[1,nseas,:,:],c[1,nseas,:,:],r[1,nseas,:,:],p[1,nseas,:,:],stderr[1,nseas,:,:],A,linreg,nlats,nlons)
linreg = FileInPrecip["PrecipAnn"][0:nyears-2,:,:]
regressmaps(m[2,nseas,:,:],c[2,nseas,:,:],r[2,nseas,:,:],p[2,nseas,:,:],stderr[2,nseas,:,:],A,linreg,nlats,nlons)

figtitlein = FigDir + Data + '_' + Version + '_DensityRegression_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

plotmap(m[2,4,:,:],m[2,0:4,:,:],m[0,4,:,:],m[0,0:4,:,:],[-0.01,-0.1],[0.01,0.1],[0.002,.02],[1,1],['Linear regression of total Precip','Linear regression of number of ' + str(day1) + '-' + str(day2) + ' day events'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

"""
# Plot regression in precip from 2-5 day events

c = np.zeros((4,nseas+1,nlats,nlons),np.float)
m = np.zeros((4,nseas+1,nlats,nlons),np.float)
r = np.zeros((4,nseas+1,nlats,nlons),np.float)
p = np.zeros((4,nseas+1,nlats,nlons),np.float)
stderr = np.zeros((4,nseas+1,nlats,nlons),np.float)

A = np.array(range(0,nyears-2))

for iseas in range(0,nseas):
        linreg = FileIn["EPrecipSeas"][0:nyears-2,iseas,:,:] + FileIn["WPrecipSeas"][0:nyears-2,iseas,:,:]
        regressmaps(m[0,iseas,:,:],c[0,iseas,:,:],r[0,iseas,:,:],p[0,iseas,:,:],stderr[0,iseas,:,:],A,linreg,nlats,nlons)

        #linreg = FileIn["EDensitySeas"][:,iseas,:,:]
        #regressmaps(m[1,iseas,:,:],c[1,iseas,:,:],r[1,iseas,:,:],p[1,iseas,:,:],stderr[1,iseas,:,:],A,linreg,nlats,nlons)

        linreg = FileInPrecip["PrecipSeas"][0:nyears-2,iseas,:,:]
        regressmaps(m[2,iseas,:,:],c[2,iseas,:,:],r[2,iseas,:,:],p[2,iseas,:,:],stderr[2,iseas,:,:],A,linreg,nlats,nlons)

#Annual 
linreg = FileIn["EPrecipAnn"][0:nyears-2,:,:] + FileIn["WPrecipAnn"][0:nyears-2,:,:]
regressmaps(m[0,nseas,:,:],c[0,nseas,:,:],r[0,nseas,:,:],p[0,nseas,:,:],stderr[0,nseas,:,:],A,linreg,nlats,nlons)
#linreg = FileIn["EDensityAnn"][:,:,:]
#regressmaps(m[1,nseas,:,:],c[1,nseas,:,:],r[1,nseas,:,:],p[1,nseas,:,:],stderr[1,nseas,:,:],A,linreg,nlats,nlons)
linreg = FileInPrecip["PrecipAnn"][0:nyears-2,:,:]
regressmaps(m[2,nseas,:,:],c[2,nseas,:,:],r[2,nseas,:,:],p[2,nseas,:,:],stderr[2,nseas,:,:],A,linreg,nlats,nlons)
figtitlein = FigDir + Data + '_' + Version + '_PrecipRegression_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

plotmap(m[2,4,:,:],m[2,0:4,:,:],m[0,4,:,:],m[0,0:4,:,:],[-0.01,-5],[0.01,5],[0.002,1],[1,1],['Linear regression of total Precip','Linear regression of precip in ' + str(day1) + '-' + str(day2) + ' day events'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)



# Plot regression of total precip from 2-5 day events

c = np.zeros((4,nseas+1,nlats,nlons),np.float)
m = np.zeros((4,nseas+1,nlats,nlons),np.float)
r = np.zeros((4,nseas+1,nlats,nlons),np.float)
p = np.zeros((4,nseas+1,nlats,nlons),np.float)
stderr = np.zeros((4,nseas+1,nlats,nlons),np.float)

A = np.array(range(0,nyears-2))

for iseas in range(0,nseas):
        linreg = FileIn["ETotalPrecipSeas"][0:nyears-2,iseas,:,:] + FileIn["WTotalPrecipSeas"][0:nyears-2,iseas,:,:]
        regressmaps(m[0,iseas,:,:],c[0,iseas,:,:],r[0,iseas,:,:],p[0,iseas,:,:],stderr[0,iseas,:,:],A,linreg,nlats,nlons)

        #linreg = FileIn["EDensitySeas"][:,iseas,:,:]
        #regressmaps(m[1,iseas,:,:],c[1,iseas,:,:],r[1,iseas,:,:],p[1,iseas,:,:],stderr[1,iseas,:,:],A,linreg,nlats,nlons)

        linreg = FileInPrecip["PrecipSeas"][0:nyears-2,iseas,:,:]
        regressmaps(m[2,iseas,:,:],c[2,iseas,:,:],r[2,iseas,:,:],p[2,iseas,:,:],stderr[2,iseas,:,:],A,linreg,nlats,nlons)

#Annual 
linreg = FileIn["ETotalPrecipAnn"][0:nyears-2,:,:] + FileIn["WTotalPrecipAnn"][0:nyears-2,:,:]
regressmaps(m[0,nseas,:,:],c[0,nseas,:,:],r[0,nseas,:,:],p[0,nseas,:,:],stderr[0,nseas,:,:],A,linreg,nlats,nlons)
#linreg = FileIn["EDensityAnn"][:,:,:]
#regressmaps(m[1,nseas,:,:],c[1,nseas,:,:],r[1,nseas,:,:],p[1,nseas,:,:],stderr[1,nseas,:,:],A,linreg,nlats,nlons)
linreg = FileInPrecip["PrecipAnn"][0:nyears-2,:,:]
regressmaps(m[2,nseas,:,:],c[2,nseas,:,:],r[2,nseas,:,:],p[2,nseas,:,:],stderr[2,nseas,:,:],A,linreg,nlats,nlons)
figtitlein = FigDir + Data + '_' + Version + '_TotalPrecipRegression_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'Day'
titlein = 'Statistics of ' + str(day1) + '-' + str(day2) + 'day events in ' + Data + " " + Version

plotmap(m[2,4,:,:],m[2,0:4,:,:],m[0,4,:,:],m[0,0:4,:,:],[-0.01,-2E12],[0.01,2E12],[0.002,4E11],[1,1],['Linear regression of total Precip','Linear regression of total precip in ' + str(day1) + '-' + str(day2) + ' day events'],titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)



Ngl.end()




