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

tbound1 = [-30,-10,-6,-3,3,6,10]
tbound2 = [-10,-6,-3,3,6,10,30]
unit = "ms"
splittype = "maxspeed"      #"days" "speed","maxspeed"
speedtspan = 4
unittitle = "m/s"
ndayplots = len(tbound1)

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


def regressmaps(m,c,r,p,stderr,A,linreg,nlats,nlons):
# Calculate linear trend
	for ilat in range(0,nlats):
		for ilon in range(0,nlons):
			m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(A,linreg[:,ilat,ilon])


def shiftlons(invar,inlons):
	nlons = inlons.shape[0]
	nlonhalf = nlons/2

	newinvar = np.zeros(invar.shape,np.float)
	newinvar[:,0:nlonhalf] = invar[:,nlonhalf:nlons]
	newinvar[:,nlonhalf:nlons] = invar[:,0:nlonhalf]
	return newinvar

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
		globmean = np.nanmean(tempplot)
		tempplot[np.where(np.isnan(tempplot))] = FillValue
		res.cnMinLevelValF       = plotmin1[iplot]          # contour levels.
		res.cnMaxLevelValF       = plotmax1[iplot]
		res.cnLevelSpacingF      = ((plotmax1[iplot]-plotmin1[iplot])/10.0)
		res.tiMainString = vartitle1[iplot] + "; mean = {:2.3g}".format(globmean) + "%"
		toplot.append(Ngl.contour_map(wks,tempplot,res))

                tempplot = plotvars2[iplot]
		globmean = np.nanmean(tempplot)
                tempplot[np.where(np.isnan(tempplot))] = FillValue
                res.cnMinLevelValF       = plotmin2[iplot]          # contour levels.
                res.cnMaxLevelValF       = plotmax2[iplot]
                res.cnLevelSpacingF      = ((plotmax2[iplot]-plotmin2[iplot])/10.0)
                res.tiMainString = vartitle2[iplot] + "; mean = {:2.3g}".format(globmean) + "%"
                toplot.append(Ngl.contour_map(wks,tempplot,res))

	
	#textres = Ngl.Resources()
	#textres.txFontHeightF = 0.015
	#Ngl.text_ndc(wks,title,0.5,0.87,textres)


	panelres = Ngl.Resources()
	panelres.nglPanelLabelBar = True
	#panelres.nglPanelYWhiteSpacePercent = 5.
	#panelres.nglPanelXWhiteSpacePercent = 5.

	panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
	if nplots > 6:
                panelres.nglPanelTop                      = 0.75
                panelres.nglPanelBottom                      = 0.1
	elif nplots > 5:
		panelres.nglPanelTop                      = 0.9
		panelres.nglPanelBottom                      = 0.1
	
	panelres.nglPanelFigureStrings            = ["a","b","c","d","e","f","g","h","i","j"]
	panelres.nglPanelFigureStringsJust        = "TopLeft"

	panelres.nglPaperOrientation = "Auto"

	plot = Ngl.panel(wks,toplot,[nplots,2],panelres)

# Get lats and lons
iloop = 0
if speedtspan == 0:
	fileadd = ""
else:
	fileadd = str(speedtspan) + "ts_"

tboundtitle = str(int(tbound1[iloop])) + '-' + str(int(tbound2[iloop]))
if splittype == "maxspeed":
	fileadd = fileadd + "MaxSpeeds_"
elif splittype == "speed":
	fileadd = fileadd + "" 
elif splittype == "days":
	fileadd = fileadd + "Sizes_"

FileIn = 'DenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version +'_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

print DirIn + FileIn
FileIn = xray.open_dataset(DirIn + FileIn)

lats = FileIn['Latitude']
lons = FileIn['Longitude']
years = FileIn['years']

nlats = len(lats)
nlons = len(lons)
nyears = len(years)

Dendays = np.zeros([ndayplots,nyears,nlats,nlons])
Precipdays = np.zeros([ndayplots,nyears,nlats,nlons])

for iloop in range(0,ndayplots):

        tboundtitle = str(int(tbound1[iloop])) + '-' + str(int(tbound2[iloop]))
        
	FileIn = 'DenDirSpd_Map_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version +'_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

	#Get lons and lats
	print DirIn + FileIn
	FileIn = xray.open_dataset(DirIn + FileIn)

	Dendays[iloop,:,:,:] = np.nansum([FileIn['EDensityAnn'],FileIn['WDensityAnn'],FileIn['SDensityAnn']],axis=0)

	if splittype == 'maxspeed' or splittype == 'speed':
		Precipdays[iloop,:,:,:] = FileIn['TPrecipAnn']
	else:	
		Precipdays[iloop,:,:,:] = np.nansum([FileIn['EPrecipAnn'], FileIn['WPrecipAnn']],axis=0)


print Dendays.shape

DenAnnAvg = np.nanmean(Dendays,axis=1)
PrecipAnnAvg = np.nanmean(Precipdays,axis=1)

DenAll = np.nansum(DenAnnAvg,axis=0)
PrecipAll = np.nansum(PrecipAnnAvg,axis=0)

DenPercent = np.zeros(DenAnnAvg.shape)
PrecipPercent = np.zeros(PrecipAnnAvg.shape)
titlesDen = []
titlesPrec = []
for iday in range(0,ndayplots):
	DenPercent[iday,:,:] = 100.0 * np.divide(DenAnnAvg[iday,:,:],DenAll) 
        PrecipPercent[iday,:,:] = 100.0 * np.divide(PrecipAnnAvg[iday,:,:],PrecipAll)

	if iday == 0:
	        titlesDen.append('< ' + str(tbound2[iday]) + unittitle)
	        titlesPrec.append('< ' + str(tbound2[iday]) + unittitle)
	elif iday == ndayplots-1:
	        titlesDen.append('> ' + str(tbound1[iday]) + unittitle)
	        titlesPrec.append('> ' + str(tbound1[iday]) + unittitle)

	titlesDen.append(str(tbound1[iday]) + ' to ' + str(tbound2[iday]) + unittitle)
        titlesPrec.append(str(tbound1[iday]) + ' to ' + str(tbound2[iday]) + unittitle)

	print np.nanmean(DenPercent[iday,:,:])
	print np.nanmean(PrecipPercent[iday,:,:])


# And now plot 

# Plot 1: Average density, and percentage easterly

figtitlein = FigDir + 'Paper4_DenPrecipClim_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[ndayplots-1]) + unit
titlein = 'Statistics of events in ' + Data + " " + Version

plotmap(DenPercent,PrecipPercent,[0.0,0.0,0.0,80.0,0.0,0.0,0.0],[1.0,2.0,10.0,100.0,10.0,2.0,1.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[30.0,30.0,30.0,30.0,30.0,30.0,30.0],titlesDen,titlesPrec,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)



Ngl.end()




