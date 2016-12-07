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
from rhwhitepackages.readwrite import shiftlons

day1 = 0
day2 = 1

FillValue = -9999
Data = "TRMM"
mapping = 'center'
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

elif day1 == 1 and day2 == 5:
        plotmin1 = [-1]
        plotmax1 = [1]
        plotspace1 = [0.2]

elif day1 == 2 and day2 == 5:
        plotmin1 = [-0.2]
        plotmax1 = [0.2]
        plotspace1 = [0.04]

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


def XrayOpen(filenamein,decodetimes=True):

        try:
                if decodetimes:
                        filein= xray.open_dataset(filenamein)
                else:
                        filein=xray.open_dataset(filenamein,decode_times=False)
        except RuntimeError:
                print filenamein
                exit("couldn't find file")
        return filein


if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2014
	PrecipDir = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/' 
	Precipin = "Regrid_" + str(sumlats) + "_" + str(sumlons) + "_SeasAnn_TRMM_" + str(startyr) + "-" + str(endyr) + "_3B42_3hrly_nonan.nc"
	if Version == 'Standard':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/TRMM_output/Standard' + str(startyr) + '/proc/'
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

	#Allin = 'DenDirSpd_Map_monthly_regrid_TRMM_'+ str(startyr) + "-" + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_Standard.nc'

elif Data == "ERAI":
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
elif Data == "CESM":
	startyr = 1990 # Don't change - tied to file names!
	endyr = 2014
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'

if day1 == 1 and day2 == 5:

	Filename1 = 'DenDirSpd_Map_Sizes_1-2day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

        Filename2 = 'DenDirSpd_Map_Sizes_2-5day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

	#Get lons and lats
	FileIn1 = XrayOpen(DirIn + Filename1)
        FileIn2 = XrayOpen(DirIn + Filename2)
	
	lats = FileIn1['Latitude']
	lons = FileIn1['Longitude']

	nlats = len(lats)
	nlons = len(lons)

	Dens = (np.nansum([FileIn1['WDensitySeas'],FileIn1['EDensitySeas'],FileIn1['SDensitySeas']],axis=0) + 
		np.nansum([FileIn2['WDensitySeas'],FileIn2['EDensitySeas'],FileIn2['SDensitySeas']],axis=0))	

else:
        Filename1 = 'DenDirSpd_Map_Ann_Sizes_' + str(day1) + '-' + str(day2) + 'day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

        FileIn1 = XrayOpen(DirIn + Filename1)

        lats = FileIn1['lat']
        lons = FileIn1['lon']

        nlats = len(lats)
        nlons = len(lons)

	Dens = FileIn1['TDensityAnn']	
        #Dens = (np.nansum([FileIn1['WDensitySeas'],FileIn1['EDensitySeas'],FileIn1['SDensitySeas']],axis=0))

nseas = Dens.shape[1]
print Dens.shape

Seasons = np.array(FileIn1["seas"])
print str(Seasons)


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

def plotmap(nplotseas,plotvarsSeas,plotmins,plotmaxs,plotspaces,vartitles,figtitle,lons,lats,minlon,maxlon,minlat,maxlat):
	
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

                for iplot in range(0,nplotseas):
                        plotvarsSeas[iplot,:,:] = shiftlons(plotvarsSeas[iplot,:,:],lons)
		lons = lonsnew

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

	print vartitles

#	plotvarsAnn[np.where(np.isnan(plotvarsAnn))] = FillValue
	for iseas in range(0,nplotseas):
		temp =plotvarsSeas[iseas,:,:]
		temp[np.where(np.isnan(temp))] = FillValue	
		res.cnMinLevelValF       = plotmins          # contour levels.
		res.cnMaxLevelValF       = plotmaxs
		res.cnLevelSpacingF      = plotspaces
		res.tiMainString = str(vartitles[iseas])
		toplot.append(Ngl.contour_map(wks,plotvarsSeas[iseas,:,:],res))

	#textres = Ngl.Resources()
	#textres.txFontHeightF = 0.015
	#Ngl.text_ndc(wks,title,0.5,0.87,textres)


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

	plot = Ngl.panel(wks,toplot,[nplotseas,1],panelres)



# And now plot 
c = np.zeros((4,nlats,nlons),np.float)
m = np.zeros((4,nlats,nlons),np.float)
r = np.zeros((4,nlats,nlons),np.float)
p = np.zeros((4,nlats,nlons),np.float)
stderr = np.zeros((4,nlats,nlons),np.float)

# Plot 1: Average density, and percentage easterly

figtitlein = FigDir + "Paper_denregress_" + Data + '_' + Version + '_Seas_' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + str(day1) + '-' + str(day2) + 'days'

for iseas in range(0,nseas):
	linreg = Dens[0:nyears,iseas,:,:]
	A = np.array(range(0,nyears-1))
	regressmaps(m[iseas,:,:],c[iseas,:,:],r[iseas,:,:],p[iseas,:,:],stderr[iseas,:,:],A,linreg,nlats,nlons)

plotmap(nseas,m,plotmin1,plotmax1,plotspace1,Seasons,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)


Ngl.end()




