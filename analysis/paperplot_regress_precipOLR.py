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
import pandas as pd
import xray
import Ngl
from scipy import stats
from scipy import interpolate
import math
from rhwhitepackages.readwrite import XrayOpen, shiftlons

Data = "TRMM"
OLRminlat = 16
OLRmaxlat = 56

daystart = 0
dayend = 1

Version = "Standard"

MinLonF = 0
MaxLonF = 359.5
MinLatF = -45
MaxLatF = 45

sumlats = 16
sumlons = 16

Esumlats = 5
Esumlons = 5

Csumlats = 16
Csumlons = 16


if daystart == 1 and dayend == 5:
        plotmin = [-15.0,-15.0,-15.0,-0.5]
        plotmax = [ 15.0, 15.0, 15.0,0.5]
        plotspace = [3.0,3.0,3.0,0.1]
elif daystart == 0 and dayend == 1:
        plotmin = [-15.0,-15.0,-15.0,-0.5]
        plotmax = [ 15.0, 15.0, 15.0,0.5]
        plotspace = [3.0,3.0,3.0,0.1]

Seas = ['DJF','MAM','JJA','SON','Ann']
nseas = 4
startyr = 1998 # Don't change - tied to file names!
endyr = 2014
Estartyr = 1980
Eendyr = 2014
Cstartyr = 1990
Cendyr = 2014
nyears = endyr - startyr-1

FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

# Get TRMM Data

DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/TRMM_output/Standard' + str(startyr) + '/proc/'
Filename = 'DenDirSpd_Map_Ann_Sizes_' + str(daystart) + '-' + str(dayend) + 'day_center_TRMM_' + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
FilenameR = 'DenDirSpd_Map_Ann_Sizes_' + str(daystart) + '-' + str(dayend) + 'day_center_TRMM_' + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_1lons_1lats.nc_regrid_OLR.nc'

TRMMIn = XrayOpen(DirIn + Filename)
TRMMInRG = XrayOpen(DirIn + FilenameR)
TRMMlat = TRMMIn['lat']
TRMMlon = TRMMIn['lon']

nlats = len(TRMMlat)
nlons = len(TRMMlon)
TRMMyears = TRMMIn['years'][0:nyears]
TRMMden = TRMMIn['TDensityAnn'][0:nyears,:,:].values
TRMMprecip = TRMMIn['TPrecipAnn'][0:nyears,:,:].values
TRMMdenOLR = TRMMInRG['TDensityAnn'][0:nyears,OLRminlat:OLRmaxlat,:].values
TRMMpreOLR = TRMMInRG['TPrecipAnn'][0:nyears,OLRminlat:OLRmaxlat,:].values

TRMMdenclim = np.nanmean(TRMMden,axis=(0,1,2))
print "TRMMdenclim", TRMMdenclim
TRMMprecipclim = np.nanmean(TRMMprecip,axis=(0,1,2))
print "annual progression:", np.nanmean(TRMMden,axis=(1,2))

print "nans: ", np.sum(np.isnan(TRMMpreOLR))

# Get ERAI Data
DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(Estartyr) + '/Precip/'


Filename = 'DenDirSpd_Map_Sizes_' + str(daystart) + '-' + str(dayend) + 'day_centre_ERAI_' + str(Estartyr) + '-' + str(Eendyr) + '_' + Version + '_regrid_' + str(Esumlons) + 'lons_' + str(Esumlats) + 'lats.nc'
FilenameR = 'DenDirSpd_Map_Sizes_' + str(daystart) + '-' + str(dayend) + 'day_centre_ERAI_' + str(Estartyr) + '-' + str(Eendyr) + '_' + Version + '_regrid_1lons_1lats.nc_regrid_OLR.nc'

ERAIn = XrayOpen(DirIn + Filename)
ERAInRG = XrayOpen(DirIn + FilenameR)
ERAyears = ERAIn['years']
if ERAyears[0] < 1998:
        for iyear in range(0,len(ERAyears)):
                if ERAyears[iyear] == 1998:
                        ERAsyear = iyear
                        break
else:
	ERAsyear = 0
if ERAyears[-1] > 2013:
	for iyear in range(ERAsyear,len(ERAyears)):
		if ERAyears[iyear] == 2013:
			ERAeyear = iyear
else:
	ERAeyear = len(ERAyears)-1	

ERAden = ERAIn['TDensityAnn'][ERAsyear:ERAeyear,...].values
ERAprecip = ERAIn['TPrecipAnn'][ERAsyear:ERAeyear,...].values
ERAdenOLR = ERAInRG['TDensityAnn'][ERAsyear:ERAeyear,OLRminlat:OLRmaxlat,...].values
ERApreOLR = ERAInRG['TPrecipAnn'][ERAsyear:ERAeyear,OLRminlat:OLRmaxlat,...].values
ERAyrs = ERAIn['years'][ERAsyear:ERAeyear]
ERAlat = ERAIn['lat']
ERAlon = ERAIn['lon']
Enlats = len(ERAlat)
Enlons = len(ERAlon)

ERAdenclim = np.mean(ERAden,axis=(0,1,2))
ERAprecipclim = np.mean(ERAprecip,axis=(0,1,2))


print ERApreOLR.shape
# Get CESM data

DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/Standard1990-2014/Precip/'

Filename = 'DenDirSpd_Map_Sizes_' + str(daystart) + '-' + str(dayend) + 'day_centre_CESM_' + str(Cstartyr) + '-' + str(Cendyr) + '_' + Version + '_regrid_' + str(Csumlons) + 'lons_' + str(Csumlats) + 'lats.nc'
FilenameR = 'DenDirSpd_Map_Sizes_' + str(daystart) + '-' + str(dayend) + 'day_centre_CESM_' + str(Cstartyr) + '-' + str(Cendyr) + '_' + Version + '_regrid_1lons_1lats.nc_regrid_OLR.nc'

CESMIn = XrayOpen(DirIn + Filename)
CESMInRG = XrayOpen(DirIn + FilenameR)
CESMyears = CESMIn['years']
if CESMyears[0] < 1998:
        for iyear in range(0,len(CESMyears)):
                if CESMyears[iyear] == 1998:
                        CESMsyear = iyear
                        break
else:   
        CESMsyear = 0
if CESMyears[-1] > 2012:
        for iyear in range(CESMsyear,len(CESMyears)):
                if CESMyears[iyear] == 2012:
                        CESMeyear = iyear
else:   
        CESMeyear = len(CESMyears)-1

CESMnyears = nyears-1

CESMden = CESMIn['TDensityAnn'][CESMsyear:CESMeyear,...].values
CESMprecip = CESMIn['TPrecipAnn'][CESMsyear:CESMeyear,...].values
CESMdenOLR = CESMInRG['TDensityAnn'][CESMsyear:CESMeyear,OLRminlat:OLRmaxlat,...].values
CESMpreOLR = CESMInRG['TPrecipAnn'][CESMsyear:CESMeyear,OLRminlat:OLRmaxlat,...].values
CESMyrs = CESMIn['years'][CESMsyear:CESMeyear]
CESMlat = CESMIn['lat']
CESMlon = CESMIn['lon']
Cnlats = len(CESMlat)
Cnlons = len(CESMlon)

CESMdenclim = np.mean(CESMden,axis=(0,1,2))
CESMprecipclim = np.mean(CESMprecip,axis=(0,1,2))


# Get OLR data
DirIn = "/home/disk/eos4/rachel/Obs/OLR/"
FileIn = "Ann_olr.mon.mean.nc" 

OLRin = XrayOpen(DirIn + FileIn)
nlatOLRT = len(OLRin["lat"])
print nlatOLRT
OLRyears = OLRin['Year']
if OLRyears[0] < 1998:
	for iyear in range(0,len(OLRyears)):
		if OLRyears[iyear] == 1998:
			OLRsyear = iyear
			break
else:
	OLRsyear = 0
if OLRyears[-1] > 2013:
        for iyear in range(0,len(OLRyears)):
                if OLRyears[iyear] == 2013:
                        OLReyear = iyear
                        break
else:
        OLReyear = len(OLRyears)-1

OLR = -1.0 * OLRin['OLRAnn'][OLRsyear:OLReyear,nlatOLRT-OLRminlat-1:nlatOLRT-OLRmaxlat-1:-1,...].values
OLRyrs = OLRin['Year'][OLRsyear:OLReyear,...]
OLRlon = OLRin['Longitude']
OLRlat = OLRin['Latitude'][nlatOLRT-OLRminlat-1:nlatOLRT-OLRmaxlat-1:-1]

Onlats = len(OLRlat)
Onlons = len(OLRlon)


figtitle = "Paper_precipmaptrends_" + Data + "-" + Version + "_" + str(startyr) + '-' + str(endyr) + '_'  + str(daystart) + '-' + str(dayend) + 'days'

print "OLR: " + str(len(OLRlon)) + " longitudes, " + str(len(OLRlat)) + " latitudes"
print "TRMM: " + str(len(TRMMlon)) + " longitudes, " + str(len(TRMMlat)) + " latitudes"

"""
#convert nans to missing values so plots work!
TRMMprecip._FillValue = -999.0
ERAprecip._FillValue = -999.0
OLR._FillValue = -999.0

TRMMprecip = np.where(np.isnan(TRMMprecip[:,:,:]),TRMMprecip._FillValue,TRMMprecip[:,:,:])
ERAprecip = np.where(np.isnan(ERAprecip[:,:,:]),ERAprecip._FillValue,ERAprecip[:,:,:])
OLR = np.where(np.isnan(OLR[:,:,:]),OLR._FillValue,OLR[:,:,:])
"""
# Regrid everything onto the lowest resolution grid, better hope there are no missing data points!

# Calculate linear trend over the 16 years 
# Calculate linear trend over the 16 years 
#yearnums = range(0,trendyears)
A = TRMMyears
#A = np.vstack([yearnums,np.ones(len(yearnums))]).T
#c = np.zeros((4,nlats,nlons),np.float)
#mTDRS = np.zeros((4,nlats*nlons),np.float)
#r = np.zeros((4,nlats,nlons),np.float)
#p = np.zeros((4,nlats,nlons),np.float)
#stderr = np.zeros((4,nlats,nlons),np.float)

A = np.vstack([TRMMyears,np.ones(len(TRMMyears))]).T
TRMMdenRS = TRMMden.reshape(nyears,nlats*nlons)
mTDRS,c = np.linalg.lstsq(A,TRMMdenRS)[0]
mTD = mTDRS.reshape(nlats,nlons)
del(c)

ERAdenRS = ERAden.reshape(nyears,Enlats*Enlons)
mEDRS,c = np.linalg.lstsq(A,ERAdenRS)[0]
mED = mEDRS.reshape(Enlats,Enlons)
del(c)

CESMA = np.vstack([CESMyrs,np.ones(len(CESMyrs))]).T
CESMdenRS = CESMden.reshape(CESMnyears,Cnlats*Cnlons)
mCDRS,c = np.linalg.lstsq(CESMA,CESMdenRS)[0]
mCD = mCDRS.reshape(Cnlats,Cnlons)
del(c)

print "TRMM:"
print np.sum(TRMMdenRS,axis=1)

TRMMprecipRS = TRMMprecip.reshape(nyears,nlats*nlons)
mTPRS,c = np.linalg.lstsq(A,TRMMprecipRS)[0]
mTP = mTPRS.reshape(nlats,nlons)
del(c)

ERAprecipRS = ERAprecip.reshape(nyears,Enlats*Enlons)
mEPRS,c = np.linalg.lstsq(A,ERAprecipRS)[0]
mEP = mEPRS.reshape(Enlats,Enlons)
del(c)

print "ERA:"
print np.sum(ERAprecipRS,axis=1)

CESMA = np.vstack([CESMyrs,np.ones(len(CESMyrs))]).T

CESMprecipRS = CESMprecip.reshape(CESMnyears,Cnlats*Cnlons)
mCPRS,c = np.linalg.lstsq(CESMA,CESMprecipRS)[0]
mCP = mCPRS.reshape(Cnlats,Cnlons)
del(c)

print "nyears: ", nyears
print "CESM"
print np.sum(CESMprecipRS,axis=1)
print "Done"
print OLR.shape
print nyears, Onlats, Onlons
OLRRS = OLR.reshape(nyears,Onlats*Onlons)
mORS,c = np.linalg.lstsq(A,OLRRS)[0]
mO = mORS.reshape(Onlats,Onlons)
del(c)
# Calculate correlation coefficients with re-gridded maps!
TRMMpreOLRRS = TRMMpreOLR.reshape(nyears,Onlats*Onlons)
TRMMdenOLRRS = TRMMdenOLR.reshape(nyears,Onlats*Onlons)
ERApreOLRRS = ERApreOLR.reshape(nyears,Onlats*Onlons)
ERAdenOLRRS = ERAdenOLR.reshape(nyears,Onlats*Onlons)
CESMpreOLRRS = CESMpreOLR.reshape(CESMnyears,Onlats*Onlons)
CESMdenOLRRS = CESMdenOLR.reshape(CESMnyears,Onlats*Onlons)


mTOLRRS = (np.linalg.lstsq(A,TRMMpreOLRRS)[0][0])
mTDOLRRS = (np.linalg.lstsq(A,TRMMdenOLRRS)[0][0])
mEOLRRS = (np.linalg.lstsq(A,ERApreOLRRS)[0][0])
mEDOLRRS = (np.linalg.lstsq(A,ERAdenOLRRS)[0][0])
mCOLRRS = (np.linalg.lstsq(CESMA,CESMpreOLRRS)[0][0])
mCDOLRRS = (np.linalg.lstsq(CESMA,CESMdenOLRRS)[0][0])


mTPOLR = mTOLRRS.reshape(Onlats,Onlons)
mEPOLR = mEOLRRS.reshape(Onlats,Onlons)

print "Pearsons correlations"

R = np.corrcoef(mORS,mEOLRRS)
print R
R = np.corrcoef(mORS,mTOLRRS)
print R

R,p = stats.pearsonr(mORS,mEOLRRS)
print "ERA precip with OLR, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.pearsonr(mORS,mTOLRRS)
print "TRMM precip with OLR, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.pearsonr(mORS,mCOLRRS)
print "CESM precip with OLR, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.pearsonr(mORS,mTDOLRRS)
print "TRMM density with OLR, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.pearsonr(mTDOLRRS,mEDOLRRS)
print "TRMM vs ERA density, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.pearsonr(mTDOLRRS,mCDOLRRS)
print "TRMM vs CESM density, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.pearsonr(mEDOLRRS,mCDOLRRS)
print "ERA vs CESM density, correlation coefficient: ", R, " with p-value of ", p


print "Spearman correlations:"
R,p = stats.spearmanr(mORS,mEOLRRS)
print "ERA precip, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.spearmanr(mORS,mTOLRRS)
print "TRMM precip, correlation coefficient: ", R, " with p-value of ", p

R,p = stats.spearmanr(mORS,mTDOLRRS)
print "TRMM density, correlation coefficient: ", R, " with p-value of ", p

# Calculate regression percentages:
print TRMMdenclim
print ERAdenclim
print CESMdenclim

mTD_perc = 100.0 * mTD/TRMMdenclim
mED_perc = 100.0 * mED/ERAdenclim
mCD_perc = 100.0 * mCD/CESMdenclim

mTD_perc[np.where(TRMMdenclim <= 0)] = -9999
mED_perc[np.where(ERAdenclim <= 0)] = -9999
mCD_perc[np.where(CESMdenclim <= 0)] = -9999

mTP_perc = 100.0 * mTP/TRMMprecipclim
mEP_perc = 100.0 * mEP/ERAprecipclim
mCP_perc = 100.0 * mCP/CESMprecipclim

mTP_perc[np.where(TRMMprecipclim <= 0)] = -9999
mEP_perc[np.where(ERAprecipclim <= 0)] = -9999
mCP_perc[np.where(CESMprecipclim <= 0)] = -9999

mO_perc = 100.0 * mO/np.nanmean(OLR,axis=0)

# Read in land-sea mask? How do we get one for TRMM? Generic one, regridded from other, with min 50% land?

# And now plot 
wkres = Ngl.Resources()
wkres.wkColorMap = "precip_diff_12lev"
wks_type = "eps"
wks = Ngl.open_wks(wks_type,FigDir + figtitle,wkres)

res = Ngl.Resources()
res.cnInfoLabelOn         = False    # Turn off informational
                                              # label.
res.pmLabelBarDisplayMode = "Always" # Turn on label bar.
res.cnLinesOn             = False    # Turn off contour lines.
 
res.nglDraw  = False
res.nglFrame = False

res.sfMissingValueV = -9999 


res.cnFillOn = True
res.cnMissingValFillColor = "white"
res.cnLineLabelsOn       = False
res.pmLabelBarDisplayMode = "Always"
res.cnLinesOn =  False


if TRMMlon[0] < 0:
	nlonhalf = nlons/2
	lonsnew = np.zeros(TRMMlon.shape,np.float)
	lonsnew[0:nlonhalf] = TRMMlon[nlonhalf:nlons]
	lonsnew[nlonhalf:nlons] = TRMMlon[0:nlonhalf] + 360.0
	TRMMlon = lonsnew

	mTD = shiftlons(mTD,TRMMlon)
	mTP = shiftlons(mTP,TRMMlon)
        mTD_perc = shiftlons(mTD_perc,TRMMlon)
	mTP_perc = shiftlons(mTP_perc,TRMMlon)
res.sfXCStartV = float(TRMMlon[0])
res.sfXCEndV = float(TRMMlon[len(TRMMlon)-1])
res.sfYCStartV = float(TRMMlat[0])
res.sfYCEndV = float(TRMMlat[len(TRMMlat)-1])

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

res.cnLevelSelectionMode = "ManualLevels" # Define your own

toplot = []

res.cnMinLevelValF       = plotmin[0]          # contour levels.
res.cnMaxLevelValF       = plotmax[0]
res.cnLevelSpacingF      = plotspace[0]
#res.tiMainString = Seas[iseas]

res.tiMainFontHeightF = 0.02
res.tiMainString = 'TRMM ' + str(daystart) + '-' + str(dayend) + 'day Density Annual regression map, % per year'
toplot.append(Ngl.contour_map(wks,mTD_perc,res))

res.sfXCStartV = float(ERAlon[0])
res.sfXCEndV = float(ERAlon[len(ERAlon)-1])
res.sfYCStartV = float(ERAlat[0])
res.sfYCEndV = float(ERAlat[len(ERAlat)-1])

res.cnMinLevelValF       = plotmin[1]          # contour levels.
res.cnMaxLevelValF       = plotmax[1]
res.cnLevelSpacingF      = plotspace[1]

res.tiMainString = 'ERA ' + str(daystart) + '-' + str(dayend) + 'day Density Annual regression map, % per year'
toplot.append(Ngl.contour_map(wks,mED_perc,res))

res.sfXCStartV = float(CESMlon[0])
res.sfXCEndV = float(CESMlon[len(CESMlon)-1])
res.sfYCStartV = float(CESMlat[0])
res.sfYCEndV = float(CESMlat[len(CESMlat)-1])

res.cnMinLevelValF       = plotmin[0]          # contour levels.
res.cnMaxLevelValF       = plotmax[0]
res.cnLevelSpacingF      = plotspace[0]

res.tiMainString = 'CESM ' + str(daystart) + '-' + str(dayend) + 'day Density Annual regression map, % per year'
toplot.append(Ngl.contour_map(wks,mCD_perc,res))

res.sfXCStartV = float(TRMMlon[0])
res.sfXCEndV = float(TRMMlon[len(TRMMlon)-1])
res.sfYCStartV = float(TRMMlat[0])
res.sfYCEndV = float(TRMMlat[len(TRMMlat)-1])

res.cnMinLevelValF       = plotmin[2]          # contour levels.
res.cnMaxLevelValF       = plotmax[2]
res.cnLevelSpacingF      = plotspace[2]

res.tiMainString = 'TRMM' + str(daystart) + '-' + str(dayend) + 'day Precip Annual regression map, % per year'
toplot.append(Ngl.contour_map(wks,mTP_perc,res))

res.sfXCStartV = float(ERAlon[0])
res.sfXCEndV = float(ERAlon[len(ERAlon)-1])
res.sfYCStartV = float(ERAlat[0])
res.sfYCEndV = float(ERAlat[len(ERAlat)-1])

res.cnMinLevelValF       = plotmin[2]          # contour levels.
res.cnMaxLevelValF       = plotmax[2]
res.cnLevelSpacingF      = plotspace[2]

res.tiMainString = 'ERA ' + str(daystart) + '-' + str(dayend) + 'day Precip Annual regression map, % per year'
toplot.append(Ngl.contour_map(wks,mEP_perc,res))

res.sfXCStartV = float(CESMlon[0])
res.sfXCEndV = float(CESMlon[len(CESMlon)-1])
res.sfYCStartV = float(CESMlat[0])
res.sfYCEndV = float(CESMlat[len(CESMlat)-1])

res.cnMinLevelValF       = plotmin[2]          # contour levels.
res.cnMaxLevelValF       = plotmax[2]
res.cnLevelSpacingF      = plotspace[2]

res.tiMainString = 'CESM ' + str(daystart) + '-' + str(dayend) + 'day Precip Annual regression map, % per year'
toplot.append(Ngl.contour_map(wks,mCP_perc,res))



res.sfXCStartV = float(OLRlon[0])
res.sfXCEndV = float(OLRlon[len(OLRlon)-1])
res.sfYCStartV = float(OLRlat[0])
res.sfYCEndV = float(OLRlat[len(OLRlat)-1])

res.cnMinLevelValF       = plotmin[3]          # contour levels.
res.cnMaxLevelValF       = plotmax[3]
res.cnLevelSpacingF      = plotspace[3]

res.tiMainString = "OLR Annual regression map"
toplot.append(Ngl.contour_map(wks,mO_perc,res))

#plot = Ngl.contour(wks,var,res)
panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
#panelres.nglPanelYWhiteSpacePercent = 5.
#panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.95
panelres.nglPanelBottom                      = 0.01

panelres.nglPanelFigureStrings            = ["a","b","c","d","e"]
panelres.nglPanelFigureStringsJust        = "BottomLeft"

#
# You can have PyNGL selection the best paper orientation for
# the shape of plots you are drawing.  This resource is for PDF or
# PS output only.
#
panelres.nglPaperOrientation = "Auto"

#labels = ["High " + fill + " U","Mid " + fill + " U","Low " + fill + " U"]
#labels2 = ["Umax Lat High","Umax Lat Mid","Umax Lat Low"]

plot = Ngl.panel(wks,toplot,[7,1],panelres)






# And now plot 
wkres2 = Ngl.Resources()
wkres2.wkColorMap = "precip_diff_12lev"
wks_type = "eps"
wks2 = Ngl.open_wks(wks_type,FigDir + figtitle + "_regrid",wkres2)

res = Ngl.Resources()
res.cnInfoLabelOn         = False    # Turn off informational
                                              # label.
res.pmLabelBarDisplayMode = "Always" # Turn on label bar.
res.cnLinesOn             = False    # Turn off contour lines.

res.nglDraw  = False
res.nglFrame = False

res.sfMissingValueV = -9999


res.cnFillOn = True
res.cnMissingValFillColor = "white"
res.cnLineLabelsOn       = False
res.pmLabelBarDisplayMode = "Always"
res.cnLinesOn =  False

res.sfXCStartV = float(OLRlon[0])
res.sfXCEndV = float(OLRlon[len(OLRlon)-1])
res.sfYCStartV = float(OLRlat[0])
res.sfYCEndV = float(OLRlat[len(OLRlat)-1])

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

res.cnLevelSelectionMode = "ManualLevels" # Define your own

toplot = []

res.cnMinLevelValF       = plotmin[1]          # contour levels.
res.cnMaxLevelValF       = plotmax[1]
res.cnLevelSpacingF      = plotspace[1]
res.tiMainString = 'TRMM' + str(daystart) + '-' + str(dayend) + 'day Precip Annual regression map'
toplot.append(Ngl.contour_map(wks2,mTPOLR,res))

res.cnMinLevelValF       = plotmin[2]          # contour levels.
res.cnMaxLevelValF       = plotmax[2]
res.cnLevelSpacingF      = plotspace[2]

res.tiMainString = 'ERA ' + str(daystart) + '-' + str(dayend) + 'day Precip Annual regression map'
toplot.append(Ngl.contour_map(wks2,mEPOLR,res))

res.cnMinLevelValF       = plotmin[3]          # contour levels. 
res.cnMaxLevelValF       = plotmax[3] 
res.cnLevelSpacingF      = plotspace[3] 
 
res.tiMainString = "OLR Annual regression map" 
toplot.append(Ngl.contour_map(wks2,mO,res)) 
 
#plot = Ngl.contour(wks,var,res) 
panelres = Ngl.Resources() 
panelres.nglPanelLabelBar = True 
#panelres.nglPanelYWhiteSpacePercent = 5. 
#panelres.nglPanelXWhiteSpacePercent = 5. 
 
panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar 
panelres.nglPanelTop                      = 0.95 
panelres.nglPanelBottom                      = 0.01 
 
#panelres.nglPanelFigureStrings            = ["a","b","c","d","e","f"] 
#panelres.nglPanelFigureStringsJust        = "BottomLeft" 
 
# 
# You can have PyNGL selection the best paper orientation for 
# the shape of plots you are drawing.  This resource is for PDF or 
# PS output only. 
# 
panelres.nglPaperOrientation = "Auto" 
 
#labels = ["High " + fill + " U","Mid " + fill + " U","Low " + fill + " U"] 
#labels2 = ["Umax Lat High","Umax Lat Mid","Umax Lat Low"] 
 
plot = Ngl.panel(wks2,toplot,[4,1],panelres) 

Ngl.end()
