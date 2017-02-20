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
import math


Data = "OLR"

Version = "Standard"
if Data == "OLR":
	plotmin1 = 100
	plotmax1 = 300
	plotspace1 = 20

	plotmin2 = -1.0
	plotmax2 = 1.0
	plotspace2 = 0.2
else:
	plotmin1 = 0
	plotmax1 = 0.5
	plotspace1 = 0.05
	
	plotmin2 = -0.01
	plotmax2 = 0.01
	plotspace2 = 0.002
	
plotmin3 = 0
plotmax3 = 1.0
plotspace3 = 0.1

plotmin4 = 0.0
plotmax4 = 0.1
plotspace4 = 0.01

MinLonF = 0
MaxLonF = 360
MinLatF = -45
MaxLatF = 45

sumlats = 10
sumlons = 20

Seas = ['DJF','MAM','JJA','SON','Ann']
nseas = 4
startyr = 1998 # Don't change - tied to file names!
endyr = 2014

FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

if Data == "TRMM":
	DirIn = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'
	FileIn = 'Regrid'+ '_' + str(sumlats) + '_' + str(sumlons) + '_TRMM_'  + str(startyr) + '-' + str(endyr) + '_seasmean_nonan.nc'
elif Data == "ERAI":
	DirIn = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
	FileIn = "SeasAnn_ERAI_Totalprecip_1980-2015_preprocess.nc"
elif Data == "CMAP":
	DirIn = "/home/disk/eos4/rachel/Obs/CMAP/"
	FileIn = "SeasAnn_precip.mon.mean.nc"
elif Data == "OLR":
        DirIn = "/home/disk/eos4/rachel/Obs/OLR/"
        FileIn = "SeasAnn_olr.mon.mean.nc" 
print DirIn + FileIn

#Get lons and lats
FileIn = xray.open_dataset(DirIn + FileIn)

lats = FileIn['Latitude']
lons = FileIn['Longitude']

nlats = len(lats)
nlons = len(lons)

if Data == "OLR":
	Density = FileIn['OLRSeas']
	DensityAnn = FileIn['OLRAnn']

else:
	Density = FileIn['PrecipSeas']
	DensityAnn = FileIn['PrecipAnn']

if Data == "TRMM":
	nyears = Density.shape[0]
	print nyears
	trendyears = 16
	startyeari = 0
	endyeari = nyears
else:
	Years = FileIn['Year']
	# take last 16 years
	nyears = len(Years)
	if nyears > 17:
		trendyears = 16
		endyeari = nyears-1
		startyeari = nyears-17
	else:
		trendyears = nyears-1
		endyeari = nyears-1
		startyeari = 0
	startyr = Years[startyeari].values
	endyr = Years[endyeari].values

print startyeari,endyeari
figtitle = Data + '_Seas_' + str(startyr) + '-' + str(endyr)
figtitle2 = Data + '_Seas_r2_p_' + str(startyr) + '-' + str(endyr)


if Data == "CMAP":
        #Convert to mm/hr from mm/day
        Density = Density[startyeari:endyeari,:,:,:]/24.0
        DensityAnn = DensityAnn[startyeari:endyeari,:,:]/24.0
else:
        Density = Density[startyeari:endyeari,:,:,:]
        DensityAnn = DensityAnn[startyeari:endyeari,:,:]
	
#convert nans to missing values so plots work!
Density._FillValue = -999.0
DensityAnn._FillValue = -999.0
DensityAnn = np.where(np.isnan(DensityAnn[:,:,:]),DensityAnn._FillValue,DensityAnn[:,:,:])
Density = np.where(np.isnan(Density[:,:,:,:]),Density._FillValue,Density[:,:,:,:])

# Calculate linear trend over the 16 years 
# Calculate linear trend over the 16 years 
yearnums = range(0,trendyears)
A = np.array(yearnums)
#A = np.vstack([yearnums,np.ones(len(yearnums))]).T
c = np.zeros((nseas+1,nlats,nlons),np.float)
m = np.zeros((nseas+1,nlats,nlons),np.float)
r = np.zeros((nseas+1,nlats,nlons),np.float)
p = np.zeros((nseas+1,nlats,nlons),np.float)
stderr = np.zeros((nseas+1,nlats,nlons),np.float)

cAnn = np.zeros((nlats,nlons),np.float)
mAnn = np.zeros((nlats,nlons),np.float)
rAnn = np.zeros((nlats,nlons),np.float)
pAnn = np.zeros((nlats,nlons),np.float)
stderrAnn = np.zeros((nlats,nlons),np.float)

for ilat in range(0,nlats):
        for ilon in range(0,nlons):
                for iseas in range(0,nseas):
                        linreg = Density[:,iseas,ilat,ilon]
                        m[iseas,ilat,ilon],c[iseas,ilat,ilon],r[iseas,ilat,ilon],p[iseas,ilat,ilon],stderr[iseas,ilat,ilon] = stats.linregress(A,linreg)


		# Calculate annual average
		iseas = nseas
		linreg = DensityAnn[:,ilat,ilon]
		m[iseas,ilat,ilon],c[iseas,ilat,ilon],r[iseas,ilat,ilon],p[iseas,ilat,ilon],stderr[iseas,ilat,ilon] = stats.linregress(A,linreg)

print m[0,:,:]

# Calculate r2
r2 = r * r

# Read in land-sea mask? How do we get one for TRMM? Generic one, regridded from other, with min 50% land?

# Define lat-lon mask for:
# Pacific: 120E - 280E
# Atlantic: 280E - 360E and 0 - 20E 
# Indian: 35E - 120E
# Calculate basin averages
Psrt = 120.0
if lons[0] < 0:
        Pend = -80.0      #280.0
else:
        Pend = 280.0

if lons[0] < 0:
        Asrt = -80.0    #280
else:
        Asrt = 280.0
Aend = 20.0

Isrt = 35.0
Iend = 120.0


Latstart = 50.0
Latend = -50.0
# Find start lons indices for these
diffs = abs(lons[1] - lons[0])/2.0
diffslat = abs(lats[1] - lats[0])/2.0
print Psrt,Pend,Asrt,Aend,Isrt,Iend

Startlati = np.intersect1d(np.where(lats >= Latstart - diffs),np.where(lats < Latstart + diffs),False)[0]
Endlati = np.intersect1d(np.where(lats >= Latend - diffs),np.where(lats < Latend + diffs),False)[0]

print Startlati,Endlati
def SplitBasin(srt,end,inputdata):
        Starti = np.intersect1d(np.where(lons >= srt - diffs),np.where(lons < srt + diffs),False)[0]
        Endi = np.intersect1d(np.where(lons > end - diffs),np.where(lons <= end + diffs),False)[0]

        if Starti > Endi:
                outden = np.concatenate([inputdata[...,Startlati:Endlati,Starti:nlons],inputdata[...,Startlati:Endlati,0:Endi]],axis = 3)
        else:
                outden = inputdata[...,Startlati:Endlati,Starti:Endi]

        return outden


PacDen = SplitBasin(Psrt,Pend,Density)
AtlDen = SplitBasin(Asrt,Aend,Density)
IndDen = SplitBasin(Isrt,Iend,Density)

AllDen = np.sum(Density[...,Startlati:Endlati,:],axis = (2,3))
AtlDen = np.sum(AtlDen,axis = (2,3))
PacDen = np.sum(PacDen,axis = (2,3))
IndDen = np.sum(IndDen,axis = (2,3))

# And now plot 
wkres = Ngl.Resources()
wkres.wkColorMap = "precip_diff_12lev"
wks_type = "eps"
wks = Ngl.open_wks(wks_type,FigDir + figtitle,wkres)
wks2 = Ngl.open_wks(wks_type,FigDir + figtitle2,wkres)

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

res.cnLevelSelectionMode = "ManualLevels" # Define your own

toplot = []
toplot2 = []

for iseas in range(0,nseas+1):
	res.cnMinLevelValF       = plotmin1          # contour levels.
	res.cnMaxLevelValF       = plotmax1
	res.cnLevelSpacingF      = plotspace1
	res.tiMainString = Seas[iseas]
	print 'seasonal means'
	if iseas < nseas:
		toplot.append(Ngl.contour_map(wks,np.nanmean(Density[:,iseas,:,:],axis=0),res))
	else:
                toplot.append(Ngl.contour_map(wks,np.nanmean(DensityAnn[:,:,:],axis=0),res))

        res.cnMinLevelValF       = plotmin2          # contour levels.
        res.cnMaxLevelValF       = plotmax2
        res.cnLevelSpacingF      = plotspace2
        res.tiMainString = Seas[iseas] + " linear regression slope"
        toplot.append(Ngl.contour_map(wks,m[iseas,:,:],res))

for iseas in range(0,nseas+1):
        res.cnMinLevelValF       = plotmin3          # contour levels.
        res.cnMaxLevelValF       = plotmax3
        res.cnLevelSpacingF      = plotspace3
        res.tiMainString = Seas[iseas] + "r^2"
        toplot2.append(Ngl.contour_map(wks2,r2[iseas,:,:],res))

        res.cnMinLevelValF       = plotmin4        # contour levels.
        res.cnMaxLevelValF       = plotmax4
        res.cnLevelSpacingF      = plotspace4
        res.tiMainString = Seas[iseas] + " pvalue"
        toplot2.append(Ngl.contour_map(wks2,p[iseas,:,:],res))



textres = Ngl.Resources()
textres.txFontHeightF = 0.015
Ngl.text_ndc(wks,"Average seasonal precip trend analysis",0.5,0.87,textres)
Ngl.text_ndc(wks2,"T test statistics for seasonal precip trend analysis",0.5,0.87,textres)


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

plot = Ngl.panel(wks,toplot,[nseas+1,2],panelres)
plot = Ngl.panel(wks2,toplot2,[nseas+1,2],panelres)

def plotdensity(indensity,title):
        figtitle = Data + "_" + title + '_' + str(startyr) + '-' + str(endyr)

        regresst = (np.zeros((5,nseas),np.float))

        for iseas in range(0,nseas):
                linreg = indensity[:,iseas]
                regresst[:,iseas] = stats.linregress(A,linreg)

        #Now plot 
        wkres = Ngl.Resources()
        wkres.wkColorMap = "MPL_BrBG"
        wks_type = "eps"
        wks = Ngl.open_wks(wks_type,FigDir + figtitle,wkres)
        res = Ngl.Resources()
        res.nglFrame = False
        res.nglDraw = False
        res.xyMarkLineMode = "Markers"
        res.xyMonoMarkLineMode = True
        res.xyMarkerColor = "red"
        res.xyMarkers = 16
        res.xyMarkerSizeF = 0.01
        res.xyYStyle = "Linear"
        res.xyXStyle = "Linear"

        res.tiYAxisString = "year"

        plot = []

	if Data == "OLR":
	        res.tiYAxisString = "OLR"
	else:
		res.tiYAxisString = "Precip"

        #print type(indensity[:,ibin])
        #print type(range(startyr,endyr+1))

        for iseas in range(0,nseas):
                res.tiMainString = Seas[iseas] + '; r2= ' + '{:5.3f}'.format(regresst[2,iseas] * regresst[2,iseas]) + '; p= ' + '{:5.3f}'.format(regresst[3,iseas])

                plot.append(Ngl.xy(wks,range(startyr,endyr+1),indensity[:,iseas],res))

        panelres = Ngl.Resources()
        panelres.nglPanelLabelBar = True
        panelres.nglPanelYWhiteSpacePercent = 5.
        panelres.nglPanelXWhiteSpacePercent = 5.

        panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
        panelres.nglPanelTop                      = 0.98
        panelres.nglPanelBottom                      = 0.02
        panelres.nglPaperOrientation = "Auto"

        txres = Ngl.Resources()
        txres.txFontHeightF = 0.015
        Ngl.text_ndc(wks,'Annual and seasonal timeseries ' + Data + " for " + title ,0.5,0.85,txres)
        
        Ngl.panel(wks,plot,[math.ceil(float(nseas)/3.0),3],panelres)

 
plotdensity(AtlDen,"Atlantic")
print "Pacific"
print type(PacDen)
plotdensity(PacDen,"Pacific")
print "Indian"
plotdensity(IndDen,"Indian")
print "All"
plotdensity(AllDen,"All")



Ngl.end()




