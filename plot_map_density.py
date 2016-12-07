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

size = 1
Data = "TRMM"

Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

if size == 0:
	plotmin1 = 0
	plotmax1 = 1000
	plotspace1 = 100

	plotmin2 = -50
	plotmax2 = 50
	plotspace2 = 10

        plotmin3 = 0
        plotmax3 = 1000
        plotspace3 = 100

        plotmin4 = -50
        plotmax4 = 50
        plotspace4 = 10
elif size == 1:
        plotmin1 = 0
        plotmax1 = 10
        plotspace1 = 1

        plotmin2 = -0.5
        plotmax2 = 0.5
        plotspace2 = 0.1

        plotmin3 = 0
        plotmax3 = 10
        plotspace3 = 1

        plotmin4 = -0.5
        plotmax4 = 0.5
        plotspace4 = 0.1

elif size == 2:
	plotmin1 = 0
	plotmax1 = 2
	plotspace1 = 0.2

	plotmin2 = -0.1
	plotmax2 = 0.1
	plotspace2 = 0.02

	plotmin3 = 0
	plotmax3 = 2
	plotspace3 = 0.2

	plotmin4 = -0.1
	plotmax4 = 0.1
	plotspace4 = 0.02


MinLonF = 0
MaxLonF = 360
MinLatF = -45
MaxLatF = 45

sumlats = 16
sumlons = 16

Seas = ['DJF','MAM','JJA','SON']
nseas = 4

# Time period for analysis
astartyr = 1998
aendyr = 2014

FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2014
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

elif Data == "ERAI":
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
elif Data == "CESM":
	startyr = 1990 # Don't change - tied to file names!
	endyr = 2014
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'


istartyr = astartyr - startyr
iendyr = aendyr - startyr
nyears = iendyr - istartyr

FileIn = 'DenDirSpd_Map_monthly_regrid_' + Data + '_'  + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + Version + '.nc'

figtitle = Data + '_2SeasDensity_' + Version + '_' + str(astartyr) + '-' + str(aendyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats' + '_' + str(size)
figtitle2 = Data + '_2SeasWestwardDensity_' + Version + '_' + str(astartyr) + '-' + str(aendyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats' + str(size)
figtitle3 = Data + '_2SeasEastwardDensity_' + Version + '_' + str(astartyr) + '-' + str(aendyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats' + str(size)
figtitle4 = Data + '_2SeasDensityPerc_' + Version + '_' + str(astartyr) + '-' + str(aendyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats' + str(size)


#Get lons and lats
print DirIn + FileIn
FileIn = xray.open_dataset(DirIn + FileIn)

lats = FileIn['Latitude']
lons = FileIn['Longitude']

nlats = len(lats)
nlons = len(lons)

Density = FileIn['DensityMapSeas'][istartyr:iendyr,...]
WWDensity = FileIn['WestwardDensityMapSeas'][istartyr:iendyr,...]
EWDensity = FileIn['EastwardDensityMapSeas'][istartyr:iendyr,...]

DensityAnn = FileIn['DensityMapAnn'][istartyr:iendyr,...]
Direction = FileIn['DirectionMapSeas'][istartyr:iendyr,...]
Speed = FileIn['SpeedMapSeas'][istartyr:iendyr,...]

print Density.shape

# Calculate linear trend over the 16 years 
yearnums = range(0,nyears)
A = np.array(yearnums)

print A.shape
#A = np.vstack([yearnums,np.ones(len(yearnums))]).T
c = np.zeros((3,nseas,nlats,nlons),np.float)
m = np.zeros((3,nseas,nlats,nlons),np.float)
r = np.zeros((3,nseas,nlats,nlons),np.float)
p = np.zeros((3,nseas,nlats,nlons),np.float)
stderr = np.zeros((3,nseas,nlats,nlons),np.float)

cAnn = np.zeros((nlats,nlons),np.float)
mAnn = np.zeros((nlats,nlons),np.float)
rAnn = np.zeros((nlats,nlons),np.float)
pAnn = np.zeros((nlats,nlons),np.float)
stderrAnn = np.zeros((nlats,nlons),np.float)

for ilat in range(0,nlats):
	for ilon in range(0,nlons):
		for iseas in range(0,nseas):
			linreg = Density[:,iseas,size,ilat,ilon]
			m[0,iseas,ilat,ilon],c[0,iseas,ilat,ilon],r[0,iseas,ilat,ilon],p[0,iseas,ilat,ilon],stderr[0,iseas,ilat,ilon] = stats.linregress(A,linreg)
			WWlinreg = WWDensity[:,iseas,size,ilat,ilon]
                        m[1,iseas,ilat,ilon],c[1,iseas,ilat,ilon],r[1,iseas,ilat,ilon],p[1,iseas,ilat,ilon],stderr[1,iseas,ilat,ilon] = stats.linregress(A,WWlinreg)
                        EWlinreg = EWDensity[:,iseas,size,ilat,ilon]
                        m[2,iseas,ilat,ilon],c[2,iseas,ilat,ilon],r[2,iseas,ilat,ilon],p[2,iseas,ilat,ilon],stderr[2,iseas,ilat,ilon] = stats.linregress(A,EWlinreg)


# Calculate r2
r2 = r * r

for ilat in range(0,nlats):
        for ilon in range(0,nlons):
                for iseas in range(0,nseas):
                        linreg = DensityAnn[:,size,ilat,ilon]
                        mAnn[ilat,ilon],cAnn[ilat,ilon],rAnn[ilat,ilon],pAnn[ilat,ilon],stderrAnn[ilat,ilon] = stats.linregress(A,linreg)


# And now plot 
wkres = Ngl.Resources()
wkres.wkColorMap = "precip_diff_12lev"
wks_type = "eps"
wks = Ngl.open_wks(wks_type,FigDir + figtitle,wkres)
wks2 = Ngl.open_wks(wks_type,FigDir + figtitle2,wkres)
wks3 = Ngl.open_wks(wks_type,FigDir + figtitle3,wkres)
wks4 = Ngl.open_wks(wks_type,FigDir + figtitle4,wkres)


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
toplot3 = []
toplot4 = []

for iseas in range(0,nseas):
	res.cnMinLevelValF       = plotmin1          # contour levels.
	res.cnMaxLevelValF       = plotmax1
	res.cnLevelSpacingF      = plotspace1
	res.tiMainString = Seas[iseas] + " events/year"
	toplot.append(Ngl.contour_map(wks,np.mean(Density[:,iseas,size,:,:],axis=0),res))
	
        res.cnMinLevelValF       = plotmin2          # contour levels.
        res.cnMaxLevelValF       = plotmax2
        res.cnLevelSpacingF      = plotspace2
        res.tiMainString = Seas[iseas] + " linear regression slope, events/year"
        toplot.append(Ngl.contour_map(wks,m[0,iseas,:,:],res))
for iseas in range(0,nseas):
        res.cnMinLevelValF       = plotmin3          # contour levels.
        res.cnMaxLevelValF       = plotmax3
        res.cnLevelSpacingF      = plotspace3
        res.tiMainString = Seas[iseas] + " Westward, events/year"
        toplot2.append(Ngl.contour_map(wks2,np.mean(WWDensity[:,iseas,size,:,:],axis=0),res))
        
	res.cnMinLevelValF       = plotmin4        # contour levels.
        res.cnMaxLevelValF       = plotmax4
        res.cnLevelSpacingF      = plotspace4
        res.tiMainString = Seas[iseas] + " Westward linear regression slope, events/year"
        toplot2.append(Ngl.contour_map(wks2,m[1,iseas,:,:],res))

for iseas in range(0,nseas):
        res.cnMinLevelValF       = plotmin3          # contour levels.
        res.cnMaxLevelValF       = plotmax3
        res.cnLevelSpacingF      = plotspace3
        res.tiMainString = Seas[iseas] + " Eastward"
        toplot3.append(Ngl.contour_map(wks3,np.mean(EWDensity[:,iseas,size,:,:],axis=0),res))

        res.cnMinLevelValF       = plotmin4        # contour levels.
        res.cnMaxLevelValF       = plotmax4
        res.cnLevelSpacingF      = plotspace4
        res.tiMainString = Seas[iseas] + " Eastward linear regression slope"
        toplot3.append(Ngl.contour_map(wks3,m[2,iseas,:,:],res))

for iseas in range(0,nseas):
        res.cnMinLevelValF       = plotmin1          # contour levels.
        res.cnMaxLevelValF       = plotmax1
        res.cnLevelSpacingF      = plotspace1
        res.tiMainString = Seas[iseas] + " events/year"
        toplot4.append(Ngl.contour_map(wks4,np.mean(Density[:,iseas,size,:,:],axis=0),res))

        res.cnMinLevelValF       = -20.0          # contour levels.
        res.cnMaxLevelValF       = 20.0
        res.cnLevelSpacingF      = 4.0
        res.tiMainString = Seas[iseas] + " regression fraction, %/year"
	globseasmean = np.mean(Density[:,iseas,size,:,:])
	fraction =  m[0,iseas,:,:] * 100.0/np.float(globseasmean)
        toplot4.append(Ngl.contour_map(wks4,fraction,res))


textres = Ngl.Resources()
textres.txFontHeightF = 0.015
if size == 0:
	spantitle = "0-1"
elif size == 1:
	spantitle = "1-2"
elif size == 2:
	spantitle = "2-6"
Ngl.text_ndc(wks,"Number of events of length " + spantitle + " days, " + Version + ' analysis',0.5,0.87,textres)
Ngl.text_ndc(wks2,"Number of westward moving events of length " + spantitle + " days, " + Version + ' analysis',0.5,0.87,textres)
Ngl.text_ndc(wks3,"Number of eastward moving events of length " + spantitle + " days, " + Version + ' analysis',0.5,0.87,textres)
Ngl.text_ndc(wks4,"Number of events of length " + spantitle + " days, " + Version + ' analysis',0.5,0.87,textres)


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

plot = Ngl.panel(wks,toplot,[4,2],panelres)
plot = Ngl.panel(wks2,toplot2,[4,2],panelres)
plot = Ngl.panel(wks3,toplot3,[4,2],panelres)
plot = Ngl.panel(wks4,toplot4,[4,2],panelres)


Ngl.end()




