
# coding: utf-8


# In[2]:

import os, errno
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xarray as xr
#import Ngl
#import math
from scipy import stats
from scipy.interpolate import griddata

from rhwhitepackages.readwrite import *
from rhwhitepackages.stats import regressmaps
from rhwhitepackages.plots import *

# plotting
# import matplotlib
import xray.plot as xplt
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plotboth = True
pthresh = 0.05

FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

MinLatF = -40
MaxLatF = 40
MinLonF = 0
MaxLonF = 360

minLats = [  2,  8, 22,-40]
maxLats = [  8, 20, 35,-25]
minLons = [200,155,150, 45]
maxLons = [280,200,180, 75]

nboxes = len(minLats)
# In[3]:

#Select lats, lons, bounds, and other parameters
anstartyr = 1999
anendyr = 2014
splittype = 'day'
tbound1 = [0.0,1.0,2.0,5.0]
tbound2 = [1.0,2.0,5.0,100.0]

#tbound1 = [1.0, 1.125, 1.5, 1.875]
#tbound2 = [1.125, 1.25, 1.625, 2.0]

nbounds = len(tbound1)
if len(tbound2) != nbounds:
    exit("tbound arrays not equal lengths")
unit = 'day'

if splittype == 'day':
    diradd = '/Sizes/'
elif splittype == 'MaxSpeeds':
    diradd = '/MaxSpeeds/'

sumlats=0
sumlons=0

minGB=0
mapping = 'center'

#minLat = [-45,-10,8,20,15]; maxLat = [-25,10,15,35,30]; minLon = [185,120,150,140,290]; maxLon = [220,160,220,220,340]  #225
minLat = [  8,  2, 22,-40]
maxLat = [  20,  8, 35,-25]
minLon = [155,200,150, 45]
maxLon = [200,280,180, 75]


#Read in files
fstartyr = 1998
fendyr = 2014
precipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"
precipClimFile = "TRMM_3B42_1998-2014_annmean.nc"
precipRegFile = 'regress_TRMM_3B42_1999-2014_annmean.nc'
precipReg2File = 'regress_TRMM_3B42_1998-2014_annmean.nc'

olrDir = '/home/disk/eos4/rachel/Obs/OLR/'
olrRegFile = 'regress_olr_1999-2013.nc'

filePrecip = xrayOpen(precipClimDir + precipClimFile)
filePrecipReg =  xrayOpen(precipClimDir + precipRegFile)
filePrecipReg2 =  xrayOpen(precipClimDir + precipReg2File)

fileOlrReg = xrayOpen(olrDir + olrRegFile)

# Open events files
iday = 1
TRMMden = getdenfilename(mapping, 'TRMM', 'Standard', 1998, 2014, iday,
                 splittype, unit, 4, minGB, tbound1, tbound2,
                 'Mon',sumlats,sumlons)
TRMM1 = geteventmapdata(iday,'Ann','TDensity',TRMMden,
                    1998, 2014,
                    MinLatF,MaxLatF,MinLonF,MaxLonF,
                    ['lat','lon'])
# 2-5 day data
iday = 2

TRMMden = getdenfilename(mapping, 'TRMM', 'Standard', 1998, 2014, iday,
                 splittype, unit, 4, minGB, tbound1, tbound2,
                 'Mon',sumlats,sumlons)
TRMM2 = geteventmapdata(iday,'Ann','TDensity',TRMMden,
                    1998, 2014,
                    MinLatF,MaxLatF,MinLonF,MaxLonF,
                    ['lat','lon'])

# ERA 1-2 day data
iday = 1
ERAden = getdenfilename(mapping, 'ERAI', 'Standard', 1980, 2014, iday,
                 splittype, unit, 4, minGB, tbound1, tbound2,
                 'Mon',sumlats,sumlons)
ERA1 = geteventmapdata(iday,'Ann','TDensity',ERAden,
                    1980, 2014,
                    MinLatF,MaxLatF,MinLonF,MaxLonF,
                    ['lat','lon'])

# ERA 2-5 day data
iday = 2
ERAden = getdenfilename(mapping, 'ERAI', 'Standard', 1980, 2014, iday,
                 splittype, unit, 4, minGB, tbound1, tbound2,
                 'Mon',sumlats,sumlons)
ERA2 = geteventmapdata(iday,'Ann','TDensity',ERAden,
                    1980, 2014,
                    MinLatF,MaxLatF,MinLonF,MaxLonF,
                    ['lat','lon'])

# Get precip data
preciplons = filePrecipReg['pcp'].coords['longitude']
preciplats = filePrecipReg['pcp'].coords['latitude']

dataOLR = fileOlrReg['olr']
lonsOLR = fileOlrReg['olr'].coords['lon']
latsOLR = fileOlrReg['olr'].coords['lat']

# need to make these arrays

if np.amax(preciplons) < 190:
    olrtemp1 = preciplons.sel(longitude=slice(np.amin(preciplons),0)) + 360.0
    olrtemp2 = preciplons.sel(longitude=slice(0,np.amax(preciplons)))
    olrnewlons = xr.concat((olrtemp2,olrtemp1),dim='longitude')
else:
    olrnewlons = preciplons

print olrnewlons

OLRnew = dataOLR.sel(lat=preciplats,method='nearest')
OLRnew2 = OLRnew.sel(lon=olrnewlons,method='nearest')

# Now mask out filePrecipReg['pcp'] based on where the sign is not opposite to
# OLRnew2

precipclimin = filePrecipReg['pcp']
origlonsin = filePrecipReg['pcp'].coords['longitude']
if np.amax(origlonsin) < 190:

    newprecip2 = precipclimin.sel(longitude=slice(np.amin(origlonsin),0))
    newprecip2.coords['longitude'] = newprecip2.coords['longitude'] + 360
    precipclim = xr.concat((
            precipclimin.sel(longitude=slice(0,np.amax(origlonsin)+1)),
                            newprecip2),
                            dim='longitude')

signs = OLRnew2[0,:,:].values * precipclim[0,:,:].values



OLRprecip = precipclim.copy(deep=True)
# Change units if necessary
if OLRprecip.long_name == 'precipitation (mm/hr)':
    OLRprecip.values = OLRprecip.values * 24.0 * 365.0
    OLRprecip['long_name'] = 'precipitation (mm/yr/yr)'
    OLRprecip['units'] = 'mm/yr/yr'
else:
    exit('i don\'t know the units, and I\'m not going to guess, sorry!')



for ilon in range(0,len(OLRprecip.longitude)):
    for ilat in range(0,len(OLRprecip.latitude)):
        if signs[ilat,ilon] > 0: 
            OLRprecip[0,ilat,ilon] = -9999  #set equal to missing data value

OLRprecip = OLRprecip.rename('mm/yr/yr')

plot = []
print OLRprecip 
figtitle = FigDir + 'paper_OLR_Precip_Regressions_1-2day_TRMM_ERA'

wkres = Ngl.Resources()
wkres.wkColorMap = "BlueDarkRed18"
wks_type = "eps"
wks = Ngl.open_wks(wks_type,figtitle,wkres)

res1 = Ngl.Resources()
initcontourplot(res1,MinLatF,MinLonF,MaxLatF,MaxLonF,OLRprecip.latitude.values,OLRprecip.longitude.values)
res1.nglMaximize = False
res1.nglDraw     = False
res1.nglFrame    = False
res1.sfMissingValueV = -9999
# play with labelbar
res1.lbOrientation   = "Vertical"
res1.pmLabelBarWidthF = 0.05
res1.lbBoxEndCapStyle = "TriangleBothEnds"
res1.lbBoxMajorExtentF = 1.0
#res1.lbLabelAutoStride = True
res1.lbAutoManage = False
# including some font heights
res1.lbLabelFontHeightF = 0.007
res1.lbTitleFontHeightF = 0.007
res1.lbTitlePosition = "bottom"
res1.lbTitleExtentF = 0.05
res1.lbTopMarginF = 0.01
res1.lbBottomMarginF = 0.01
res1.tiMainFontHeightF = 0.015

# set position
res1.vpWidthF  = 0.75
res1.vpHeightF = 0.22
res1.vpXF      = .1
res1.vpYF      = .95

res1.cnMinLevelValF       = -35
res1.cnMaxLevelValF       = 35
res1.cnLevelSpacingF = 10

# turn off grideines
res1.mpGridAndLimbOn = False

print OLRprecip.values.shape

# add label bar title for units
res1.lbTitleOn = True
res1.lbTitleString = "mm/yr~S~2"

res1.tiMainPosition = 'Left'
res1.tiMainOffsetXF = -0.03
res1.tiMainString = 'a.'
plot = (Ngl.contour_map(wks,OLRprecip[0,:,:].values,res1))

# Add boxes
gsres = Ngl.Resources()
gsres.gsLineColor       = "Black"
txres               = Ngl.Resources()
txres.txFont        = "helvetica-bold"
txres.txFontHeightF = 0.012
txres.txJust        = "BottomCenter"

for i in range(0,nboxes):
    Ngl.add_polyline(wks,plot,[minLons[i],minLons[i]],[minLats[i],maxLats[i]],gsres)
    Ngl.add_polyline(wks,plot,[maxLons[i],maxLons[i]],[minLats[i],maxLats[i]],gsres)
    Ngl.add_polyline(wks,plot,[minLons[i],maxLons[i]],[minLats[i],minLats[i]],gsres)
    Ngl.add_polyline(wks,plot,[minLons[i],maxLons[i]],[maxLats[i],maxLats[i]],gsres)

    Ngl.add_text(wks,plot,str(i+1),0.5*(minLons[i] + maxLons[i]),minLats[i] +
                    0.25 * (maxLats[i] - minLats[i]) ,txres)

# Add label 
    Ngl.add_text(wks,plot,str(i+1),0.5*(minLons[i] + maxLons[i]),minLats[i] +
                    0.25 * (maxLats[i] - minLats[i]) ,txres)


Ngl.draw(plot)

# Add other plots
# TRMM line plot
#Now plot plots
res2 = Ngl.Resources()
res2.nglMaximize = False
res2.nglFrame = False
res2.nglDraw = False
res2.xyMarkLineMode = "Markers"
res2.xyMonoMarkLineMode = True
res2.xyMarkerColor = "blue"
res2.xyMarkers = 5
res2.xyMarkerSizeF = 0.005
res2.xyYStyle = "Linear"
res2.xyXStyle = "Linear"

#res2.tiMainFontHeightF = 0.035
#res2.tiYAxisFontHeightF = 0.03
res2.tiXAxisOn = False
#res2.tmXBLabelFontHeightF = 0.03
#res2.tmYLLabelFontHeightF = 0.03
#res2.tmYLMode = "Automatic"
#res2.tmYLFormat = "@6^g"
res2.tiYAxisString      = "TRMM"
res2.tiYAxisFontColor   = "blue"
# Turn off right hand axis labels
res2.tmYROn       = False
res2.tmYRLabelsOn = False
res2.tmYRMinorOn  = False

res2.vpWidthF = 0.3
res2.vpHeightF = 0.3
res2.vpXF      = .1
res2.vpYF      = .65

print TRMM1.year[0].values
res2.trXMinF     = ERA1.year[0].values - 1
res2.trXMaxF     = ERA1.year[-1].values + 1

res2.tiMainPosition = 'Left'
res2.tiMainOffsetXF = -0.03
res2.tiMainString = 'b.'

plot = (Ngl.xy(wks,TRMM1.year.values,TRMM1.values,res2))

Ngl.add_text(wks,plot,'b',ERA1.year.values[0],Ngl.get_float(plot,"tmYLMajorLengthF"),txres)

# Add 2-5 day data

res22 = Ngl.Resources()
res22.nglDraw = False
res22.nglFrame = False
res22.nglMaximize = False
res22.vpHeightF = res2.vpHeightF
res22.vpWidthF = res2.vpWidthF
res22.vpXF = res2.vpXF
res22.vpYF = res2.vpYF
res22.trXMinF = res2.trXMinF
res22.trXMaxF = res2.trXMaxF

#
# Turn off the bottom and left tickmarks, since we will use the ones from
# the first plot.
#
res22.tmXBOn       = False
res22.tmXBLabelsOn = False
res22.tmXBMinorOn  = False
res22.tmYLOn       = False
res22.tmYLLabelsOn = False
res22.tmYLMinorOn  = False

#
# Turn on the right Y labels and tickmarks and move the axis string to
# the right.
#
res22.tmYRLabelsOn = True
res22.tmYROn       = True
res22.tmYUseLeft   = False
res22.tmYRFormat   = "f"      # Gets rid of unnecessary trailing zeros

#
# Move the Y axis string to the right.
#
res22.tiYAxisString      = "ERA-Interim"
res22.tiYAxisSide        = "Right"
res22.tiYAxisFontColor   = "red"
res22.tiXAxisFontHeightF = Ngl.get_float(plot,"tiXAxisFontHeightF")
res22.xyMarkerColor = "red"
res22.xyMarkers = '.'
res22.xyMarkLineMode = "Markers"
res22.xyMonoMarkLineMode = True
#
# Make sure the font heights and tickmark lengths are the same as
# the first plot.
#
res22.tmYRLabelFontHeightF = Ngl.get_float(plot,"tmYLLabelFontHeightF")
res22.tmYRMajorLengthF     = Ngl.get_float(plot,"tmYLMajorLengthF")

Ngl.draw(plot)

if plotboth:    # if want to include 1-2 AND 2-5 day events
    plot2 = (Ngl.xy(wks,ERA1.year.values,ERA1.values,res22))
    Ngl.draw(plot2)

# Repeat for ERA-I data next to first plot
res2.vpXF      = .6
res2.vpYF      = .65
res22.vpXF      = .6
res22.vpYF      = .65


res2.trXMinF     = ERA1.year[0].values - 1
res2.trXMaxF     = ERA1.year[-1].values + 1
res22.trXMinF     = res2.trXMinF
res22.trXMaxF     = res2.trXMaxF

res2.tiMainString = 'c.'
plot = (Ngl.xy(wks,TRMM2.year.values,TRMM2.values,res2))

res22.tmYRLabelFontHeightF = Ngl.get_float(plot,"tmYLLabelFontHeightF")
res22.tmYRMajorLengthF     = Ngl.get_float(plot,"tmYLMajorLengthF")

A = range(1980,2001+1)
linreg = ERA2.sel(year=slice(1980,2001)).values
print A
print linreg.shape
regress = stats.linregress(A,linreg)
print regress

Ngl.draw(plot)

if plotboth:
    plot2 = (Ngl.xy(wks,ERA2.year.values,ERA2.values,res22))
    Ngl.draw(plot2)


# Maximize and draw
#psres = True
#Ngl.maximize_plot(wks,plot,psres) 
Ngl.frame(wks)

quit()

