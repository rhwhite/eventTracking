# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os, errno
import numpy as np
import numpy.ma as ma
import netCDF4
from netCDF4 import Dataset
import datetime as dt
import re
import sys
import Ngl
import xray
import math


startyr = 1998 # Don't change - tied to file names!
endyr = 2014

anstartyr = 1998 #year for analysis start
anendyr = 2014 #year for analysis end
nyears = anendyr - anstartyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

starttsteps = 0
endtsteps = 46752 # 16 years
anntsteps = 2920 # timesteps per year

minevent = 100000

DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'

DirP = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/'

FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '_5th_nanto25.nc'

FileP = 'TRMM_1998-2014_3B42_3hrly_nonan.nc'

FileO = 'DensityMap_' + str(startyr) + '-' + str(endyr) + '_5th_nanto25.nc'

filetimespan = "3hrly"

tempfile = xray.open_dataset(DirP + FileP)
lons = tempfile['longitude'].values
lats = tempfile['latitude'].values

nlons = len(lons)
nlats = len(lats)


#Create new datasets
RDenMap = np.zeros((nyears,4,nlats,nlons),np.int)

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('years', nyears)
ncfile.createDimension('size', 4)
ncfile.createDimension('lon', nlons)
ncfile.createDimension('lat', nlats)


DensityMap = ncfile.createVariable('DensityMap','f4',('years','size','lat','lon'),fill_value=-9999)
Longitude = ncfile.createVariable('Longitude','f4',('lon'),fill_value=-9999)
Latitude = ncfile.createVariable('Latitude','f4',('lat'),fill_value=-9999)
Years = ncfile.createVariable('years','f4',('years'),fill_value=-9999)

setattr(DensityMap,'Extra Info','Size based on timespan: 0 is < 8 (1 day), 1: < 16 (2 days) < 48 (6 days), 2: > 48)')
 
Longitude[:] = lons
Latitude[:] = lats


datain = xray.open_dataset(DirI + FileI1)

print(datain.coords)

tspan=datain['timespan'].values
tmean = datain['tmean']
tstart = datain['tstart'].values
xmin1 = datain['xmin']

xcenter = datain['xcentermean'].values
ycenter = datain['ycentermean'].values

nevents = tspan.shape[0]
print nevents

ievent = 1

curtime = tstart[starttsteps]
n = 0
iyear = 0
for n in range(0,nevents):
	curtime = tstart[n]
	if curtime == starttsteps:
		mints[iyear] = n
		break

for iyear in range(0,nyears):
	Years[iyear] = anstartyr + iyear
	for n in range(int(mints[iyear]),nevents):
		curtime = tstart[n]
		if curtime == (iyear + 1) * anntsteps:
			if iyear < nyears-1:
				mints[iyear+1] = n
			
			maxts[iyear] = n  #switched to n instead of n-1 because python doesn't use the last index in a range
			break

print mints
print maxts

for iyear in range(0,nyears):
	print mints[iyear]
	print maxts[iyear]
	for ievent in range(int(mints[iyear]),int(maxts[iyear])):

		if (ievent % 100000 == 0):
			print "ievent: " + str(ievent)	
		if tspan[ievent] <= 8:
			RDenMap[iyear,0,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1

		elif tspan[ievent] < 16:
                        RDenMap[iyear,1,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
                elif tspan[ievent] < 48:
                        RDenMap[iyear,2,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1
		else:
                        RDenMap[iyear,3,int(round(ycenter[ievent])),int(round(xcenter[ievent]))] += 1

DensityMap[:,:,:,:] = RDenMap

datain.close()






"""
#Now plot log-log plots
wkres = Ngl.Resources()
wkres.wkColorMap = "MPL_BrBG"
wks_type = "eps"
wks = Ngl.open_wks(wks_type,FigDir + figtitle1,wkres)
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

	
if (plotdensity):
	res.tiYAxisString = "frequency"
else:
	res.tiYAxisString = "number"

for ibin in range(0,histnbins):
	res.tiMainString = 'Tspan annual values for bin '  + str(tspan_bins_plot[ibin])

	plot.append(Ngl.xy(wks,range(startyr,endyr+1),tspanplot[:,ibin],res))

panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.98
panelres.nglPanelBottom                      = 0.02
panelres.nglPaperOrientation = "Auto"

print 'about to panel'


Ngl.panel(wks,plot,[math.ceil(float(histnbins)/3.0),3],panelres)

print 'panelled'

#Now plot log-log plots
wkres = Ngl.Resources()
wks = Ngl.open_wks(wks_type,FigDir + figtitle2,wkres)
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

if (plotdensity):
        res.tiYAxisString = "frequency"
else:
        res.tiYAxisString = "number"

for ibin in range(0,histnbins):
	res.tiMainString = 'Sspan annual values for bin '  + str(sspan_bins_plot[ibin])
        plot.append(Ngl.xy(wks,range(startyr,endyr+1),sspanplot[:,ibin],res))

panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.98
panelres.nglPanelBottom                      = 0.02
panelres.nglPaperOrientation = "Auto"

print 'about to panel'


Ngl.panel(wks,plot,[math.ceil(float(histnbins)/3.0),3],panelres)

print 'panelled'



#Now plot log-log plots
wkres = Ngl.Resources()
wks = Ngl.open_wks(wks_type,FigDir + figtitle3,wkres)
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

if (plotdensity):
        res.tiYAxisString = "frequency"
else:
        res.tiYAxisString = "number"

for ibin in range(0,histnbins):
        res.tiMainString = 'Precip annual values for bin '  + str(precip_bins_plot[ibin])
        plot.append(Ngl.xy(wks,range(startyr,endyr+1),precipplot[:,ibin],res))

panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.98
panelres.nglPanelBottom                      = 0.02
panelres.nglPaperOrientation = "Auto"

print 'about to panel'


Ngl.panel(wks,plot,[math.ceil(float(histnbins)/3.0),3],panelres)

print 'panelled'



"""
Ngl.end()









