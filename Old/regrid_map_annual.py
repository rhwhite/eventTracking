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

sumlats = 5
sumlons = 10

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


FileI = 'DensityMap_' + str(startyr) + '-' + str(endyr) + '_5th_nanto25.nc'
FileO = 'DensityMap_regrid' + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_5th_nanto25.nc'

filetimespan = "3hrly"

file = xray.open_dataset(DirI + FileI)
lons = file['Longitude'].values
lats = file['Latitude'].values
years = file['years'].values
DenMap = file['DensityMap'].values

nlons = len(lons)
nlats = len(lats)

nlonsnew = nlons/sumlons
nlatsnew = nlats/sumlats

print nlonsnew
print nlatsnew

Latsnew = np.zeros(nlatsnew,np.float)
Lonsnew = np.zeros(nlonsnew,np.float)

RDenMapnew = np.zeros((nyears,4,nlatsnew,nlonsnew),np.int)

#'regrid' by summing over large boxes (averaging lons and lats)
inlat = 0
for ilats in range(0,nlats,sumlats):
	inlon = 0
	Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlats])
#	print np.mean(lats[ilats:ilats+sumlats])
	for ilons in range(0,nlons,sumlons):
		Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlons])
#		print np.mean(lons[ilons:ilons+sumlons])
		RDenMapnew[:,:,inlat,inlon] = np.sum(DenMap[:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(2,3),dtype=np.int)		
		inlon += 1
	inlat += 1

file.close()

ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('years', nyears)
ncfile.createDimension('size', 4)
ncfile.createDimension('lon', nlonsnew)
ncfile.createDimension('lat', nlatsnew)


ODensityMap = ncfile.createVariable('DensityMap','f4',('years','size','lat','lon'),fill_value=-9999)
OLongitude = ncfile.createVariable('Longitude','f4',('lon'),fill_value=-9999)
OLatitude = ncfile.createVariable('Latitude','f4',('lat'),fill_value=-9999)
OYears = ncfile.createVariable('years','f4',('years'),fill_value=-9999)
Diffs1DensityMap = ncfile.createVariable('DensityMap_Last-First','f4',('size','lat','lon'),fill_value=-9999)
Diffs2DensityMap = ncfile.createVariable('DensityMap_Last4-First4','f4',('size','lat','lon'),fill_value=-9999)
Diffs3DensityMap = ncfile.createVariable('DensityMap_Last8-First8','f4',('size','lat','lon'),fill_value=-9999)

DiffsPC2DensityMap = ncfile.createVariable('DensityMap_Last4-First4PC','f4',('size','lat','lon'),fill_value=-9999)


setattr(ODensityMap,'Extra Info','Size based on timespan: 0 is < 8 (1 day), 1: < 16 (2 days) < 48 (6 days), 2: > 48)')
setattr(ODensityMap,'History','Summation over large gridboxes: ' + str(sumlats) + ' lat and ' + str(sumlons) + ' lon')
setattr(Diffs1DensityMap,'Extra Info','Last year minus first year (excluding 1998 and 2014)')
setattr(Diffs2DensityMap,'Extra Info','Last 4 years minus first 4 years excluding first and last')
setattr(Diffs3DensityMap,'Extra Info','Last 7 years minus first 7 years excluding first and last')
setattr(DiffsPC2DensityMap,'Extra Info','Percentage of Last 4 years minus first 4 years, relative to first 4 years')



OYears[:] = years
OLatitude[:] = Latsnew
OLongitude[:] = Lonsnew
ODensityMap[:,:,:,:] = RDenMapnew

Diffs1DensityMap[:,:,:] = RDenMapnew[nyears-2,:,:,:] - RDenMapnew[1,:,:,:]
Diffs2DensityMap[:,:,:] = np.sum(RDenMapnew[nyears-5:nyears-1,:,:,:],axis=0) - np.sum(RDenMapnew[1:5,:,:,:],axis=0)
Diffs3DensityMap[:,:,:] = np.sum(RDenMapnew[nyears-7:nyears-1,:,:,:],axis=0) - np.sum(RDenMapnew[1:7,:,:,:],axis=0)

DiffsPC2DensityMap[:,:,:] = 100.0 * (np.sum(RDenMapnew[nyears-5:nyears-1,:,:,:],axis=0,dtype=np.float) - np.sum(RDenMapnew[1:5,:,:,:],axis=0,dtype=np.float))/np.sum(RDenMapnew[1:5,:,:,:],axis=0,dtype=np.float)




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









