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

Data = "TRMM"

Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

print Version

sumlats = 40 #10 16 is approximately 4 x 4 degrees for 0.25 degree data
sumlons = 40 #20

startyr = 1990 # Don't change - tied to file names!
endyr = 2014

nmonths = 12

seas = ['JFD','MAM','JJA','SON']
nseas = 4

plotdensity = False

minevent = 100000

if Data == "TRMM":

	startyr = 1998 # Don't change - tied to file names!
	endyr = 2014
	anstartyr = 1998 #year for analysis start
	anendyr = 2014 #year for analysis end
	
	if Version == 'Standard':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'

	elif Version == '5th_nanto25':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'

	elif Version == '5th_nantozero':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'

	elif Version == '7thresh':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
	elif Version == '6th_from6':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
	elif Version == '5th_from48':
		DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'
		DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'
	else:
        	sys.exit('unexpected Version')
	FileI = 'DenDirSpd_Map_monthly_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

elif Data == "ERAI":
        anstartyr = 1990 #year for analysis start
        anendyr = 2011 #year for analysis end
        

	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	DirO = DirI
	FileI = 'DenDirSpd_Map_monthly_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

elif Data == "CESM":
        anstartyr = 1990 #year for analysis start
        anendyr = 2011 #year for analysis end
        

        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'
        DirO = DirI
	FileI = 'DenDirSpd_Map_monthly_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

FileO = 'DenDirSpd_Map_monthly_regrid_' + Data + "_"  + str(startyr) + '-' + str(endyr) + '_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats_' + Version + '.nc'



nyears = anendyr - anstartyr + 1
mints = np.zeros(nyears)
maxts = np.zeros(nyears)



filetimespan = "3hrly"

print DirI + FileI
file = xray.open_dataset(DirI + FileI)
lons = file['Longitude'].values
lats = file['Latitude'].values
years = file['years'].values
print years
months = file['months'].values

nlons = len(lons)
nlats = len(lats)

print nlats, nlons
nlonsnew = math.floor(nlons/sumlons)
nlatsnew = math.floor(nlats/sumlats)

nlats2 = np.int(nlatsnew * sumlats)
nlons2 = np.int(nlonsnew * sumlons)

print nlons2
print nlats2

Latsnew = np.zeros(nlatsnew,np.float)
Lonsnew = np.zeros(nlonsnew,np.float)

# Create Outfile structure
ncfile = Dataset(DirO + FileO, 'w')
ncfile.createDimension('years', nyears) 
ncfile.createDimension('seas', nseas)
ncfile.createDimension('size', 4)
ncfile.createDimension('lon', nlonsnew)
ncfile.createDimension('lat', nlatsnew)
OLongitude = ncfile.createVariable('Longitude','f4',('lon'),fill_value=-9999)
OLatitude = ncfile.createVariable('Latitude','f4',('lat'),fill_value=-9999)
OYears = ncfile.createVariable('years','f4',('years'),fill_value=-9999)
OSeas = ncfile.createVariable('seas',str,('seas'))

OYears[:] = years
OSeas = seas

for variable in ["DensityMap","EastwardDensityMap","WestwardDensityMap","DirectionMap","DistanceMap","SpeedMap"]:
	VarMap = file[variable].values
	if variable in ["DensityMap","EastwardDensityMap","WestwardDensityMap"]:
		VarMapnew = np.zeros((nyears,nmonths,4,nlatsnew,nlonsnew),np.int)
	else:
                VarMapnew = np.zeros((nyears,nmonths,4,nlatsnew,nlonsnew),np.float64)

	#'regrid' by summing over large boxes (averaging lons and lats)
	inlat = 0
	for ilats in range(0,nlats2,sumlats):
		inlon = 0
		Latsnew[inlat] = np.mean(lats[ilats:ilats+sumlats])
	#	print np.mean(lats[ilats:ilats+sumlats])
		for ilons in range(0,nlons2,sumlons):
			Lonsnew[inlon] = np.mean(lons[ilons:ilons+sumlons])
	#		print np.mean(lons[ilons:ilons+sumlons])
			if variable in ["DensityMap","EastwardDensityMap","WestwardDensityMap"]:
				VarMapnew[:,:,:,inlat,inlon] = np.sum(VarMap[:,:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(3,4),dtype=np.int)		
			else:
				VarMapnew[:,:,:,inlat,inlon] = np.nanmean(VarMap[:,:,:,ilats:ilats+sumlats,ilons:ilons+sumlons],axis=(3,4),dtype=np.float64)	
			inlon += 1
		inlat += 1
	# Create seasonal and or means!

	if variable in ["DensityMap","EastwardDensityMap","WestwardDensityMap"]:
	        VarMapSeas = np.zeros((nyears,nseas,4,nlatsnew,nlonsnew),np.int)

		VarMapSeas[:,0,:,:,:] = np.nansum([VarMapnew[:,0,:,:,:],VarMapnew[:,1,:,:,:],VarMapnew[:,11,:,:,:]],axis=0)
		VarMapSeas[:,1,:,:,:] = np.nansum([VarMapnew[:,2,:,:,:],VarMapnew[:,3,:,:,:],VarMapnew[:,4,:,:,:]],axis=0)
		VarMapSeas[:,2,:,:,:] = np.nansum([VarMapnew[:,5,:,:,:],VarMapnew[:,6,:,:,:],VarMapnew[:,7,:,:,:]],axis=0)
		VarMapSeas[:,3,:,:,:] = np.nansum([VarMapnew[:,8,:,:,:],VarMapnew[:,9,:,:,:],VarMapnew[:,10,:,:,:]],axis=0)

	        VarMapAnn = np.nansum(VarMapnew,axis=1)

	else:
	        VarMapSeas = np.zeros((nyears,nseas,4,nlatsnew,nlonsnew),np.float64)

                VarMapSeas[:,0,:,:,:] = np.nanmean([VarMapnew[:,0,:,:,:],VarMapnew[:,1,:,:,:],VarMapnew[:,11,:,:,:]],axis=0)
                VarMapSeas[:,1,:,:,:] = np.nanmean([VarMapnew[:,2,:,:,:],VarMapnew[:,3,:,:,:],VarMapnew[:,4,:,:,:]],axis=0)
                VarMapSeas[:,2,:,:,:] = np.nanmean([VarMapnew[:,5,:,:,:],VarMapnew[:,6,:,:,:],VarMapnew[:,7,:,:,:]],axis=0)
                VarMapSeas[:,3,:,:,:] = np.nanmean([VarMapnew[:,8,:,:,:],VarMapnew[:,9,:,:,:],VarMapnew[:,10,:,:,:]],axis=0)
	
		VarMapAnn = np.nanmean(VarMapnew,axis=1)
	
	OVarMap = ncfile.createVariable(variable + "Ann",'f4',('years','size','lat','lon'),fill_value=-9999)
	OVarMapSeas = ncfile.createVariable(variable + "Seas",'f4',('years','seas','size','lat','lon'),fill_value=-9999)

	Diffs3VarMap = ncfile.createVariable(variable + "_Last8-First8",'f4',('size','lat','lon'),fill_value=-9999)

	SeasDiffs3VarMap = ncfile.createVariable(variable + "Seas_Last8-First8",'f4',('seas','size','lat','lon'),fill_value=-9999)


	setattr(OVarMap,'Extra Info','Size based on timespan: 0 is < 8 (1 day), 1: < 16 (2 days) < 48 (6 days), 2: > 48)')
	setattr(OVarMap,'History','Summation over large gridboxes: ' + str(sumlats) + ' lat and ' + str(sumlons) + ' lon')
	setattr(Diffs3VarMap,'Extra Info','Last 7 years minus first 7 years excluding first and last')
	setattr(SeasDiffs3VarMap,'Extra Info','Last 7 years minus first 7 years excluding first and last')


	OVarMap[:,:,:,:] = VarMapAnn
	OVarMapSeas[:,:,:,:,:] = VarMapSeas
	Diffs3VarMap[:,:,:] = np.nanmean(VarMapAnn[nyears-7:nyears-1,:,:,:],axis=0) - np.nanmean(VarMapAnn[1:7,:,:,:],axis=0)
	SeasDiffs3VarMap[:,:,:,:] = np.nanmean(VarMapSeas[nyears-7:nyears-1,:,:,:,:],axis=0) - np.nanmean(VarMapSeas[1:7,:,:,:,:],axis=0)

# After last variable, write lats and lons, which will be identical for all variables and just need to be written once
OLatitude[:] = Latsnew
OLongitude[:] = Lonsnew



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









