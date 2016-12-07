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
import math
from scipy import stats
from rhwhitepackages.readwrite import shiftlons
from rhwhitepackages.readwrite import XrayOpen
import argparse
import resource

rsrcV = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit starts as  :', soft
print 'Hard limit starts as :', hard
resource.setrlimit(rsrcV, (50000000000, hard)) #limit memory usage
#                          137438953472
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit changed to :', soft

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--minlat',type=int,nargs='?',default=-45,help='min lat')
parser.add_argument('--maxlat',type=int,nargs='?',default=45,help='max lat')
parser.add_argument('--minlon',type=int,nargs='?',default=0,help='min lon')
parser.add_argument('--maxlon',type=int,nargs='?',default=360,help='max lon')
parser.add_argument('--sumlats',type=int,nargs='?',default=1,help='number of lats to sum over')
parser.add_argument('--sumlons',type=int,nargs='?',default=1,help='number of lons to sum over')


args = parser.parse_args()

print args.Data

Data = args.Data[0]
Version = args.Version[0]
startyr = args.startyr[0]
endyr = args.endyr[0]
MinLonF = args.minlon
MaxLonF = args.maxlon
MinLatF = args.minlat
MaxLatF = args.maxlat
sumlats = args.sumlats
sumlons = args.sumlons

FillValue = -9999

# Time period for analysis
nyears = endyr - startyr + 1
print nyears
mapping = "center"
FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

if Data == "TRMM":
        PrecipDir = "/home/disk/eos4/rachel/Obs/TRMM/"

	if Version == "ERAIgd":
		Filename = "Regrid_" + str(sumlats) + "_" + str(sumlons) + "_regrid2ERAI_TRMM_3B42_1998-2014_annclim.nc" 
	else:
	        Filename = "TRMM_3B42_1998-2014_annclim.nc" 

elif Data == "CESM":
	PrecipDir = "/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/"
	Filename = "ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc"
else:
	print "not set up for " + Data + " yet!"

FileInPrecip = XrayOpen(PrecipDir + Filename)

if Data in ['TRMM']:
	latin = FileInPrecip['latitude'].sel(latitude=slice(MinLatF,MaxLatF))
	lonin = FileInPrecip['longitude']
	invar = 'pcp'
        if latin[0] > latin[1]:
                PrecipIn = FileInPrecip[invar][:,::-1,:].sel(latitude=slice(MinLatF,MaxLatF)) #["PrecipClimAnn"]
        else:
                PrecipIn = FileInPrecip[invar].sel(latitude=slice(MinLatF,MaxLatF)) #["PrecipClimAnn"]

elif Data in ['CESM']:
        latin = FileInPrecip['lat'].sel(lat=slice(MinLatF,MaxLatF))
        lonin = FileInPrecip['lon']
        invar = 'PRECT'

        if latin[0] > latin[1]:
                PrecipIn = FileInPrecip[invar][:,::-1,:].sel(lat=slice(MinLatF,MaxLatF)) #["PrecipClimAnn"]
        else:
                PrecipIn = FileInPrecip[invar].sel(lat=slice(MinLatF,MaxLatF)) #["PrecipClimAnn"]

try:
	print np.amax(PrecipIn)
	if PrecipIn.units == "mm/hr":
		pass
	elif PrecipIn.units == "mm/day":
                print 'converting'
		PrecipIn = PrecipIn / 24.0 #convert to mm/hr
	elif PrecipIn.units == "m/s":
		PrecipIn = PrecipIn * 1000.0 * 60.0 * 60.0	# convert from m/s to mm/hr
	elif PrecipIn.units == "":
		if np.amax(PrecipIn) < 200.0:
			print "guessing we don't need to convert precip units! You should check this!"
		else:
			print "guessing that we need to convert precip units - you should check this!"
                        PrecipIn = PrecipIn / 24.0 #convert to mm/hr

	else:
		exit("unexpected unit in Precip file:" + PrecipIn.units)
except AttributeError:
	if np.amax(PrecipIn) < 1.0:
		print "guessing we need to convert precip units! You should check this!"
		PrecipIn = PrecipIn * 24.0 #convert to mm/day
	else:
		print "guessing that we don't need to convert precip units - you should check this!"

ndays = PrecipIn.shape[0]
print ndays

def plotmap(plotvars,plotmin1,plotmax1,title,figtitle,lons,lats,minlon,maxlon,minlat,maxlat):
	nplots = plotvars.shape[0]
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

	
	temp,lons = Ngl.add_cyclic(plotvars[0],lons)	# add to lons
	nlons = len(lons)
	# if lons start negative, shift everything over so there isn't a line down the middle of the Pacific
	#if lons[0] < 0:
	#	nlonhalf = nlons/2
	#	lonsnew = np.zeros(lons.shape,np.float)
	#	lonsnew[0:nlonhalf] = lons[nlonhalf:nlons]
	#	lonsnew[nlonhalf:nlons] = lons[0:nlonhalf] + 360.0
	#	lons = lonsnew
	#
	#	for iplot in range(0,nplots):
	#		plotvars[iplot] = shiftlons(plotvars[iplot],lons)
	#else:
	#	lonsnew = lons

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


	for iplot in range(0,nplots):
		print iplot
		#tempplot = plotvars[iplot,:,:]
		tempplot = Ngl.add_cyclic(plotvars[iplot])
		res.cnMinLevelValF       = plotmin1[iplot]          # contour levels.
		res.cnMaxLevelValF       = plotmax1[iplot]
		res.cnLevelSpacingF      = ((plotmax1[iplot]-plotmin1[iplot])/10.0)
		if iplot == 0:
			res.tiMainString = title
		else:
			res.tiMainString = ""
		toplot.append(Ngl.contour_map(wks,tempplot,res))

	
	textres = Ngl.Resources()
	textres.txFontHeightF = 0.015
	#Ngl.text_ndc(wks,title,0.5,0.87,textres)


	panelres = Ngl.Resources()
	panelres.nglPanelLabelBar = True
	#panelres.nglPanelYWhiteSpacePercent = 5.
	#panelres.nglPanelXWhiteSpacePercent = 5.

	panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
	if nplots > 5:
		panelres.nglPanelTop                      = 0.8
		panelres.nglPanelBottom                      = 0.15
	else:
                panelres.nglPanelTop                      = 0.95
                panelres.nglPanelBottom                      = 0.01


	#panelres.nglPanelFigureStrings            = ["a","b","c","d","e","f"]
	#panelres.nglPanelFigureStringsJust        = "BottomLeft"

	panelres.nglPaperOrientation = "Auto"

	plot = Ngl.panel(wks,toplot,[nplots,1],panelres)
	Ngl.destroy(wks)
# Get lats and lons
iday = 0

titles = []

PrecipIn[np.where(np.isnan(PrecipIn))] = FillValue
PrecipIn = np.where(PrecipIn < 0.00001,FillValue,PrecipIn)

title = (Data + ' annual mean Precip, mm/hr, 1998-2014')
figtitlein = FigDir + Data + '_annualmean_1998-2014' 

clims1,clims2 = [0.0],[0.5]

plotmap(PrecipIn[iday:iday+1,:,:],clims1,clims2,title,figtitlein,lonin,latin,MinLonF,MaxLonF,MinLatF,MaxLatF)



Ngl.end()




