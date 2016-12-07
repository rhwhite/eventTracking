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
from scipy import interpolate
from rhwhitepackages.readwrite import shiftlons
from rhwhitepackages.readwrite import XrayOpen
from rhwhitepackages.stats import regressmaps
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
parser.add_argument('--splittype',metavar='splittype',type=str,nargs='?',default='day',help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs='?',default=4,help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound1',type=int,nargs='+',help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=int,nargs="+",help='upper bounds')
parser.add_argument('--unit',type=str,nargs='?',default='day',help='units of split type')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--startyr',metavar='startyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--endyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--minlat',type=int,nargs='?',default=-45,help='min lat')
parser.add_argument('--maxlat',type=int,nargs='?',default=45,help='max lat')
parser.add_argument('--minlon',type=int,nargs='?',default=0,help='min lon')
parser.add_argument('--maxlon',type=int,nargs='?',default=360,help='max lon')
parser.add_argument('--sumlats',type=int,nargs='?',default=-1,help='number of lats to sum over')
parser.add_argument('--sumlons',type=int,nargs='?',default=-1,help='number of lons to sum over')
parser.add_argument('--minGB',type=int,nargs='?',default=-1,help='minimum number of gridboxes to count as an event')


args = parser.parse_args()

print args.Data

splittype = args.splittype
speedtspan = args.speedtspan
tbound1 = args.tbound1
tbound2 = args.tbound2
unit = args.unit
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
minGB = args.minGB

if splittype=='day':
	splitname = 'Sizes'
elif splittype =='maxspeed':
	splitname = 'MaxSpeeds_' + str(speedtspan) + 'ts'

nbounds = len(tbound1)

FillValue = -9999

# Time period for analysis
nyears = endyr - startyr + 1
mapping = "center"
FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' + Version + str(startyr) + '/proc/'

if Data == "TRMM":
	Fstartyr = 1998
	Fendyr = 2014

        PrecipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"

	if sumlats > 0:
		PrecipClimFile = "Regrid_" + str(sumlats) + "_" + str(sumlons) + '_TRMM_3B42_1998-2014_annmean.nc'
	else:
                PrecipClimFile = 'TRMM_3B42_1998-2014_annmean.nc' 

	if Version == '5th_nanto25':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
	elif Version == '5th_nantozero':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
	elif Version == '7thresh':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
	elif Version == '6th_from6':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
	elif Version == '5th_from48':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'
elif Data == "TRMMERAIgd":
        if Version in ["ERAIgd","2mmhr"]:
		PrecipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"
                if sumlats > 0:
			PrecipClimFile = "Regrid_" + str(sumlats) + "_" + str(sumlons) + "_regrid2ERAI_TRMM_3B42_1998-2014_annmean.nc"
		else:
			PrecipClimFile = "regrid2ERAI_TRMM_3B42_1998-2014_annmean.nc"
elif Data == "ERAI":
	PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
	if sumlats > 0:
		PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ncra_ERAI_Totalprecip_1980-2015_preprocess.nc' #'SeasAnn_ERAI_Totalprecip_1980-2015_preprocess.nc'
	else:
                PrecipClimFile = 'ERAI_Totalprecip_1980-2015_annmean.nc'

elif Data == "ERA20C":
        PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERA_20C/'
        if sumlats > 0:
		PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ERA_20C_Ann_Totalprecip_' + str(startyr) + '-' + str(endyr) + '.nc'
	else:
                PrecipClimFile = 'ERA_20C_Ann_Totalprecip_' + str(startyr) + '-' + str(endyr) + '.nc'

elif Data == "CESM":
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '/proc/'
	
	PrecipClimDir = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
	if sumlats > 0:
		PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'
        else:
                PrecipClimFile = 'ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'


FileInPrecip = XrayOpen(PrecipClimDir + PrecipClimFile)

if Data in ['TRMM']:
	lats = FileInPrecip['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
	invar = 'pcp'
	lons = FileInPrecip['lon']
elif Data in ['TRMMERAIgd']:
        lats = FileInPrecip['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
        invar = 'pcp'

elif Data in ["ERAI"]:
        lats = FileInPrecip['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
	invar = 'tpnew'

elif Data in ["ERA20C"]:
	lats = FileInPrecip['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
	invar = 'tpnew'

elif Data in ["CESM"]:
        lats = FileInPrecip['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
	invar = 'PRECT'


if lats[0] > lats[1]:
	PrecipIn = FileInPrecip[invar][:,::-1,:].sel(lat=slice(MinLatF,MaxLatF),time=slice(str(startyr),str(endyr))) #["PrecipClimAnn"]
else:
	PrecipIn = FileInPrecip[invar].sel(lat=slice(MinLatF,MaxLatF),time=slice(str(startyr),str(endyr))) #["PrecipClimAnn"]

try:
	print np.amax(PrecipIn)
	if PrecipIn.units == "mm/hr":
		print 'converting'
		PrecipIn = PrecipIn * 24.0 #convert to mm/day
	elif PrecipIn.units == "mm/day":
		pass	
	else:
		error("unexpected unit in Precip file")
except AttributeError:
	if np.amax(PrecipIn) < 1.0:
		print "guessing we need to convert precip units! You should check this!"
		PrecipIn = PrecipIn * 24.0 #convert to mm/day
	else:
		print "guessing that we don't need to convert precip units - you should check this!"


def getERAdata(ERAfile):
	ERAPrecipFile = XrayOpen(ERAfile)
	ERAlats = ERAPrecipFile['latitude']

	if ERAlats[0] > ERAlats[1]:
		ERAPrecipIn = ERAPrecipFile['tpnew'][:,::-1,:].sel(latitude=slice(MinLatF,MaxLatF))
	else:
		ERAPrecipIn = ERAPrecipFile['tpnew'].sel(latitude=slice(MinLatF,MaxLatF))

	ERAlats = ERAPrecipIn.latitude
	ERAlons = ERAPrecipIn.longitude
	ERAyears = ERAPrecipIn.time['time.year'].values

	return(ERAPrecipIn, ERAlats, ERAlons, ERAyears)

ERAIPrecipClimDir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
ERAIPrecipClimFile = 'ERAI_Totalprecip_1980-2015_annmean.nc'

ERAIPrecipIn, ERAIlats, ERAIlons, ERAIyears = getERAdata(ERAIPrecipClimDir + ERAIPrecipClimFile)

ERA20PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERA_20C/'
ERA20PrecipClimFile = 'ERA_20C_Ann_Totalprecip_1980-2011.nc'

ERA20PrecipIn, ERA20lats, ERA20lons, ERA20years = getERAdata(ERA20PrecipClimDir + ERA20PrecipClimFile)

ERAInyears = ERAIyears.shape[0]
ERA20nyears = ERA20years.shape[0]


def regressmaps(m,c,r,p,stderr,A,linreg,nlats,nlons):
# Calculate linear trend
	for ilat in range(0,nlats):
		for ilon in range(0,nlons):
			m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(A,linreg[:,ilat,ilon])

def plotline(nlines,datain,yearsin,nsep,title,Datatitle,figtitlein,LineTitles):

        wkres = Ngl.Resources()
        wkres.wkColorMap = "MPL_BrBG"
        wks_type = "eps"
        wksL = Ngl.open_wks(wks_type,figtitlein,wkres)
        res = Ngl.Resources()
        res.nglFrame = False
        res.nglDraw = False
        res.xyMarkLineMode = "Markers"
        res.xyMonoMarkLineMode = True
        res.xyMarkerColor = "blue"
        res.xyMarkers = 16
        res.xyMarkerSizeF = 0.01
        res.xyYStyle = "Linear"
        res.xyXStyle = "Linear"

        res.tiMainFontHeightF = 0.035
        res.tiYAxisFontHeightF = 0.03
        res.tiXAxisOn = False
        res.tmXBLabelFontHeightF = 0.03
        res.tmYLLabelFontHeightF = 0.03
        res.tmYLMode = "Automatic"
        res.tmYLFormat = "@6^g"

        res.vpWidthF = 0.9
        res.vpHeightF = 0.9


        plot = []

	nyearsin = len(yearsin)
        yearnums = range(0,nyearsin)
        A = np.array(yearnums)
        regress = (np.zeros((5,nsep),np.float))

	for iline in range(0,nlines):
		for ibin in range(0,nsep):
			if ibin == 0:
				res.tiYAxisString = LineTitles[iline]
			else:
				res.tiYAxisString = ""

			linreg = datain[iline][ibin][:]
			regress[:,ibin] = stats.linregress(A,linreg)

			if ibin == nsep -1:
				res.tiMainString = '{:^80}'.format('          >' + '{:2.1g}'.format(tbound1[ibin]) + unit + ' thr: ' + str(FracThreshold[ibin]) +  '; r2 = ' + '{:5.3f}'.format(regress[2,ibin]*regress[2,ibin]) + '; p=' + '{:5.3f}'.format(regress[3,ibin]) + "           ")
			else:
				res.tiMainString = '{:^80}'.format('{:2.1g}'.format(tbound1[ibin]) + '-' + '{:2.1g}'.format(tbound2[ibin]) + unit + ' thr: ' + str(FracThreshold[ibin]) + '; r2 = ' + '{:5.3f}'.format(regress[2,ibin]*regress[2,ibin]) + '; p=' + '{:5.3f}'.format(regress[3,ibin]))

			plot.append(Ngl.xy(wksL,yearsin,datain[iline,ibin,:],res))

        panelres = Ngl.Resources()
        panelres.nglPanelYWhiteSpacePercent = 8.
        panelres.nglPanelXWhiteSpacePercent = 0.0

        if nlines > 3:
                panelres.nglPanelTop                      = 0.8
                panelres.nglPanelBottom                      = 0.15
        else:
                panelres.nglPanelTop                      = 0.95
                panelres.nglPanelBottom                      = 0.01

        #panelres.nglScale = False
        panelres.nglMaximize = True

        txres = Ngl.Resources()
        txres.txFontHeightF = 0.012
        Ngl.text_ndc(wksL,'Annual timeseries of global number of events of various scales from' + Datatitle,0.5,0.94,txres)
        Ngl.panel(wksL,plot,[nlines,float(nsep)],panelres)


def plotmap(plotvars1,plotvars2,plotmin1,plotmax1,plotmin2,plotmax2,vartitle1,vartitle2,title,figtitle,lons,lats,minlon,maxlon,minlat,maxlat):
	nplots = plotvars1.shape[0]
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
	if lons[0] <= 0:
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
        #temp,lonsnew = Ngl.add_cyclic(plotvars1[0], lons)

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
	        tempplot,lonsnew = Ngl.add_cyclic(plotvars1[iplot], lons)

		tempplot[np.where(np.isnan(tempplot))] = FillValue
		res.cnMinLevelValF       = plotmin1[iplot]          # contour levels.
		res.cnMaxLevelValF       = plotmax1[iplot]
		res.cnLevelSpacingF      = ((plotmax1[iplot]-plotmin1[iplot])/10.0)
		res.tiMainString = vartitle1[iplot] + "; mean = {:2.3g}".format(np.nanmean(tempplot))
		toplot.append(Ngl.contour_map(wks,tempplot,res))

                tempplot = plotvars2[iplot]
                tempplot[np.where(np.isnan(tempplot))] = FillValue
                res.cnMinLevelValF       = plotmin2[iplot]          # contour levels.
                res.cnMaxLevelValF       = plotmax2[iplot]
                res.cnLevelSpacingF      = ((plotmax2[iplot]-plotmin2[iplot])/10.0)
		res.tiMainString = vartitle2[iplot] + "; mean = {:2.3g}".format(np.nanmean(tempplot))
                toplot.append(Ngl.contour_map(wks,tempplot,res))

	
	textres = Ngl.Resources()
	textres.txFontHeightF = 0.015
	Ngl.text_ndc(wks,title,0.5,0.87,textres)


	panelres = Ngl.Resources()
	panelres.nglPanelLabelBar = True
	#panelres.nglPanelYWhiteSpacePercent = 5.
	#panelres.nglPanelXWhiteSpacePercent = 5.

	panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
	if nplots > 3:
		panelres.nglPanelTop                      = 0.8
		panelres.nglPanelBottom                      = 0.15
	else:
                panelres.nglPanelTop                      = 0.95
                panelres.nglPanelBottom                      = 0.01


	#panelres.nglPanelFigureStrings            = ["a","b","c","d","e","f"]
	#panelres.nglPanelFigureStringsJust        = "BottomLeft"

	panelres.nglPaperOrientation = "Auto"

	plot = Ngl.panel(wks,toplot,[nplots,2],panelres)

# Get lats and lons
iday = 0
if minGB > 0:
	fileadd = '_min' + str(minGB) + 'GB'
else:
	fileadd = ''

if sumlats > 0:
	FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
else:
	FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd +  '.nc'

FileIn = XrayOpen(DirIn + FileIn)

lats = FileIn['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
lons = FileIn['lon'].values
years = FileIn['year'].sel(year=slice(startyr,endyr)).values

nlats = len(lats)
nlons = len(lons)
nyears = len(years)

difflat = lats[1]-lats[0]
difflon = lons[1]-lons[0]

if splittype == 'day':
	arraynbounds = nbounds + 3
else:
	arraynbounds = nbounds

mmperevent = np.zeros([arraynbounds,nyears,nlats,nlons])      # +1 for total sum as first plot
eventspermm = np.zeros([arraynbounds,nyears,nlats,nlons])      # +1 for total sum as first plot
Precipdays = np.zeros([arraynbounds,nyears,nlats,nlons])   # + 1 for total precip as first plot
Events = np.zeros([arraynbounds,nyears,nlats,nlons])
RegressPrecip = np.zeros([arraynbounds,nlats,nlons])
RegressEvents = np.zeros([arraynbounds,nlats,nlons])

temp = np.zeros([nlats,nlons],float)

A = PrecipIn.time['time.year'].values


titlesP = []
titlesE = []
titlesF = []

titlesPrg = []
titlesErg = []

for iday in range(0,nbounds):
	tbounds = str(tbound1[iday]) + '-' + str(tbound2[iday])
	if sumlats > 0:
		FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + tbounds + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
	else:
                FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + tbounds + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd + '.nc'

	#Get lons and lats
	FileIn = xray.open_dataset(DirIn + FileIn)

	mmperevent[iday,:,:,:] = FileIn['TPrecip']/FileIn['TDensity'].sel(lat=slice(MinLatF,MaxLatF),year=slice(startyr,endyr)) 

        eventspermm[iday,:,:,:] = FileIn['TDensity']/FileIn['TPrecip'].sel(lat=slice(MinLatF,MaxLatF),year=slice(startyr,endyr))


	Precipdays[iday,:,:,:] = FileIn['TPrecip'].sel(lat=slice(MinLatF,MaxLatF),year=slice(startyr,endyr))
	Events[iday,:,:,:] = FileIn['TDensity'].sel(lat=slice(MinLatF,MaxLatF),year=slice(startyr,endyr))

	regressmaps(RegressPrecip[iday,:,:],temp,temp,temp,temp,A,Precipdays[iday,:,:,:],nlats,nlons)
	titlesPrg.append('Lin reg. of precip in ' + tbounds + ' day events')
        regressmaps(RegressEvents[iday,:,:],temp,temp,temp,temp,A,Events[iday,:,:,:],nlats,nlons)
        titlesErg.append('Lin reg. of number of ' + tbounds + ' day events')

        titlesP.append('Precip in ' + tbounds + ' day events')
        titlesE.append('Number of ' + tbounds + ' day events')
	titlesF.append('Fraction of precip in '+ tbounds + ' day events')

#DenAnnAvg = np.nanmean(Dendays,axis=1)
Precipdays[nbounds,:,:,:] = PrecipIn * 365.25
PrecipAllDays = np.nansum(Precipdays[0:nbounds,:,:,:],axis=0)	# Total precip captured in events

PrecipNC = PrecipIn * 365.25 - PrecipAllDays	# Annual: total precip NOT captured in events

Events[nbounds,:,:,:] = np.nansum(Events[0:nbounds,:,:,:],axis=0)

PrecipAnnAvg = np.nanmean(Precipdays,axis=1)	# mean over years
EventsAnnAvg = np.nanmean(Events,axis=1)
EventspermmAvg = np.nanmean(eventspermm,axis=1)

PrecipAll = np.nansum(PrecipAnnAvg[0:nbounds,:,:],axis=0)

PrecipAnom = PrecipIn - np.nanmean(PrecipIn,axis=0) 	# mm/day
PrecipCAnom = np.nanmean(PrecipNC,axis=0) - PrecipNC


PrecipAllPercent = 100.0 * np.divide(PrecipAll,np.nanmean(PrecipIn*365.0,axis=0))	# convert PrecipIn from mm/day to mm/yr before dividing
PrecipAllFrac = np.divide(PrecipAll,np.nanmean(PrecipIn*365.0,axis=0))       # convert PrecipIn from mm/day to mm/yr before dividing

# Mask out region with >=50% from each type of event

# Regressions of totals from events
SumPrecipPercent = 100.0 * np.nansum(PrecipAll)/np.nansum(PrecipIn * 365.0)
PrecipPercent = np.zeros(PrecipAnnAvg.shape)
PrecipFrac = np.zeros(PrecipAnnAvg.shape)
PrecipMask = np.zeros([nbounds,nlats,nlons])

RegressPrecip[nbounds,:,:] = np.nansum(RegressPrecip[0:nbounds,:,:],axis=0)
titlesPrg.append('Linear reg. of total precip in events ' + str(startyr) + '-' + str(endyr))
titlesErg.append('blank plot')


# Regressions of totals
regressmaps(RegressPrecip[nbounds+1,:,:],temp,temp,temp,temp,A,PrecipIn*365.25,nlats,nlons)
titlesPrg.append('Linear reg. of total precip ' + str(startyr) + '-' + str(endyr))
titlesP.append('Total precip ' + str(startyr) + '-' + str(endyr))
regressmaps(RegressEvents[nbounds+1,:,:],temp,temp,temp,temp,A,np.nansum(Events,axis=0),nlats,nlons)
titlesErg.append('Linear reg. of total events ' + str(startyr) + '-' + str(endyr))
titlesE.append('Total events ' + str(startyr) + '-' + str(endyr))

# Regressions of total NOT captured by events
regressmaps(RegressPrecip[nbounds+2,:,:],temp,temp,temp,temp,A,PrecipNC,nlats,nlons)
titlesPrg.append('Linear reg. of non-event precip ' + str(startyr) + '-' + str(endyr))
titlesErg.append('blank plot')


arrayindex = 0
#	#DenPercent[0,:,:] = DenAll / (difflat * difflon)
#	PrecipPercent[0,:,:] = PrecipAllPercent	# Percent of total TRMM precip falling in these events
#	#titlesDen.append("Total annual event density, events/(yr deg~S1~2 )")
#	#titlesPrec.append("Percentage of total precipitation captured in events")
#	arrayindex += 1


FracThreshold = [0.5,0.25,0.2,0.1]
tboundlist = []
for iday in range(0,nbounds):
        tboundlist.append(str(tbound1[iday]) + '-' + str(tbound2[iday]))
	#DenPercent[iday+arrayindex,:,:] = 100.0 * (np.where(DenAll > 0, np.divide(DenAnnAvg[iday,:,:],DenAll),0.0)) 
        PrecipPercent[iday,:,:] = 100.0 * (np.where(PrecipAll > 0, np.divide(PrecipAnnAvg[iday,:,:],PrecipAll),0.0))
	PrecipFrac[iday,:,:] = np.where(PrecipAll > 0, np.divide(PrecipAnnAvg[iday,:,:],PrecipAll),0.0)
	PrecipMask[iday,:,:] = np.where(PrecipFrac[iday,:,:]<FracThreshold[iday],0.0,1.0)

print tboundlist
PrecipMask = xray.DataArray(PrecipMask,coords=[('days',tboundlist),('lat',lats),('lon',lons)])
print np.sum(PrecipMask,axis=(1,2))

ERAIPrecipMask = PrecipMask.sel(method = 'nearest',lat = ERAIlats,lon=ERAIlons)

# And now plot 

figtitlein = FigDir + 'PrecipFrac_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit + '_min' + fileadd
titlein = 'Events in ' + Data + " " + Version

clims1,clims2,clims3,clims4 = [-0.0,0.0,0.0,0.0,0.0,-40.0,-40.0],[2500.0,1000.0,500.0,500.0,40.0,40.0,40.0],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[1.0,1.0,1.0,1.0,1.0,1.0,60.0]

plotvar1 = np.copy(PrecipAnnAvg[0:nbounds,:,:])
plotvar2 = np.copy(PrecipFrac[0:nbounds,:,:])
plotmap(plotvar1,plotvar2,clims1,clims2,clims3,clims4 ,titlesP,titlesF,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

"""
figtitlein = FigDir + 'ConstructPrecip_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit + '_min' + fileadd
titlein = 'Events in ' + Data + " " + Version

clims1,clims2,clims3,clims4 = [-40.0,-20.0,-10.0,-40.0,-40.0,-40.0,-40.0],[40.0,20.0,10.0,40.0,40.0,40.0,40.0],[-60.0,-1.0,-0.1,-60.0,-60.0,-60.0,-60.0],[60,1.0,0.1,60.0,60.0,60.0,60.0]

plotvar1 = np.copy(RegressPrecip[:,:,:])
plotvar2 = np.copy(RegressEvents[:,:,:])
plotmap(plotvar1,plotvar2,clims1,clims2,clims3,clims4 ,titlesPrg,titlesErg,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)
"""

"""
figtitlein = FigDir + 'AnnualPrecipDensity_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit + '_min' + fileadd
titlein = 'Events in ' + Data + " " + Version
clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,0.0,0.0,0.0],[1500.0,1000.0,500.0,500.0,3000.0,500.0],[0.0,0.0,0.0,0.0,0.0,0.0],[5000,50,5,1.0,5000,60.0]

plotvar1 = np.copy(PrecipAnnAvg[0:nbounds+1,:,:])
plotvar2 = np.copy(EventsAnnAvg[0:nbounds+1,:,:])
plotmap(plotvar1,plotvar2,clims1,clims2,clims3,clims4 ,titlesP,titlesE,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)
"""

figtitlein = FigDir + 'ConstructPrecipReg_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit + '_min' + fileadd
titlein = 'Events in ' + Data + " " + Version
clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,0.0,0.0,0.0],[1500.0,1000.0,500.0,500.0,3000.0,500.0],[0.0,0.0,0.0,0.0,0.0,0.0],[5000,50,5,1.0,5000,60.0]

createPrg = np.zeros(PrecipAnnAvg.shape)
titlesCPrg = []
for ibound in range(0,nbounds):
        tbounds = str(tbound1[ibound]) + '-' + str(tbound2[ibound])
	createPrg[ibound,:,:] = RegressPrecip[nbounds,:,:] * (PrecipFrac[ibound,:,:])


	titlesCPrg.append('Estimated lin reg. of precip in ' + tbounds + ' day events')
clims1,clims2 = [-40.0,-20.0,-10.0,-10.0,-40.0,0.0],[40.0,20.0,10.0,10.0,40.0,10.0]
clims3 = clims1
clims4 = clims2
titlesCPrg.append("blank plot")
plotvar1 = np.copy(RegressPrecip[0:nbounds+1,:,:])
plotvar2 = np.copy(createPrg[0:nbounds+1,:,:])
plotmap(plotvar1,plotvar2,clims1,clims2,clims3,clims4 ,titlesPrg,titlesCPrg,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)


# Now plot annual Precip anomaly from TRMM: PrecipAnom[:,:,:] 
# * fraction of rain incorporated in events (PrecipAllPercent)
# * fraction of rain in x-y day events (PrecipPercent[iday,:,:])
# * events per mm/hr:  EventspermmAvg[iday,:,:]

globmeaneventanom = np.zeros([5,nbounds,nyears],float)

Linetitle = []

for iday in range(0,nbounds):
	EventAnom = 365.25 * PrecipAnom * PrecipFrac[iday,:,:] * EventspermmAvg[iday,:,:]	# convert from events/day to events/year
	globmeaneventanom[0,iday,:] = np.nansum(EventAnom,axis=(1,2))
	Linetitle.append('constructed number of events')

        CaptEventAnom = PrecipCAnom * PrecipFrac[iday,:,:] * EventspermmAvg[iday,:,:]
        globmeaneventanom[1,iday,:] = np.nansum(CaptEventAnom,axis=(1,2))
	Linetitle.append('number of captured events')
	#globmeaneventanom[2,iday,:] = np.nansum(EventAnom,axis=(1,2)) + np.nansum(CaptEventAnom,axis=(1,2))

        globmeaneventanom[2,iday,:] = np.nansum(Precipdays[iday,:,:,:],axis=(1,2))
	Linetitle.append('precip in events')

	globmeaneventanom[3,iday,:] = np.nansum(PrecipMask[iday,:,:].values * PrecipIn*365.0,axis=(1,2))
	Linetitle.append('constructed precip in events')

	
        globmeaneventanom[4,iday,:] = np.nansum(ERAIPrecipMask[iday,:,:].values * ERAIPrecipIn.sel(time=slice(str(startyr),str(endyr))) *365.0,axis=(1,2))
        Linetitle.append('constructed ERAI precip in events')

globmeaneventanom[3,nbounds-1,:] = 365.25 * np.nansum(PrecipAnom,axis=(1,2))
globmeaneventanom = np.where(np.isfinite(globmeaneventanom),globmeaneventanom,0)
#print globmeaneventanom
figtitlein = FigDir + 'ConstructPrecipTrends_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit + '_min' + fileadd + '_FracThreshVariable'
titlein = 'Events in ' + Data + " " + Version

plotline(5,globmeaneventanom,years,nbounds,'Constructed trends','TRMM',figtitlein,Linetitle)


ERAIglobmeananom = np.zeros([2,nbounds,ERAInyears],float)
ERAILinetitle = []

interpERAIPrecipIn = ERAIPrecipIn.sel(latitude = lats,longitude=lons,method='nearest')

interpERAIPrecipAnom = interpERAIPrecipIn - np.nanmean(interpERAIPrecipIn,axis=0)     # mm/day


for iday in range(0,nbounds):
        ERAIglobmeananom[0,iday,:] = np.nansum(ERAIPrecipMask[iday,:,:].values * ERAIPrecipIn.sel(time=slice('1980','2015')) *365.0,axis=(1,2))
        ERAILinetitle.append('constructed ERAI precip in events')

        EventAnom = 365.25 * interpERAIPrecipAnom * PrecipFrac[iday,:,:] * EventspermmAvg[iday,:,:]       # convert from events/day to events/year
        ERAIglobmeananom[1,iday,:] = np.nansum(EventAnom,axis=(1,2))
        ERAILinetitle.append('constructed ERAI number of events')



ERAIglobmeananom = np.where(np.isfinite(ERAIglobmeananom),ERAIglobmeananom,0)
figtitlein = FigDir + 'ConstructPrecipTrendsERAI_' + Data + '_' + Version + '_' + str(startyr) + '-' + str(endyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit + '_min' + fileadd + '_FracThreshVariable'
titlein = 'Events in ' + Data + " " + Version

plotline(2,ERAIglobmeananom,ERAIyears,nbounds,'Constructed trends','ERAI',figtitlein,ERAILinetitle)



Ngl.end()




