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
from scipy import stats

DataM = ["TRMM","ERA20C"]
VersionM = ["Standard","20Cstd"]
nrows = len(DataM)
print nrows

splittype = "day"
unit = "day"
tbound1 = [0,1,2,5]
tbound2 = [1,2,5,100]

nbounds = len(tbound1)

mapping = 'center'

MinLonF = 0
MaxLonF = 360
MinLatF = -20
MaxLatF = 20

sumlat = 1
sumlon = 1

test = 0

regresst = (np.zeros((5,nbounds),np.float))

toplot = []
toplotx = []

for idata in range(0,len(DataM)):
	Data = DataM[idata]
	Version = VersionM[idata]
	print Data, Version
	if Data == "CESM":
		Fstartyr = 1990
		Fendyr = 2014

		startyr = 1990
		endyr = 2014
		anendyr = 2011
		anstartyr = 1990
	elif Data == "TRMM":
		Fstartyr = 1998
		Fendyr = 2014
		startyr = 1998
		endyr = 2015

		anstartyr = 1998
		anendyr = 2014

	elif Data == "ERAI":
		Fstartyr = 1980
		Fendyr = 2014
		startyr = 1980
		endyr = 2014

		anstartyr = 1980
		anendyr = 2014
	elif Data == "ERA20C":
		Fstartyr = 1980
		Fendyr = 2011
		startyr = 1980
		endyr = 2011

		anstartyr = 1980
		anendyr = 2010

	if test == 1:
		endyr = startyr
	nyears = anendyr - anstartyr +1

	toplotx.append(range(anstartyr,anendyr+1))

	mints = np.zeros(nyears)
	maxts = np.zeros(nyears)

	plotdensity = False

	starttsteps = 0
	anntsteps = 2920 # timesteps per year

	minevent = 100000

	DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' + Version + str(startyr) + '/proc/'
	FileI1 = 'All_Precip_Sizes_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

	FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

	filetimespan = "3hrly"

	# Get lats and lons
	ibound = 0

	if splittype == "maxspeed":
		fileadd = "MaxSpeeds_" + str(speedtspan) + "ts_"
	elif splittype == "speed":
		fileadd = "Speeds_"
	elif splittype == "day":
		fileadd = "Sizes_"

	if tbound1[ibound] < 0:
		tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))
	else:
		tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))


        File = 'DenDirSpd_Map_Ann_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version +'.nc'

	print DirIn + File
	FileIn = xray.open_dataset(DirIn + File)

	lats = FileIn['lat']
	lons = FileIn['lon']
	years = FileIn['years'].sel(years = slice(anstartyr,anendyr))
	print years
	if lats[-1] < lats[0]:
		print("latitudes are not south to north!")

	nlats = len(lats)
	nlons = len(lons)
	nyears = len(years)
	print nyears

	Dendays = np.zeros([nbounds,nyears])   
	Precipdays = np.zeros([nbounds,nyears]) 

	for ibound in range(0,nbounds):
		if splittype == "maxspeed":
			fileadd = "MaxSpeeds_" + str(speedtspan) + "ts_"
		elif splittype == "speed":
			fileadd = "Speeds_"
		elif splittype == "day":
			fileadd = "Sizes_"

		if tbound1[ibound] < 0:
			tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))
		else:
			tboundtitle = str(int(tbound1[ibound])) + '-' + str(int(tbound2[ibound]))

		FileIn = 'DenDirSpd_Map_Ann_' + fileadd + tboundtitle + unit + '_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version +'.nc'

		FileIn = xray.open_dataset(DirIn + FileIn)

		if lats[-1] < lats[0]:
                        Dendays[ibound,:] = np.sum(FileIn['TDensity'].sel(lat=slice(MaxLatF,MinLatF),years=slice(anstartyr,anendyr)),axis=(1,2))
		else:
			Dendays[ibound,:] = np.sum(FileIn['TDensity'].sel(lat=slice(MinLatF,MaxLatF),years=slice(anstartyr,anendyr)),axis=(1,2))
	

	#	Dendays[ibound,:] = np.sum(FileIn['TDensityAnn'],axis=(1,2))
	#        Precipdays[ibound,:,:,:] = FileIn['TPrecipAnn']

	print Dendays
	toplot.append(Dendays)

print toplot

def plotnow(nlines,datain,yearsin,nsep,title,Datatitle,figtitlein):
	
	#Now plot plots
	wkres = Ngl.Resources()
	wkres.wkColorMap = "MPL_BrBG"
	wks_type = "eps"
	wks = Ngl.open_wks(wks_type,figtitlein,wkres)
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

	if (plotdensity):
		res.tiYAxisString = "frequency"
	else:
		res.tiYAxisString = "number of events"

	plotin = datain
	if title in ["Footprint"]:
		unit = "km~S1~2"
		convert = 1.0
	elif title in ["Total Precip"]:
		unit = "m~S1~3"
		convert = 1.0
	elif title in ["UniSpatSpan", "AvSpatSpan"]:
		unit = "gridboxes"
		convert = 1.0
	elif title in ["Timespan","TimeSpan"]:
		unit = "days"
		convert = 8.0
	else:
		unit = "unspecified"


	for iline in range(0,nlines):
		yearsplot = yearsin[iline]
		print yearsplot
#		nyearsin = len(yearsplot)
#		yearnums = range(0,nyearsin)
		A = np.array(yearsplot)
		regress = (np.zeros((5,nsep),np.float))

		for ibin in range(0,nsep):
			if ibin == 0:
				if (plotdensity):
					res.tiYAxisString = "frequency"
				else:
					res.tiYAxisString = "number of events"
			else:
				res.tiYAxisString = ""

			linreg = plotin[iline][ibin][:]
			print A.shape
			print linreg.shape
			regress[:,ibin] = stats.linregress(A,linreg)

			if ibin == nsep -1:
				res.tiMainString = '{:^80}'.format('          >' + '{:2.1g}'.format(tbound1[ibin]) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]) + "           ")
			else:
				res.tiMainString = '{:^80}'.format('{:2.1g}'.format(tbound1[ibin]) + '-' + '{:2.1g}'.format(tbound2[ibin]) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]))

			plot.append(Ngl.xy(wks,yearsplot,plotin[iline][ibin,:],res))

	panelres = Ngl.Resources()
	panelres.nglPanelLabelBar = True
	panelres.nglPanelYWhiteSpacePercent = 8.
	panelres.nglPanelXWhiteSpacePercent = 0.0

	panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
	panelres.nglPanelTop                      = 1.0
	panelres.nglPanelBottom                      = 0.00
	panelres.nglPanelLeft			= 0.0
	panelres.nglPanelRight			= 1.0
	panelres.nglPaperOrientation = "Portrait"
	panelres.nglScale = False
	panelres.nglMaximize = True
	
	#txres = Ngl.Resources()
	#txres.txFontHeightF = 0.012
	#Ngl.text_ndc(wks,'Annual timeseries of global number of events of various scales from' + Datatitle,0.5,0.94,txres)

	Ngl.panel(wks,plot,[nlines,float(nbounds)],panelres)

	print 'panelled'


# Now plot


# Time span
variable = "Timespan"
vartitle = "annualtrends"
if plotdensity:
	addtitle = '_dens'
else:
        addtitle = '_numb'

figtitle = "Paper_TRMM_ERA20C" + str(MinLatF) + "to" + str(MaxLatF) + "N_" + vartitle + addtitle

plotnow(2,toplot,toplotx,nbounds,variable,Data,FigDir + figtitle)

FileIn.close()
Ngl.end()









