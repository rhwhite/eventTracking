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

Data = "CESM"

day1 = [0,1,2,5]
day2 = [1,2,5,100]

ndayplots = len(day1)

mapping = 'centre'
Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

MinLonF = 0
MaxLonF = 360
MinLatF = -20
MaxLatF = 20

sumlats = 1
sumlons = 1

test = 0

regresst = (np.zeros((5,ndayplots),np.float))


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
	endyr = 2014

	anstartyr = 1998
	anendyr = 2015

elif Data == "ERAI":
	Fstartyr = 1980
        Fendyr = 2014
        startyr = 1980
        endyr = 2014

	anstartyr = 1980
	anendyr = 2015

if test == 1:
	endyr = startyr
nyears = anendyr - anstartyr +1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

starttsteps = 0
anntsteps = 2920 # timesteps per year

minevent = 100000

if Data == "TRMM":
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
        	FileI1 = 'Precip_Sizes_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
elif Data == "ERAI":
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
        FileI1 = 'Precip_Sizes_ERAI' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
elif Data == "CESM":
	DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	FileI1 = 'Precip_Sizes_CESM' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

else:
	print("unexpected data type")
	exit()

FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

filetimespan = "3hrly"

# Get lats and lons
iday = 0
File = 'DenDirSpd_Map_Sizes_' + str(day1[iday]) + '-' + str(day2[iday]) + 'day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
print DirIn + File
FileIn = xray.open_dataset(DirIn + File)

lats = FileIn['lat']
lons = FileIn['lon']
years = FileIn['years'].sel(years = slice(anstartyr,anendyr))
print years
if lats[-1] < lats[0]:
	exit("latitudes are not south to north!")

nlats = len(lats)
nlons = len(lons)
nyears = len(years)
print nyears

Dendays = np.zeros([ndayplots,nyears])    # +1 for total sum as first plot
Precipdays = np.zeros([ndayplots,nyears]) # + 1 for total precip as first plot

print(FileIn['TDensityAnn'])


for iday in range(0,ndayplots):
        FileIn = 'DenDirSpd_Map_Sizes_' + str(day1[iday]) + '-' + str(day2[iday]) + 'day_' + mapping + '_' + Data + "_" + str(startyr) + '-' + str(endyr) + '_' + Version + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'

        #Get lons and lats
        print DirIn + FileIn
        FileIn = xray.open_dataset(DirIn + FileIn)

        Dendays[iday,:] = np.sum(FileIn['TDensityAnn'].sel(lat=slice(MinLatF,MaxLatF),years=slice(anstartyr,anendyr)),axis=(1,2))
#	Dendays[iday,:] = np.sum(FileIn['TDensityAnn'],axis=(1,2))
#        Precipdays[iday,:,:,:] = FileIn['TPrecipAnn']

print Dendays

def plotnow(nlines,datain,nyearsin,nsep,title,Datatitle,figtitlein):
	
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
	print plotin.shape
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

	yearnums = range(0,nyearsin)
	A = np.array(yearnums)
	regress = (np.zeros((5,nsep),np.float))

	for ibin in range(0,nsep):
		if ibin == 0:
			if (plotdensity):
				res.tiYAxisString = "frequency"
			else:
				res.tiYAxisString = "number of events"
		else:
			res.tiYAxisString = ""

		linreg = plotin[ibin][:]
		regress[:,ibin] = stats.linregress(A,linreg)

		if ibin == nsep -1:
			res.tiMainString = '{:^80}'.format('          >' + '{:2.1g}'.format(day1[ibin]) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]) + "           ")
		else:
			res.tiMainString = '{:^80}'.format('{:2.1g}'.format(day1[ibin]) + '-' + '{:2.1g}'.format(day2[ibin]) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]))

		plot.append(Ngl.xy(wks,range(startyr,endyr+1),plotin[ibin,:],res))

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

	Ngl.panel(wks,plot,[nlines,float(ndayplots)],panelres)

	print 'panelled'


# Now plot


# Time span
variable = "Timespan"
vartitle = "annualtrends"
if plotdensity:
	addtitle = '_dens'
else:
        addtitle = '_numb'

figtitle = "Paper_frommap_" + str(MinLatF) + "to" + str(MaxLatF) + "N_" + Data + "_" + Version + '_' + vartitle + addtitle

plotnow(1,Dendays,nyears,ndayplots,variable,Data,FigDir + figtitle)

FileIn.close()
Ngl.end()









