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
import matplotlib.pyplot as plt


day1 = [0,1,2,5]
day2 = [1,2,5,100]
logiczonal = 1
Regions = ["Region0-5N+S","Region5-10N+S","Region10-20N+S","Region30-50N+S"]

speedtspan = 4

nregions = len(Regions)
ndayplots = len(day1)

#Version = 'standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
Version = 'Standard'
#Version = '5th_from48'

Data = "TRMM"

R = 6371000     # radius of Earth in m

if logiczonal == 1 and speedtspan == 0:
	speedbinedges = [-35.0,-32.5,-30.0,-27.5,-25.0,-22.5,-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,-0.1,0.1,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25.0,27.5,30.0,32.5,35.0]
else:
        speedbinedges = [-22.5,-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,-0.1,0.1,2.5,5,7.5,10,12.5,15,17.5,20,22.5]

histnbins = len(speedbinedges) - 1
speedbins = np.zeros(histnbins,np.float)
for ibin in range(0,histnbins-1):
	speedbins[ibin] = (speedbinedges[ibin] + speedbinedges[ibin+1]) / 2.0


sspanbinedges = [0,10,20,40,80,200,1000]
pspanbinedges = [0,800,1600,6400,20000,40000,200000]
regresst = (np.zeros((5,histnbins),np.float))

plotdensity = False

starttsteps = 0
endtsteps = 46752 # 16 years
anntsteps = 2920 # timesteps per year

minevent = 100000

FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

if Data == "TRMM":
	Fstartyr = 1998
	Fendyr = 2014
	mult = 3
	startyr = 1998
	endyr = 2014
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
elif Data == "ERAI":
	Fstartyr = 1980
	Fendyr = 2014

	startyr = 1980
	endyr = 2014
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
elif Data == "CESM":
	DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	FileI = 'DenDirSpd_Map_monthly_CESM_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

else:
	print("unexpected data type")
	exit()

nyears = endyr - startyr + 1


#Get lons and lats
iday = 0
FileInLats = '/home/disk/eos4/rachel/Obs/TRMM/3hrly/SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc'
print FileInLats
FileIn = xray.open_dataset(FileInLats)

lats = FileIn['Latitude'].values
lons = FileIn['Longitude'].values

nlats = len(lats)
nlons = len(lons)

plotspeed = np.zeros([ndayplots,histnbins])
plotzonalspeed = np.zeros([ndayplots,histnbins])

daylengths = []

for iregion in range(0,nregions):
	fig1 = plt.figure()

	region = Regions[iregion]

	for iday in range(0,ndayplots):
		daylengths.append(str(day1[iday]) + 'to' + str(day2[iday]))
		FileIn = region + '_Precip_Sizes_' + str(day1[iday]) + '-' + str(day2[iday]) + 'day_' + Data + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

		#Get lons and lats
		print DirIn + FileIn
		FileIn = xray.open_dataset(DirIn + FileIn)
		if speedtspan == 0:
			maxspeed = FileIn['xmaxtravel'].values
		else:
			maxspeed = FileIn['xmaxspeed_' + str(speedtspan) + 'ts'].values
		plotzonalspeed[iday,:],zonalspeed_bin = np.histogram(maxspeed,bins=speedbinedges,density=True)

		toplot = xray.DataArray(maxspeed,dims=['events'])
		plt.subplot(2,2,iday+1)
		if len(maxspeed) > 0:
			xray.plot.hist(toplot,bins=speedbinedges)

		plt.title(str(day1[iday]) + ' to ' + str(day2[iday]) + ' day events')

	if speedtspan == 0:
		plt.suptitle(region + ' events, maximum zonal speed')
		figtitle = region + "_histogram_maxzonal"
	else:
                plt.suptitle(region + ' events, maximum zonal speed over ' + str(speedtspan * mult) + 'hours')
                figtitle = region + "_histogram_maxzonal_" + str(speedtspan * mult) + '_hrs'


	fig1.tight_layout(pad=2)	# Put more white space around plots
	plt.savefig(FigDir + figtitle + ".eps",format='eps',dpi=1200., facecolor='w', edgecolor='w')
	#toplot = xray.DataArray(plotspeed,coords=[daylengths,speedbins],dims=['daylengths','bins'])


	#toplotzonal = xray.DataArray(plotzonalspeed,coords=[daylengths,speedbins],dims=['daylengths','bins'])

	#xray.plot.hist(toplot)

#	plt.show()

def get_bar(x,y,dx,ymin,bar_width_perc=0.6):
  dxp = (dx * bar_width_perc)/2.
  xbar = numpy.array([x-dxp,x+dxp,x+dxp,x-dxp,x-dxp])
  ybar = numpy.array([ ymin, ymin,    y,    y, ymin])
  return xbar,ybar

def plothist(datain,binsin,histnbins,title,Datatitle,figtitlein):

	nvalues = datain.shape[0]
        plotin = np.zeros([nvalues,histnbins])

        #Now plot log-log plots
        wkres = Ngl.Resources()
        wkres.wkColorMap = "MPL_BrBG"
        wks_type = "eps"
        wks = Ngl.open_wks(wks_type,figtitlein,wkres)
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

        res.tiYAxisString = "fraction"


def plotdensity(indensity,title):
	figtitle = Data + "_" + Version + '_' + title + "_" + str(histnbins) + 'bins_number_' + str(startyr) + '-' + str(endyr) + '.nc'

	regresst = (np.zeros((5,nbins),np.float))

	for ibin in range(0,nbins):
		linreg = indensity[:,ibin]
		regresst[:,ibin] = stats.linregress(A,linreg)

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

	res.tiYAxisString = "number"


	#print type(indensity[:,ibin])
	#print type(range(startyr,endyr+1))
	
	for ibin in range(0,nbins):
		res.tiMainString = 'Timespan ' +  str(ibin) +'; r2= ' + '{:5.3f}'.format(regresst[2,ibin] * regresst[2,ibin]) + '; p= ' + '{:5.3f}'.format(regresst[3,ibin])

		plot.append(Ngl.xy(wks,range(startyr,endyr+1),indensity[:,ibin],res))

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
	Ngl.text_ndc(wks,'Annual timeseries of global number of events of various timescales from ' + Data + " for " + title ,0.5,0.85,txres)

	Ngl.panel(wks,plot,[math.ceil(float(nbins)/3.0),3],panelres)


Ngl.end()









