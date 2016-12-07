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
day2 = [1,2,5,10]
logiczonal = 1
Regions = ["Mid","Tr"]

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

histnbins = 17 # number of bins

speedbinedges = [-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,-0.1,0.1,2.5,5,7.5,10,12.5,15,17.5,20]
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
	if region == "Tr":
		regionname = "Tropical (15S-15N)"
	elif region == "Mid":
		regionname = "Mid-latitudes (15-45N/S)"

	for iday in range(0,ndayplots):
		daylengths.append(str(day1[iday]) + 'to' + str(day2[iday]))
		FileIn = 'All' + region + '_Precip_Sizes_' + str(day1[iday]) + '-' + str(day2[iday]) + 'day_' + Data + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

		#Get lons and lats
		print DirIn + FileIn
		FileIn = xray.open_dataset(DirIn + FileIn)

		ystart = FileIn['ycenterstart']
		xstart = FileIn['xcenterstart'] 
		yend = FileIn['ycenterend']
		xend = FileIn['xcenterend'] 
		tspan = FileIn['timespan'].values

		# Calculate distance using the Haversine formula, so we can do this using numpy

		lats1 = np.radians(lats[ystart[:]])
		lats2 = np.radians(lats[yend[:]])
		lons1 = np.radians(lons[xstart[:]])
		lons2 = np.radians(lons[xend[:]])

		a = (np.power(np.sin((lats2-lats1)/2),2) + np.cos(lats2) * np.cos(lats1) * np.power(np.sin((lons2-lons1)/2),2))
		c = 2.0 * np.arctan(np.sqrt(a),np.sqrt(1-a))
		distance = R * c

		latsmean = (lats1 + lats2) / 2.0

		az = (np.cos(latsmean) * np.cos(latsmean) * np.power(np.sin((lons2-lons1)/2),2))
		cz = 2.0 * np.arctan(np.sqrt(az),np.sqrt(1-az))

		distancez = R * cz

		pi = 3.14
		# Calculate speed
		angle = np.arctan2((lats2-lats1),(lons2-lons1))
		ones = np.zeros(lons1.shape)
		ones[...] = 1.0
		negones = ones * -1.0

		direction = np.where(lons2>=lons1,ones,negones)	# True where 
		speed = direction * distance / (tspan*3.0*60.0*60.0)
		zonalspeed = direction * distancez/(tspan*3.0*60.0*60.0)

	#	speed = np.where(tspan <1,0.0,speed)
	#        zonalspeed = np.where(tspan <1,0.0,zonalspeed)

		plotspeed[iday,:], speed_bin = np.histogram(speed,bins=speedbinedges,density=True)
		plotzonalspeed[iday,:],zonalspeed_bin = np.histogram(zonalspeed,bins=speedbinedges,density=True)

		if logiczonal == 1:
			toplot = xray.DataArray(zonalspeed,dims=['events'])
		else:
			toplot = xray.DataArray(speed,dims=['events'])

		plt.subplot(2,2,iday+1)
		xray.plot.hist(toplot,bins=speedbinedges)

		plt.title(str(day1[iday]) + ' to ' + str(day2[iday]) + ' day events')

	if logiczonal ==1:
		plt.suptitle(regionname + ' events, Zonal speed')
		figtitle = "All_" + region + "_histogram_zonal"
	else:
		plt.suptitle(regionname + ' events, Total speed')
		figtitle = "All_" + region + "_histogram_speed"

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


def SplitBasin(srt,end,density):
	Starti = np.intersect1d(np.where(lons >= srt - diffs),np.where(lons < srt + diffs),False)[0]
	Endi = np.intersect1d(np.where(lons > end - diffs),np.where(lons <= end + diffs),False)[0]

	if Starti > Endi:
		outden = np.concatenate([density[:,:,:,:,Starti:nlons],density[:,:,:,:,0:Endi]],axis = 4)
	else:
		outden = density[:,:,:,:,Starti:Endi]

	return outden

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









