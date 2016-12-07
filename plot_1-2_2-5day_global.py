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
from scipy import stats
import math

day1 = 1
day2 = 2

# Time period for analysis
astartyr = 1998
aendyr = 2014
nyears = aendyr - astartyr + 1
print nyears

Data = "TRMM"
mapping = 'centre'
Version = 'Standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
#Version = '6th_from6'
#Version = '5th_from48'

#anntyears
if Data == "TRMM":
	anntsteps = 8 * 365

day1 = 1
day2 = 2


plotmin1 = 0
plotmax1 = 3
plotspace1 = 0.25

plotmin2 = -0.1
plotmax2 = 0.1
plotspace2 = 0.02

plotmin3 = 0
plotmax3 = 3
plotspace3 = 0.25

plotmin4 = -0.1
plotmax4 = 0.1
plotspace4 = 0.02

MinLonF = 0
MaxLonF = 360
MinLatF = -45
MaxLatF = 45

yearnums = range(astartyr,aendyr+1)
A = np.array(yearnums)

def regressmaps(m,c,r,p,stderr,A,linreg,nlats,nlons):
# Calculate linear trend
	for ilat in range(0,nlats):
		for ilon in range(0,nlons):
			m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(A,linreg[:,ilat,ilon])


def plotdensity(indensity,A,title,titlerows,titlecols):
        figtitle = title + Data + "_" + Version + "_" + str(startyr) + '-' + str(endyr) + '_AnnualReg'

        nrows = len(titlerows)
        ncols = len(titlecols)
        nplots = nrows * ncols
        titles = []
        for trow in titlerows:
                for tcol in titlecols:
                        titles.append(trow + ' ' + tcol)

        regresst = (np.zeros((5,nplots),np.float))
	print indensity.shape
        if nplots > 1:
                for iplot in range(0,nplots):
                        linreg = indensity[iplot,:]
			regresst[:,iplot] = stats.linregress(A,linreg)
        else:
                regresst[:,0] = stats.linregress(A,indensity)

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
        res.tiMainFontHeightF = "0.015"
        plot = []

        if nplots > 1:
                for iplot in range(0,nplots):
                        res.tiMainString = titles[iplot] +'; r2= ' + '{:5.3f}'.format(regresst[2,iplot] * regresst[2,iplot]) + '; p= ' + '{:5.3f}'.format(regresst[3,iplot])
                        res.tiYAxisString = titlecols[iplot - (iplot//ncols)*ncols]

                        plot.append(Ngl.xy(wks,range(startyr,endyr+1),indensity[iplot,:],res))
        else:
                res.tiMainString = titles +'; r2= ' + '{:5.3f}'.format(regresst[2,0] * regresst[2,0]) + '; p= ' + '{:5.3f}'.format(regresst[3,0])

                plot.append(Ngl.xy(wks,range(startyr,endyr),indensity,res))

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
        Ngl.text_ndc(wks,'Annual timeseries for ' + title ,0.5,0.85,txres)

        Ngl.panel(wks,plot,[math.ceil(float(nplots)/ncols),ncols],panelres)



FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
if Data == "TRMM":
	startyr = 1998 # Don't change - tied to file names!
	endyr = 2014
	if Version == 'Standard':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Standard/Precip/'
	elif Version == '5th_nanto25':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
	elif Version == '5th_nantozero':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
	elif Version == '7thresh':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
	elif Version == '6th_from6':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
	elif Version == '5th_from48':
		DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'

	# When this is updated need to shift seasons back!!!
elif Data == "ERAI":
	DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
elif Data == "CESM":
	startyr = 1990 # Don't change - tied to file names!
	endyr = 2014
	DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'


FileIn = "Precip_Sizes_" + str(day1) + '-' + str(day2) + "day_" + str(startyr) + "-" + str(endyr) + "_" + Version + ".nc"

#Get lons and lats
print DirIn + FileIn
FileInday = xray.open_dataset(DirIn + FileIn)

tstart = FileInday["tstart"].values
timespan = FileInday["timespan"].values
totalprecip = FileInday["totalprecip"].values
gridboxspan = FileInday["gridboxspan"].values

nevents = tstart.shape[0]

mints = np.zeros([nyears],np.double)
maxts = np.zeros([nyears],np.double)

plotden = np.zeros([nyears],np.int)
plotprecip = np.zeros([nyears],np.double)
plotgridboxspan = np.zeros([nyears],np.double)
plottimespan = np.zeros([nyears],np.double)

plotdenevent = np.zeros([nyears],np.double)
plotprecipevent = np.zeros([nyears],np.double)
plotgridevent = np.zeros([nyears],np.double)
plottimeevent = np.zeros([nyears],np.double)

plotdengrid = np.zeros([nyears],np.double)
plotprecipgrid = np.zeros([nyears],np.double)
plotgridgrid = np.zeros([nyears],np.double)
plottimegrid = np.zeros([nyears],np.double)

iyear = 0
for iyear in range(0,nyears):
	for n in range(int(mints[iyear]),nevents):
		curtime = tstart[n]
		if curtime >= (iyear + 1) * anntsteps:
			if iyear < nyears-1:
				mints[iyear+1] = n
			maxts[iyear] = n #swiched from n-1, because python in ranges doesn't include the last
			break

print maxts-mints


for iyear in range(0,nyears):
	iday = 0

	plotden[iyear] = maxts[iyear] - mints[iyear]
	plotprecip[iyear] = np.sum(totalprecip[mints[iyear]:maxts[iyear]])
	plotgridboxspan[iyear] = np.sum(gridboxspan[mints[iyear]:maxts[iyear]])
	plottimespan[iyear] = np.sum(timespan[mints[iyear]:maxts[iyear]])

	plotdenevent[iyear] = plotden[iyear]/plotden[iyear]
	plotprecipevent[iyear] = plotprecip[iyear]/plotden[iyear]
	plotgridevent[iyear] = plotgridboxspan[iyear]/plotden[iyear]
	plottimeevent[iyear] = plottimespan[iyear]/plotden[iyear]	

	plotdengrid[iyear] = np.nan
	plotprecipgrid[iyear] = plotprecip[iyear]/plotgridboxspan[iyear]
	plotgridgrid[iyear] = plotgridboxspan[iyear]/plotgridboxspan[iyear]
	plottimegrid[iyear] = np.nan

	#plotsizestep[iyear] = np.sum(gridboxspan[mints[iyear]:maxts[iyear]])/np.sum(timespan[mints[iyear]:maxts[iyear]])

titlerows = ['in total','per event','per gridbox']
titlecols = ['number of events','total precip mm', 'total gridboxes', 'total timespan']

#       plotdensity(np.array([Den0E,Den1E,Den2E,Den3E,Den4E,Den0W,Den1W,Den2W,Den3W,Den4W,Den0E + Den0W,Den1W + Den1E,Den2W+Den2E,Den3W+Den3E,Den4W+Den4E]),basin + '_' + str(day1) + '-' + str(day2) + 'days',titlerows,titlecols)

plotdensity(np.array([plotden,plotprecip,plotgridboxspan,plottimespan,plotdenevent,plotprecipevent,plotgridevent,plottimeevent,plotdengrid,plotprecipgrid,plotgridgrid,plottimegrid]),A,'Global_' + str(day1) + '-' + str(day2),titlerows,titlecols)



Ngl.end()




