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



Basins = np.array(["All","Pacific","Atlantic"])
Regions = np.array(["Mid","Tr"])
nbasins = Basins.shape[0]
nregions = Regions.shape[0]
for ibasin in range(0,nbasins):
	for iregion in range(0,nregions):
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
			SizesIn = Basins[ibasin] + Regions[iregion] + "_Precip_Sizes_1-2day" 
		elif Data == "ERAI":
			DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
		elif Data == "CESM":
			startyr = 1990 # Don't change - tied to file names!
			endyr = 2014
			DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(startyr) + '-' + str(endyr) + '/Precip/'


		FileIn1to2 = Basins[ibasin] + Regions[iregion] + "_Precip_Sizes_1-2day_" + Data + str(startyr) + "-" + str(endyr) + "_" + Version + ".nc"

                FileIn2to5 = Basins[ibasin] + Regions[iregion] + "_Precip_Sizes_2-5day_" + Data + str(startyr) + "-" + str(endyr) + "_" + Version + ".nc"

		#Get lons and lats
		print DirIn + FileIn1to2
		FileIn1to2 = xray.open_dataset(DirIn + FileIn1to2)
		FileIn2to5 = xray.open_dataset(DirIn + FileIn2to5)

		tstart12 = FileIn1to2["tstart"].values
                timespan12 = FileIn1to2["timespan"].values
                totalprecip12 = FileIn1to2["totalprecip"].values
                gridboxspan12 = FileIn1to2["gridboxspan"].values

		tstart25 = FileIn2to5["tstart"].values	
                timespan25 = FileIn2to5["timespan"].values
		totalprecip25 = FileIn2to5["totalprecip"].values
		gridboxspan25 = FileIn2to5["gridboxspan"].values

		nevents = [tstart12.shape[0],tstart25.shape[0]]

		mints = np.zeros([2,nyears],np.double)
		maxts = np.zeros([2,nyears],np.double)

		plotden = np.zeros([2,nyears],np.int)
		plotprecip = np.zeros([2,nyears],np.double)
		ploteventprecip = np.zeros([2,nyears],np.double)
		plotgridprecip = np.zeros([2,nyears],np.double)
                plottimespanevent = np.zeros([2,nyears],np.double)
                plotsizestep = np.zeros([2,nyears],np.double)


		iyear = 0
		for iyear in range(0,nyears):
			for n in range(int(mints[0,iyear]),nevents[0]):
				curtime12 = tstart12[n]
				if curtime12 >= (iyear + 1) * anntsteps:
					if iyear < nyears-1:
						mints[0,iyear+1] = n
					maxts[0,iyear] = n #swiched from n-1, because python in ranges doesn't include the last
					break
			for n in range(int(mints[1,iyear]),nevents[1]):
				curtime25 = tstart25[n]
				if curtime25 >= (iyear + 1) * anntsteps:
					if iyear < nyears-1:
						mints[1,iyear+1] = n
					maxts[1,iyear] = n #swiched from n-1, because python in ranges doesn't include the last
					break

		print maxts[0]-mints[0]
		print maxts[1]-mints[1]


		for iyear in range(0,nyears):
			iday = 0
			plotden[iday,iyear] = maxts[iday,iyear] - mints[iday,iyear]

			plotprecip[iday,iyear] = np.sum(totalprecip12[mints[iday,iyear]:maxts[iday,iyear]])
			ploteventprecip[iday,iyear] = plotprecip[iday,iyear]/plotden[iday,iyear]
			plotgridprecip[iday,iyear] = plotprecip[iday,iyear]/(np.sum(gridboxspan12[mints[iday,iyear]:maxts[iday,iyear]]))
			plottimespanevent[iday,iyear] = np.sum(timespan12[mints[iday,iyear]:maxts[iday,iyear]])/plotden[iday,iyear]
			plotsizestep[iday,iyear] = np.sum(gridboxspan12[mints[iday,iyear]:maxts[iday,iyear]])/np.sum(timespan12[mints[iday,iyear]:maxts[iday,iyear]])


			iday = 1
			plotden[iday,iyear] = maxts[iday,iyear] - mints[iday,iyear]

			plotprecip[iday,iyear] = np.sum(totalprecip25[mints[iday,iyear]:maxts[iday,iyear]])
			ploteventprecip[iday,iyear] = plotprecip[iday,iyear]/plotden[iday,iyear]
			plotgridprecip[iday,iyear] = plotprecip[iday,iyear]/(np.sum(gridboxspan25[mints[iday,iyear]:maxts[iday,iyear]]))
                        plottimespanevent[iday,iyear] = np.sum(timespan25[mints[iday,iyear]:maxts[iday,iyear]])/plotden[iday,iyear]
                        plotsizestep[iday,iyear] = np.sum(gridboxspan25[mints[iday,iyear]:maxts[iday,iyear]])/np.sum(timespan25[mints[iday,iyear]:maxts[iday,iyear]])


		titlerows = ['1-2days','2-5days']
		titlecols = ['number of events','total precip mm', 'precip per event mm/event', 'precip per gridbox, mm/gridbox','timespan per event','size per timestep']

	#       plotdensity(np.array([Den0E,Den1E,Den2E,Den3E,Den4E,Den0W,Den1W,Den2W,Den3W,Den4W,Den0E + Den0W,Den1W + Den1E,Den2W+Den2E,Den3W+Den3E,Den4W+Den4E]),basin + '_' + str(day1) + '-' + str(day2) + 'days',titlerows,titlecols)
		print Regions[iregion]

		plotdensity(np.array([plotden[0,:],plotprecip[0,:],ploteventprecip[0,:],plotgridprecip[0,:],plottimespanevent[0,:],plotsizestep[0,:],plotden[1,:],plotprecip[1,:],ploteventprecip[1,:],plotgridprecip[1,:],plottimespanevent[1,:],plotsizestep[1,:]]),A,Basins[ibasin] + '_' + Regions[iregion] + '_',titlerows,titlecols)



Ngl.end()




