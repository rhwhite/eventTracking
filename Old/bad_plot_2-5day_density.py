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

#Version = 'standard'
#Version = '5th_nanto25'
#Version = '5th_nantozero'
#Version = '7thresh'
Version = 'Standard'
#Version = '5th_from48'

Data = "TRMM"

day1 = 1
day2 = 2

centering = "centre"
plotdensity = False

starttsteps = 0
endtsteps = 46752 # 16 years
anntsteps = 2920 # timesteps per year

minevent = 100000

if Data == "TRMM":
	Fstartyr = 1998
	Fendyr = 2014

	startyr = 1998
	endyr = 2014
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
		FileIP = 'DenDirSpd_Map_' + str(day1) + '-' + str(day2) + 'Day_' + centering + '_' + Data + '_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
                FileI = 'DenDirSpd_Map_noP' + str(day1) + '-' + str(day2) + 'Day_' + centering + '_' + Data + '_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

elif Data == "ERAI":
	Fstartyr = 1980
	Fendyr = 2014

	startyr = 1980
	endyr = 2014
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
	FileIP = 'DenDirSpd_Map_monthly_ERAI_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
        FileI = 'DenDirSpd_Map_monthly_ERAI_noP' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

elif Data == "CESM":
	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	FileI = 'DenDirSpd_Map_monthly_CESM_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'

else:
	print("unexpected data type")
	exit()

nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)


FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

filetimespan = "3hrly"

print DirI + FileI
datain = xray.open_dataset(DirI + FileI)
datainP = xray.open_dataset(DirI + FileIP)

#print(datain.coords)

lats = datain['Latitude']
lons = datain['Longitude']

nlons = len(lons)

years = datain['years']
seasons = datain['Seasons']

# Read in land-sea mask? How do we get one for TRMM? Generic one, regridded from other, with min 50% land?

# Define lat-lon mask for:
# Pacific: 120E - 280E
# Atlantic: 280E - 360E and 0 - 20E 
# Indian: 35E - 120E
Psrt = 120.0
if lons[0] < 0:
	Pend = -80.0      #280.0
else:
	Pend = 280.0

if lons[0] < 0:
	Asrt = -80.0	#280
else:
	Asrt = 280.0
Aend = 30.0

Isrt = 30.0
Iend = 120.0
# Find start lons indices for these
diffs = abs(lons[1] - lons[0])/2.0

def SplitBasin(srt,end,density):
	Starti = np.intersect1d(np.where(lons >= srt - diffs),np.where(lons < srt + diffs),False)[0]
	Endi = np.intersect1d(np.where(lons > end - diffs),np.where(lons <= end + diffs),False)[0]
	ndims = len(density.shape)
	
	if Starti > Endi:
		outden = np.concatenate([density[...,Starti:nlons],density[...,0:Endi]],axis = ndims-1)
	else:
		outden = density[...,Starti:Endi]

	return outden


# Calculate trend analysis

yearnums = range(0,nyears-1)
A = np.array(yearnums)

def plotdensity(indensity,title,titlerows,titlecols):
	figtitle = Data + "_" + Version + '_' + title + "_" + str(startyr) + '-' + str(endyr) + '_AnnualReg'

	nrows = len(titlerows)
	ncols = len(titlecols)
	nplots = nrows * ncols
        titles = []
        for trow in titlerows:
                for tcol in titlecols:
                        titles.append(trow + ' ' + tcol)

	regresst = (np.zeros((5,nplots),np.float))
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

	#print type(indensity[:,ibin])
	#print type(range(startyr,endyr+1))

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


# Call SplitBasin for Pacific:

# Get relevant data:
densityE = datain['EDensity'].values
densityW = datain['WDensity'].values

TotalPrecipE = datainP['ETotalPrecip'].values
TotalPrecipW = datainP['WTotalPrecip'].values

#AverageSizeE = 0.001 * 0.001 * datain['ESize'].values / datain['ETSpan'] 
#AverageSizeW = 0.001 * 0.001 * datain['WSize'].values / datain['WTSpan']

AverageSizeE = datain['ESize'].values / datain['ETSpan'].values
AverageSizeW = datain['WSize'].values / datain['WTSpan'].values

AvgPrecipE = datainP['EPrecip'].values / datain['ETSpan'].values
AvgPrecipW = datainP['WPrecip'].values / datain['WTSpan'].values

TimeE = datain['ETSpan'].values
TimeW = datain['WTSpan'].values

#density2E[density2E > 5E35] = 0.0
#density2W[density2W > 5E35] = 0.0

nyears = densityE.shape[0]
nseas = densityE.shape[1]

ndims = len(densityE.shape)

nevents3 = (SplitBasin(lons[0],lons[nlons-1],densityE) != 0).sum(3).sum(2).sum(1)
nevents4 = (SplitBasin(lons[0],lons[nlons-1],densityW) != 0).sum(3).sum(2).sum(1)

print nevents3 
print nevents4

for basin in ("All","Atlantic","Pacific","Indian"):

	if basin == "All":
		srt = lons[0]
		end = lons[nlons-1]
	if basin == "Atlantic":
		srt = Asrt
		end = Aend
	elif basin == "Pacific":
		srt = Psrt
		end = Pend
	elif basin == "Indian":
		srt = Isrt
		end = Iend

	DenE = np.nansum(SplitBasin(srt,end,densityE),axis = (1,ndims-2,ndims-1))
	AvPE = np.nansum(SplitBasin(srt,end,AvgPrecipE),axis = (1,ndims-2,ndims-1))/((SplitBasin(srt,end,densityE) != 0).sum(3).sum(2).sum(1))
	AvSE = np.nansum(SplitBasin(srt,end,AverageSizeE),axis = (1,ndims-2,ndims-1))/((SplitBasin(srt,end,densityE) != 0).sum(3).sum(2).sum(1))
	TPE = np.nansum(SplitBasin(srt,end,TotalPrecipE),axis = (1,ndims-2,ndims-1))

	DenW = np.nansum(SplitBasin(srt,end,densityW),axis = (1,ndims-2,ndims-1))
	AvPW = np.nansum(SplitBasin(srt,end,AvgPrecipW),axis = (1,ndims-2,ndims-1))/((SplitBasin(srt,end,densityW) != 0).sum(3).sum(2).sum(1))
	AvSW = np.nansum(SplitBasin(srt,end,AverageSizeW),axis = (1,ndims-2,ndims-1))/((SplitBasin(srt,end,densityW) != 0).sum(3).sum(2).sum(1))
	TPW = np.nansum(SplitBasin(srt,end,TotalPrecipW),axis = (1,ndims-2,ndims-1))

	DenT = DenE + DenW
	AvPT = AvPE + AvPW
	AvST = AvSE + AvSW
	TPT = TPE + TPW

	#Plot each basin separately
	print basin

	titlerows = ['Easterly','Westerly','All']
	titlecols = ['number of events','avg event precip mm/hr', 'average event size, gridboxes', 'total precip, m^3']

#	plotdensity(np.array([Den0E,Den1E,Den2E,Den3E,Den4E,Den0W,Den1W,Den2W,Den3W,Den4W,Den0E + Den0W,Den1W + Den1E,Den2W+Den2E,Den3W+Den3E,Den4W+Den4E]),basin + '_' + str(day1) + '-' + str(day2) + 'days',titlerows,titlecols)

        plotdensity(np.array([DenE,AvPE,AvSE,TPE,DenW,AvPW,AvSW,TPW,DenT,AvPT,AvST,TPT]),basin + '_' + str(day1) + '-' + str(day2) + 'days_3',titlerows,titlecols)
datain.close()


Ngl.end()









