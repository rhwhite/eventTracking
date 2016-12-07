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

Data = "ERAI"

histnbins = 6 # number of bins
tspanbinedges = [0,8,16,24,40,72,100]
sspanbinedges = [0,10,20,40,80,200,1000]
pspanbinedges = [0,800,1600,6400,20000,40000,200000]
regresst = (np.zeros((5,histnbins),np.float))

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
		FileIP = 'DenDirSpd_Map_monthly_' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'
                FileI = 'DenDirSpd_Map_monthly_noP' + str(startyr) + '-' + str(endyr) + '_' + Version + '.nc'


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

density = datain['DensityMap'].values
densityEasterly = datain['WestwardDensityMap']
lats = datain['Latitude']
lons = datain['Longitude']

nlons = len(lons)

years = datain['years']
months = datain['months']

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
Aend = 20.0

Isrt = 35.0
Iend = 120.0
# Find start lons indices for these
diffs = abs(lons[1] - lons[0])/2.0

print Psrt,Pend,Asrt,Aend,Isrt,Iend

def SplitBasin(srt,end,density):
	Starti = np.intersect1d(np.where(lons >= srt - diffs),np.where(lons < srt + diffs),False)[0]
	Endi = np.intersect1d(np.where(lons > end - diffs),np.where(lons <= end + diffs),False)[0]

	if Starti > Endi:
		outden = np.concatenate([density[:,:,:,:,Starti:nlons],density[:,:,:,:,0:Endi]],axis = 4)
	else:
		outden = density[:,:,:,:,Starti:Endi]

	return outden

# Call SplitBasin for Pacific:
PacDen = SplitBasin(Psrt,Pend,density)
AtlDen = SplitBasin(Asrt,Aend,density)
IndDen = SplitBasin(Isrt,Iend,density)

AllDen = np.sum(density,axis = (3,4))
AtlDen = np.sum(AtlDen,axis = (3,4))
PacDen = np.sum(PacDen,axis = (3,4))
IndDen = np.sum(IndDen,axis = (3,4))


# Create annual means
AllDen = np.sum(AllDen,axis = 1)
AtlDen = np.sum(AtlDen,axis = 1)
PacDen = np.sum(PacDen,axis = 1)
IndDen = np.sum(IndDen,axis = 1)

datain.close()
# Calculate trend analysis

nbins = AtlDen.shape[1]

yearnums = range(0,nyears)
A = np.array(yearnums)

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


#Plot each basin separately
print "Atlantic"
print type(AtlDen)
print type(PacDen)
plotdensity(AtlDen,"Atlantic")
print "Pacific"
print type(PacDen)
plotdensity(PacDen,"Pacific")
print "Indian"
plotdensity(IndDen,"Indian")
print "All"
plotdensity(AllDen,"All")


Ngl.end()









