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

Fstartyr = 1998
Fendyr = 2014

startyr = 2000
endyr = 2000
nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

starttsteps = 0
endtsteps = 46752 # 16 years
anntsteps = 2920 # timesteps per year

minevent = 100000

if Data == "TRMM":
	if Version == '6th_from6' or Version == '5th_from48' or Version == 'Standard':
        	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/' + Version + '/Precip/'
        	FileI1 = 'Precip_Sizes_' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
elif Data == "ERAI":
        DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/ERAI_output/' + Version + str(startyr) + '/Precip/'
        FileI1 = 'Precip_Sizes_ERAI' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'
elif Data == "CESM":
	DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '-' + str(Fendyr) + '/Precip/' 
	FileI1 = 'Precip_Sizes_CESM' + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + '.nc'

else:
	print("unexpected data type")
	exit()

FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

filetimespan = "3hrly"

print DirI + FileI1
datain = xray.open_dataset(DirI + FileI1)

print(datain.coords)

tspan=datain['timespan']
tmean = datain['tmean']
tstart = datain['tstart'].values
sspan = datain['uniquegridboxspan']
totalprecip = datain['totalprecip']
xmin1 = datain['xmin']

print totalprecip.shape

nevents = tspan.shape[0]

print nevents

curtime = tstart[starttsteps]
n = 0
iyear = 0

if startyr != Fstartyr:
       starttsteps = anntsteps * (startyr - Fstartyr)


for n in range(0,nevents):
	curtime = tstart[n]
	if curtime == starttsteps:
		mints[iyear] = n
		break

for iyear in range(0,nyears):
	for n in range(int(mints[iyear]),nevents):
		curtime = tstart[n]
		if curtime == starttsteps + (iyear + 1) * anntsteps:
			if iyear < nyears-1:
				mints[iyear+1] = n
			
			maxts[iyear] = n #swiched from n-1, because python in ranges doesn't include the last
			break

print mints
print maxts
print (maxts - mints)


# Sort events based on event size and bin data
for iyear in range(0,nyears):
	precipyear = totalprecip[mints[iyear]:maxts[iyear]]

	sortedprecip = np.sort(precipyear,axis=0,kind='quicksort')
	
	globalprecip = np.sum(precipyear)

	print globalprecip

	# Work backwards to find thresholds
	# first threshold: 75% of total:
	thresh75 = globalprecip * 0.75
	thresh50 = globalprecip * 0.5
	thresh25 = globalprecip * 0.25

	runningtotal = globalprecip
	# ievent is last event in this year
	neventsyr = len(precipyear)
	ievent = neventsyr - 1
	while runningtotal > thresh75:
		runningtotal = runningtotal - sortedprecip[ievent]
		ievent = ievent - 1
	nevents75 = ievent
        while runningtotal > thresh50:
                runningtotal = runningtotal - sortedprecip[ievent]
                ievent = ievent - 1
        nevents50 = ievent

        while runningtotal > thresh25:
                runningtotal = runningtotal - sortedprecip[ievent]
                ievent = ievent - 1
        nevents25 = ievent

	print nevents25
	print nevents50-nevents25
	print nevents75-nevents50
	print neventsyr-nevents75

	print sortedprecip[nevents25]
        print sortedprecip[nevents50]
        print sortedprecip[nevents75]
	print sortedprecip[neventsyr]

	print np.mean(sortedprecip[0:nevents25])
        print np.mean(sortedprecip[nevents25:nevents50])
        print np.mean(sortedprecip[nevents50:nevents75])
        print np.mean(sortedprecip[nevents75:nevents])


	exit()
	#Bin data 


	pstart = 0

	pspan_bin_spec = [pstart]
	pspan_bin_spec.append(100.0)
	pspan_bin_spec.append(500.0)
	pspan_bin_spec.append(5000.0)
	pspan_bin_spec.append(10000000.0)

histnbins = len(pspan_bin_spec)
print pspan_bin_spec

#Loop through years, binning data for each year
precipplot = np.zeros([histnbins])

iyear = 0
precipplot[:],precip_bin = np.histogram(totalprecip[mints[iyear]:maxts[iyear]],bins=pspan_bin_spec,density=plotdensity)


datain.close()
# Calculate total precip
globalprecip = np.sum(totalprecip[mints[iyear]:maxts[iyear]])




yearnums = range(0,nyears)
A = np.array(yearnums)
regresst = (np.zeros((5,histnbins),np.float))
regresss = (np.zeros((5,histnbins),np.float))
regressp = (np.zeros((5,histnbins),np.float))



for ibin in range(0,histnbins):
	linreg = tspanplot[:,ibin]
	regresst[:,ibin] = stats.linregress(A,linreg)
        linreg = sspanplot[:,ibin]
        regresss[:,ibin] = stats.linregress(A,linreg)
        linreg = precipplot[:,ibin]
        regressp[:,ibin] = stats.linregress(A,linreg)


#Now plot log-log plots
wkres = Ngl.Resources()
wkres.wkColorMap = "MPL_BrBG"
wks_type = "eps"
wks = Ngl.open_wks(wks_type,FigDir + figtitle1,wkres)
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

	
if (plotdensity):
	res.tiYAxisString = "frequency"
else:
	res.tiYAxisString = "number"


for ibin in range(0,histnbins):
	res.tiMainString = 'Timespan '  + str(tspan_bin[ibin]/8.0) + '-' + str(tspan_bin[ibin+1]/8.0) + 'days; r2= ' + '{:5.3f}'.format(regresst[2,ibin] * regresst[2,ibin]) + '; p= ' + '{:5.3f}'.format(regresst[3,ibin])

	plot.append(Ngl.xy(wks,range(startyr,endyr+1),tspanplot[:,ibin],res))

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
Ngl.text_ndc(wks,'Annual timeseries of global number of events of various timescales from ERAI',0.5,0.85,txres)
print 'about to panel'


Ngl.panel(wks,plot,[math.ceil(float(histnbins)/3.0),3],panelres)

print 'panelled'

#Now plot log-log plots
wkres = Ngl.Resources()
wks = Ngl.open_wks(wks_type,FigDir + figtitle2,wkres)
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

if (plotdensity):
        res.tiYAxisString = "frequency"
else:
        res.tiYAxisString = "number"

for ibin in range(0,histnbins):
	res.tiMainString = 'Uni Sspan annual, bin '  + str(tspan_bins_plot[ibin]) + '; r2= ' + '{:5.3f}'.format(regresss[2,ibin] * regresss[2,ibin]) + '; p= ' + '{:5.3f}'.format(regresss[3,ibin])
        plot.append(Ngl.xy(wks,range(startyr,endyr+1),sspanplot[:,ibin],res))

panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.98
panelres.nglPanelBottom                      = 0.02
panelres.nglPaperOrientation = "Auto"

print 'about to panel'


Ngl.panel(wks,plot,[math.ceil(float(histnbins)/3.0),3],panelres)

print 'panelled'



#Now plot log-log plots
wkres = Ngl.Resources()
wks = Ngl.open_wks(wks_type,FigDir + figtitle3,wkres)
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

if (plotdensity):
        res.tiYAxisString = "frequency"
else:
        res.tiYAxisString = "number"

for ibin in range(0,histnbins):
        res.tiMainString = 'Precip annual for bin '  + str(precip_bins_plot[ibin]) + '; r2= ' + '{:5.3f}'.format(regressp[2,ibin] * regressp[2,ibin]) + '; p= ' + '{:5.3f}'.format(regressp[3,ibin])
        plot.append(Ngl.xy(wks,range(startyr,endyr+1),precipplot[:,ibin],res))

panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.98
panelres.nglPanelBottom                      = 0.02
panelres.nglPaperOrientation = "Auto"

print 'about to panel'


Ngl.panel(wks,plot,[math.ceil(float(histnbins)/3.0),3],panelres)

print 'panelled'




Ngl.end()









