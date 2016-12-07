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
tspanbinedges = [0,1.5,2.5,3.5,4.5,5.5,8.5]#6.5,7.5,8.5]
tspanbinedges2 = [0,8.5,16.5,40.5,48.5,56.5,96.5]#64.5,72.5,80.5]
usspanbinedges = [0,20,40,60,100,200,10000]#,400,1000,10000]
sspanbinedges = [0,40,80,120,400,800,40000]#1600,4000,40000]
avsspanbinedges = [0,10,20,30,50,70,150]#150,350,10000]

pspanbinedges = [0,800,1600,2400,3200,6400,20000,40000,200000]
regresst = (np.zeros((5,histnbins),np.float))


if Data == "CESM":
	Fstartyr = 1990
	Fendyr = 2014

	startyr = 1990
	endyr = 2011 
elif Data == "TRMM":
	Fstartyr = 1998
	Fendyr = 2014
	startyr = 1998
	endyr = 2014
elif Data == "ERAI":
	Fstartyr = 1980
        Fendyr = 2014
        startyr = 1980
        endyr = 2014


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
usspan = datain['uniquegridboxspan']
sspan = datain['gridboxspan']
totalprecip = datain['totalprecip']
xmin1 = datain['xmin']

avsspan = np.divide(sspan,tspan)

nevents = tspan.shape[0]

print nevents

curtime = tstart[starttsteps]
n = 0
iyear = 0
for n in range(0,nevents):
	curtime = tstart[n]
	if curtime == starttsteps:
		mints[iyear] = n
		break

for iyear in range(0,nyears):
	for n in range(int(mints[iyear]),nevents):
		curtime = tstart[n]
		if curtime == (iyear + 1) * anntsteps:
			if iyear < nyears-1:
				mints[iyear+1] = n
			
			maxts[iyear] = n #swiched from n-1, because python in ranges doesn't include the last
			break

print mints
print maxts
print (mints - maxts)

#Bin data 

if len(tspanbinedges) > 1:
	tspan_bin_spec = tspanbinedges
else:
	tstart = 0.9
	tmult = 30/float(histnbins)
	tpower = 1 + 2.5/float(histnbins)
	trange = int(100/histnbins)
	tmax = 500

	tspan_bin_spec = [tstart]
	for binidx in range(1,histnbins):
		tspan_bin_spec.append(tspan_bin_spec[binidx-1] + tmult * ((math.pow(tpower,binidx))))

	tspan_bin_spec.append(tmax)

if len(sspanbinedges) > 1:
	sspan_bin_spec = sspanbinedges
else:
	sstart = 0.9
	spower = 1 + 10.0/float(histnbins)
	smult = 75.0/float(histnbins)
	smax = 200000

	sspan_bin_spec = [sstart]
	for binidx in range(1,histnbins):
	        sspan_bin_spec.append(smult * (math.pow(spower,binidx)))

	sspan_bin_spec.append(smax)

if len(usspanbinedges) > 1:
        usspan_bin_spec = usspanbinedges
else:
   	sys.error("undefined bin edges for usspanbinegdes")

if len(avsspanbinedges) > 1:
        avsspan_bin_spec = avsspanbinedges
else:
        sys.error("undefined bin edges for usspanbinegdes")



pstart = 0
ppower = 1 + 13.0/float(histnbins)
pmult = 500.0/float(histnbins)
prange = int(100000/histnbins)
pmax = 2000000

pspan_bin_spec = [pstart]
for binidx in range(1,histnbins):
        pspan_bin_spec.append(pmult * (math.pow(ppower,binidx)))

pspan_bin_spec.append(pmax)


def plotnow(datain,mintsin,maxtsin,nyearsin,span_bins_spec,histnbins,title,Datatitle,figtitlein):
	
	plotin = np.zeros([nyearsin,histnbins])

	for iyear in range(0,nyearsin):
		plotin[iyear,:],span_bin = np.histogram(datain[mintsin[iyear]:maxtsin[iyear]],bins=span_bins_spec,density=plotdensity)

	span_bins_plot = np.zeros(histnbins)

	for ibin in range(0,histnbins):
		span_bins_plot[ibin] = 0.5 *(span_bin[ibin] + span_bin[ibin+1])

	if title == "SpatSpan" or title == "UniSpatSpan" or title == "AvSpatSpan":
		unit = "gridboxes"
		convert = 1.0
	elif title == "TimeSpan":
		unit = "days"
		convert = 8.0
	else:
		unit = "unspecified"

	yearnums = range(0,nyearsin)
	A = np.array(yearnums)
	regress = (np.zeros((5,histnbins),np.float))

	for ibin in range(0,histnbins):
		linreg = plotin[:,ibin]
		regress[:,ibin] = stats.linregress(A,linreg)

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

	res.tiYAxisString = "year"

	plot = []

		
	if (plotdensity):
		res.tiYAxisString = "frequency"
	else:
		res.tiYAxisString = "number"


	for ibin in range(0,histnbins):
		res.tiMainString = title  + ' ' + '{:2.1f}'.format(span_bins_spec[ibin]/convert) + '-' + '{:2.1f}'.format(span_bins_spec[ibin+1]/convert) + unit + '; r2= ' + '{:5.3f}'.format(regress[2,ibin] * regress[2,ibin]) + '; p= ' + '{:5.3f}'.format(regress[3,ibin])

		plot.append(Ngl.xy(wks,range(startyr,endyr+1),plotin[:,ibin],res))

	panelres = Ngl.Resources()
	panelres.nglPanelLabelBar = True
	panelres.nglPanelYWhiteSpacePercent = 5.
	panelres.nglPanelXWhiteSpacePercent = 5.

	panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
	panelres.nglPanelTop                      = 0.98
	panelres.nglPanelBottom                      = 0.02
	panelres.nglPaperOrientation = "Auto"

	txres = Ngl.Resources()
	txres.txFontHeightF = 0.012
	Ngl.text_ndc(wks,'Annual timeseries of global number of events of various scales from' + Datatitle,0.5,0.94,txres)
	print 'about to panel'


	Ngl.panel(wks,plot,[math.ceil(float(histnbins)/3.0),3],panelres)

	print 'panelled'


# Now plot


# Time span
variable = "TimeSpan"
vartitle = "0-1day"
if plotdensity:
        figtitle = Data + "_" + Version + '_' + variable + vartitle + '_' + str(histnbins) + 'bins_density.nc'
else:
        figtitle = Data + "_" + Version + '_' + variable + vartitle +'_' + str(histnbins) + 'bins_number_' + str(startyr) + '-' + str(endyr) + '.nc'

plotnow(tspan,mints,maxts,nyears,tspan_bin_spec,histnbins,variable,Data,FigDir + figtitle)


# Time span
tspan_bin_spec2 = tspanbinedges2

variable = "TimeSpan"
vartitle = "0-7day"
if plotdensity:
        figtitle = Data + "_" + Version + '_' + variable + vartitle + '_' + str(histnbins) + 'bins_density.nc'
else:
        figtitle = Data + "_" + Version + '_' + variable + vartitle +'_' + str(histnbins) + 'bins_number_' + str(startyr) + '-' + str(endyr) + '.nc'

plotnow(tspan,mints,maxts,nyears,tspan_bin_spec2,histnbins,variable,Data,FigDir + figtitle)



# Spatial span
variable = "AvSpatSpan"
if plotdensity:
        figtitle = Data + "_" + Version + '_' + variable + '_' + str(histnbins) + 'bins_density.nc'
else:
        figtitle = Data + "_" + Version + '_' + variable + '_' + str(histnbins) + 'bins_number_' + str(startyr) + '-' + str(endyr) + '.nc'

plotnow(avsspan,mints,maxts,nyears,avsspan_bin_spec,histnbins,variable,Data,FigDir + figtitle)

# Unique Spatial span
if plotdensity:
        figtitle = Data + "_" + Version + '_UniSpatSpan' + str(histnbins) + 'bins_density.nc'
else:
        figtitle = Data + "_" + Version + '_UniSpatSpan_' + str(histnbins) + 'bins_number_' + str(startyr) + '-' + str(endyr) + '.nc'

plotnow(usspan,mints,maxts,nyears,usspan_bin_spec,histnbins,"UniSpatSpan",Data,FigDir + figtitle)


datain.close()

Ngl.end()









