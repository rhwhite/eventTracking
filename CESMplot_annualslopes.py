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

Data = "CESM"

test = 0
histnbins = 4 # number of bins
tspanbinedges = [0,7.9,15.9,39.9,500]#6.5,7.5,8.5]
avsspanbinedges = [0,30000,50000,200000,5000000]#150,350,10000]

precipbinedges = [0,2.0E9,6.0E9,6.0E10,1.0E16]
regresst = (np.zeros((5,histnbins),np.float))


if Data == "CESM":
	Fstartyr = 1990
	Fendyr = 2014

	startyr = 1990
	endyr = 2014 
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

if test == 1:
	endyr = startyr
nyears = endyr - startyr + 1

mints = np.zeros(nyears)
maxts = np.zeros(nyears)

plotdensity = False

starttsteps = 0
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

totprecip = np.divide(totalprecip,1000.0)	# Convert to m3 from mm m2
avsspan = np.divide(sspan,tspan*1000.0*1000.0)	# Find average and convert to km2

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

print (maxts - mints)

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
"""
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
   	sys.error("undefined bin edges for usspanbinedges")
"""
if len(avsspanbinedges) > 1:
        avsspan_bin_spec = avsspanbinedges
else:
        sys.error("undefined bin edges for usspanbinedges")

if len(precipbinedges) > 1:
        precip_bin_spec = precipbinedges
else:
        sys.error("undefined bin edges for precipbinedges")


def plotnow(nlines,datain,mintsin,maxtsin,nyearsin,span_bins_spec,histnbins,title,Datatitle,figtitlein):
	
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

	res.trXMaxF = 2014
	res.trYMaxF = 1997 

	plot = []

	if (plotdensity):
		res.tiYAxisString = "frequency"
	else:
		res.tiYAxisString = "number of events"

	for iline in range(0,nlines):
	        plotin = np.zeros([nyearsin,histnbins])

		for iyear in range(0,nyearsin):
			plotin[iyear,:],span_bin = np.histogram(datain[iline][mintsin[iyear]:maxtsin[iyear]],bins=span_bins_spec[iline][:],density=plotdensity)

		span_bins_plot = np.zeros(histnbins)

		for ibin in range(0,histnbins):
			span_bins_plot[ibin] = 0.5 *(span_bin[ibin] + span_bin[ibin+1])

		if title[iline] in ["Footprint"]:
			unit = "km~S1~2"
			convert = 1.0
		elif title[iline] in ["Total Precip"]:
			unit = "m~S1~3"
			convert = 1.0
		elif title[iline] in ["UniSpatSpan", "AvSpatSpan"]:
			unit = "gridboxes"
			convert = 1.0
		elif title[iline] in ["Timespan","TimeSpan"]:
			unit = "days"
			convert = 8.0
		else:
			unit = "unspecified"

		yearnums = range(0,nyearsin)
		A = np.array(yearnums)
		regress = (np.zeros((5,histnbins),np.float))

		for ibin in range(0,histnbins):
			if ibin == 0:
				if (plotdensity):
					res.tiYAxisString = "frequency"
				else:
					res.tiYAxisString = "number of events"
			else:
				res.tiYAxisString = ""

			linreg = plotin[:,ibin]
			regress[:,ibin] = stats.linregress(A,linreg)

			if ibin == histnbins -1:
				res.tiMainString = '{:^80}'.format('          >' + '{:2.1g}'.format(span_bins_spec[iline][ibin]/convert) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]) + "           ")
			else:
				res.tiMainString = '{:^80}'.format('{:2.1g}'.format(span_bins_spec[iline][ibin]/convert) + '-' + '{:2.1g}'.format(span_bins_spec[iline][ibin+1]/convert) + unit + '; p=' + '{:5.3f}'.format(regress[3,ibin]))

			plot.append(Ngl.xy(wks,range(startyr,endyr+1),plotin[:,ibin],res))

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

	Ngl.panel(wks,plot,[nlines,float(histnbins)],panelres)

	print 'panelled'


# Now plot


# Time span
variable = ["Timespan","Footprint","Total Precip"]
vartitle = "annualtrends"
if plotdensity:
        figtitle = "CESM2Paper_" + Data + "_" + Version + '_' + vartitle + '_dens'
else:
        figtitle = "CESM2Paper_" + Data + "_" + Version + '_' + vartitle + '_numb'

plotarray = [tspan,avsspan,totprecip]
plotarraybins = [tspan_bin_spec,avsspan_bin_spec,precip_bin_spec]
ntypes = len(plotarray)
print ntypes
plotnow(ntypes,plotarray,mints,maxts,nyears,plotarraybins,histnbins,variable,Data,FigDir + figtitle)


datain.close()

Ngl.end()









