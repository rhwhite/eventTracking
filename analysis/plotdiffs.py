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


histnbins = 10 # number of bins

startyr = 1998
endyr = 2014 

plotdensity = True

mintsteps1 = 5842	# 11688 is 4 years 
maxtsteps1 = 8762 # number of timesteps (3hrly) to include; 2920 is one year, 23376 8 years,  46720 is 16 years

mintsteps2 = 8763  # 35064 is 4 years from end 
maxtsteps2 = 11683 # number of timesteps (3hrly) to include; 2920 is one year # 43800 is 15 years, 46752 is 16 years with leap year


minevent = 100000

DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Precip/'
FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + 'tspan.nc'
FileI2 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '.nc'

FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
if (plotdensity):
	figtitle1 = 'logs_precip_tspan_sspan_' + str(histnbins) + 'bins_diffs_' + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)  + '_density2.nc'
else:
        figtitle1 = 'logs_precip_tspan_sspan_' + str(histnbins) + 'bins_diffs_' + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)  + '2.nc'
 
FigDir + figtitle1

filetimespan = "3hrly"


datain = xray.open_dataset(DirI + FileI1)
datain2 = xray.open_dataset(DirI + FileI2)

print(datain.coords)

tspan=datain['timespan']
tmean = datain['tmean']
tstart = datain2['tstart'].values
sspan = datain2['span']
totalprecip = datain2['totalprecip']
xmin1 = datain['xmin']
xmin2 = datain2['xmin']

if (xmin1 != xmin2).any():
	sys.exit('xmin from two files not equal')

nevents = tspan.shape[0]

curtime = tmean[mintsteps1]
mintimeidx1 = 0

while curtime < mintsteps1:
	mintimeidx1 = mintimeidx1 + 1
	curtime = tstart[mintimeidx1]

print mintimeidx1
maxtimeidx1 = mintimeidx1
curtime = tmean[maxtimeidx1]

while curtime < maxtsteps1:
	maxtimeidx1 = maxtimeidx1 + 1
	curtime = tstart[maxtimeidx1]

print 'maximum event of 1st is: ', maxtimeidx1

mintimeidx2 = 0
curtime = tstart[mintimeidx2]

while curtime < mintsteps2:
        mintimeidx2 = mintimeidx2 + 1
        curtime = tstart[mintimeidx2]

print mintimeidx2
maxtimeidx2 = mintimeidx2
curtime = tmean[maxtimeidx2]

while curtime < maxtsteps2:
        maxtimeidx2 = maxtimeidx2 + 1
        curtime = tstart[maxtimeidx2]

print 'maximum event of 2nd is: ', maxtimeidx2



#Bin data 

tstart = 0.9
tmult = 30/float(histnbins)
tpower = 1 + 2.5/float(histnbins)
trange = int(100/histnbins)
print tpower
tmax = 500

sstart = 0.9
spower = 1 + 10.0/float(histnbins)
smult = 75.0/float(histnbins)
#print srange
smax = 200000

pstart = 0
ppower = 1 + 13.0/float(histnbins)
pmult = 500.0/float(histnbins)

prange = int(100000/histnbins)
print ppower
pmax = 2000000

tspan_bin_spec = [tstart]
sspan_bin_spec = [sstart]
pspan_bin_spec = [pstart]
for binidx in range(1,histnbins):
	tspan_bin_spec.append(tspan_bin_spec[binidx-1] + tmult * ((math.pow(tpower,binidx))))
        sspan_bin_spec.append(smult * (math.pow(spower,binidx)))
        pspan_bin_spec.append(pmult * (math.pow(ppower,binidx)))


#tspan_bin_spec = range(tstart,tstart + trange * (histnbins),trange) 
#sspan_bin_spec = range(sstart,sstart + srange * (histnbins),srange)
#pspan_bin_spec = range(pstart,pstart + prange * (histnbins),prange)

tspan_bin_spec.append(tmax)
sspan_bin_spec.append(smax)
pspan_bin_spec.append(pmax)

print tspan_bin_spec
print sspan_bin_spec
print pspan_bin_spec


tspanplot0,tspan_bin = np.histogram(tspan,bins=tspan_bin_spec,density=plotdensity)
sspanplot0,sspan_bin = np.histogram(sspan,bins=sspan_bin_spec,density=plotdensity)
precipplot0,precip_bin = np.histogram(totalprecip,bins=pspan_bin_spec,density=plotdensity)


tspanplot1,tspan_bin = np.histogram(tspan[mintimeidx1:maxtimeidx1],bins=tspan_bin_spec,density=plotdensity)
sspanplot1,sspan_bin = np.histogram(sspan[mintimeidx1:maxtimeidx1],bins=sspan_bin_spec,density=plotdensity)
precipplot1,precip_bin = np.histogram(totalprecip[mintimeidx1:maxtimeidx1],bins=pspan_bin_spec,density=plotdensity)


tspanplot2,tspan_bin = np.histogram(tspan[mintimeidx2:maxtimeidx2],bins=tspan_bin_spec,density=plotdensity)
sspanplot2,sspan_bin = np.histogram(sspan[mintimeidx2:maxtimeidx2],bins=sspan_bin_spec,density=plotdensity)
precipplot2,precip_bin = np.histogram(totalprecip[mintimeidx2:maxtimeidx2],bins=pspan_bin_spec,density=plotdensity)

tspanplot1N,tspan_bin = np.histogram(tspan[mintimeidx1:maxtimeidx1],bins=tspan_bin_spec,density='false')
sspanplot1N,sspan_bin = np.histogram(sspan[mintimeidx1:maxtimeidx1],bins=sspan_bin_spec,density='false')
precipplot1N,precip_bin = np.histogram(totalprecip[mintimeidx1:maxtimeidx1],bins=pspan_bin_spec,density='false')


tspanplot2N,tspan_bin = np.histogram(tspan[mintimeidx2:maxtimeidx2],bins=tspan_bin_spec,density='false')
sspanplot2N,sspan_bin = np.histogram(sspan[mintimeidx2:maxtimeidx2],bins=sspan_bin_spec,density='false')
precipplot2N,precip_bin = np.histogram(totalprecip[mintimeidx2:maxtimeidx2],bins=pspan_bin_spec,density='false')





tspanplot3 = tspanplot2 - tspanplot1
sspanplot3 = sspanplot2 - sspanplot1
precipplot3 = precipplot2 - precipplot1

tspanplot4 = 100.0 * (tspanplot2 - tspanplot1)/(tspanplot1)
sspanplot4 = 100.0 * (sspanplot2 - sspanplot1)/(sspanplot1)
precipplot4 = 100.0 * (precipplot2 - precipplot1)/(precipplot1)

tspanplot4N = 100.0 * (tspanplot2N - tspanplot1N)/(tspanplot1N)
sspanplot4N = 100.0 * (sspanplot2N - sspanplot1N)/(sspanplot1N)
precipplot4N = 100.0 * (precipplot2N - precipplot1N)/(precipplot1N)



log_tspanplot0 = np.log(tspanplot0)
log_tspanplot1 = np.log(tspanplot1)
log_tspanplot2 = np.log(tspanplot2)

log_sspanplot0 = np.log(sspanplot0)
log_sspanplot1 = np.log(sspanplot1)
log_sspanplot2 = np.log(sspanplot2)

log_precipplot0 = np.log(precipplot0)
log_precipplot1 = np.log(precipplot1)
log_precipplot2 = np.log(precipplot2)



#Get absolute values, take logs, then add negatives back in.
tspanplot3abs = abs(tspanplot3)
sspanplot3abs = abs(sspanplot3)
precipplot3abs = abs(precipplot3)

tspanfactors = tspanplot3/tspanplot3abs
sspanfactors = sspanplot3/sspanplot3abs
precipfactors = precipplot3/precipplot3abs

logtspanplot3abs = np.log(tspanplot3abs)
logsspanplot3abs = np.log(sspanplot3abs)
logprecipplot3abs = np.log(precipplot3abs)

log_tspanplot3 = logtspanplot3abs * tspanfactors
log_sspanplot3 = logsspanplot3abs * sspanfactors
log_precipplot3 = logprecipplot3abs * precipfactors


log_tspanplot3_2 = log_tspanplot2 - log_tspanplot1
log_sspanplot3_2 = log_sspanplot2 - log_sspanplot1
log_precipplot3_2 = log_precipplot2 - log_precipplot1


tspan_bins_plot = np.zeros(histnbins)
sspan_bins_plot = np.zeros(histnbins)
precip_bins_plot = np.zeros(histnbins)

tspanplot4[np.isinf(tspanplot4)] = np.nan
sspanplot4[np.isinf(sspanplot4)] = 0
precipplot4[np.isinf(precipplot4)] = 0
log_tspanplot1[np.isinf(tspanplot4)] = np.nan
log_tspanplot3[np.isinf(tspanplot4)] = np.nan
log_tspanplot3_2[np.isinf(tspanplot4)] = np.nan


for ibin in range(0,histnbins):
        tspan_bins_plot[ibin] = 0.5 *(tspan_bin[ibin] + tspan_bin[ibin+1])
        sspan_bins_plot[ibin] = 0.5 *(sspan_bin[ibin] + sspan_bin[ibin+1])
        precip_bins_plot[ibin] = 0.5 *(precip_bin[ibin] + precip_bin[ibin+1])


log_tspan_bins_plot = np.log(tspan_bins_plot)
log_sspan_bins_plot = np.log(sspan_bins_plot)
log_precip_bins_plot = np.log(precip_bins_plot)

print tspanplot1
print tspanplot2
print tspanplot3

print sspanplot1
print sspanplot2
print sspanplot3

print precipplot1
print precipplot2
print precipplot3

print log_tspan_bins_plot
print log_sspan_bins_plot
print log_precip_bins_plot

datain.close()
datain2.close()

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
res.xyMarkerColor = "black"
res.xyMarkers = 16
res.xyMarkerSizeF = 0.01
res.xyYStyle = "Linear"
res.xyXStyle = "Linear"

plot = []

if (plotdensity):
	res.tiYAxisString = "log frequency"
        res.tiMainString = 'Absolute differences '  + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)  +  '.nc'
else:
	res.tiMainString = 'Absolute differences '  + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)  +  '.nc'
	res.tiYAxisString = "log number"

print 'titledone'
res.tiMainString = 'Absolute values from 1998-2014 data'

res.tiXAxisString = "timespan (in 3 hourly increments) all data"
plot.append(Ngl.xy(wks,log_tspan_bins_plot,log_tspanplot0,res))

res.tiXAxisString = "spatial span (in TRMM gridboxes) all data"
plot.append(Ngl.xy(wks,log_sspan_bins_plot,log_sspanplot0,res))

res.tiXAxisString = "precip total (in mm) all data"
plot.append(Ngl.xy(wks,log_precip_bins_plot,log_precipplot0,res))

#res.tiXAxisString = "timespan (in 3 hourly increments) first half"
#plot.append(Ngl.xy(wks,log_tspan_bins_plot,log_tspanplot1,res))

#res.tiXAxisString = "spatial span (in TRMM gridboxes) first half"
#plot.append(Ngl.xy(wks,log_sspan_bins_plot,log_sspanplot1,res))

#res.tiXAxisString = "precip total (in mm) first half"
#plot.append(Ngl.xy(wks,log_precip_bins_plot,log_precipplot1,res))

#print 'first 3'
res.tiMainString = 'Absolute differences '  + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)


res.tiXAxisString = "timespan (in 3 hourly increments) diffs"
plot.append(Ngl.xy(wks,log_tspan_bins_plot,log_tspanplot3,res))

res.tiXAxisString = "spatial span (in TRMM gridboxes) diffs"
plot.append(Ngl.xy(wks,log_sspan_bins_plot,log_sspanplot3,res))

res.tiXAxisString = "precip total (in mm) second half"
plot.append(Ngl.xy(wks,log_precip_bins_plot,log_precipplot3,res))


#print 'next 3'

#res.tmYMajorGrid = True

#res.tiXAxisString = "log timespan (in 3 hourly increments)"
#plot.append(Ngl.xy(wks,log_tspan_bins_plot,log_tspanplot3,res))

#res.tiXAxisString = "log spatial span (in TRMM gridboxes)"
#plot.append(Ngl.xy(wks,log_sspan_bins_plot,log_sspanplot3,res))

#res.tiXAxisString = "log precip total (in mm)"
#plot.append(Ngl.xy(wks,log_precip_bins_plot,log_precipplot3,res))


print 'third 3'

res.tiYAxisString = "percentage change in density"

res.tiMainString = '% change (density)' + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)

res.tiXAxisString = "log of timespan (in 3 hourly increments)"
plot.append(Ngl.xy(wks,log_tspan_bins_plot,tspanplot4,res))

res.tiXAxisString = "log of spatial span (in TRMM gridboxes)"
plot.append(Ngl.xy(wks,log_sspan_bins_plot,sspanplot4,res))

res.tiXAxisString = "log of precip total (in mm)"
plot.append(Ngl.xy(wks,log_precip_bins_plot,precipplot4,res))


res.tiYAxisString = "percentage change in number"
res.tiMainString = '% change (number) ' + str(mintsteps2) + 'to' + str(maxtsteps2) + '-' + str(mintsteps1) + 'to' + str(maxtsteps1)

res.tiXAxisString = "log of timespan (in 3 hourly increments)"
plot.append(Ngl.xy(wks,log_tspan_bins_plot,tspanplot4,res))

res.tiXAxisString = "log of spatial span (in TRMM gridboxes)"
plot.append(Ngl.xy(wks,log_sspan_bins_plot,sspanplot4,res))

res.tiXAxisString = "log of precip total (in mm)"
plot.append(Ngl.xy(wks,log_precip_bins_plot,precipplot4,res))



panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.98
panelres.nglPanelBottom                      = 0.02
panelres.nglPaperOrientation = "Auto"

print 'about to panel'


Ngl.panel(wks,plot,[4,3],panelres)

print 'panelled'


Ngl.end()



