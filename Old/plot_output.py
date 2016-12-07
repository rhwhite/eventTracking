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

histnbins = 40 # number of bins

startyr = 1998
endyr = 2014 

plotdensity = False

mintsteps = 0	# 
maxtsteps = 23360 # number of timesteps (3hrly) to include; 2920 is one year # 43800 is 15 years, 46720 is 16 years

minevent = 100000

DirI = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Precip/'
FileI1 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + 'tspan.nc'
FileI2 = 'Precip_Sizes_' + str(startyr) + '-' + str(endyr) + '.nc'

DirO = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/Precip/'
FileO = 'hist_precip_tspan_sspan_bins_' + str(mintsteps) + '-' + str(maxtsteps) + '_' + str(plotdensity) + '.nc'

FigDir = FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'
figtitle1 = 'precip_tspan_sspan_' + str(histnbins) + 'bins_' + str(mintsteps) + '-' + str(maxtsteps) + '_' + str(plotdensity) + '.nc' 
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

curtime = tmean[mintsteps]
mintimeidx = 0

while curtime < mintsteps:
	mintimeidx = mintimeidx + 1
	curtime = tstart[mintimeidx]


print mintimeidx
maxtimeidx = mintimeidx
curtime = tmean[maxtimeidx]
print curtime

while curtime < maxtsteps:
	maxtimeidx = maxtimeidx + 1
	curtime = tstart[maxtimeidx]

print nevents
print 'maximum event is: ', maxtimeidx

#Bin data 

try:
    os.remove(DirO + FileO)
except OSError:
    pass
tstart = 2
trange = 2

sstart = 1
srange = 750

pstart = 0
prange = 10000

tspan_bin_spec = range(tstart,tstart + trange * (histnbins+1),trange) 
sspan_bin_spec = range(sstart,sstart + srange * (histnbins+1),srange)
pspan_bin_spec = range(pstart,pstart + prange * (histnbins+1),prange)


tspanplot,tspan_bin = np.histogram(tspan[mintimeidx:maxtimeidx],bins=tspan_bin_spec,density=plotdensity)
sspanplot,sspan_bin = np.histogram(sspan[mintimeidx:maxtimeidx],bins=sspan_bin_spec,density=plotdensity)
precipplot,precip_bin = np.histogram(totalprecip[mintimeidx:maxtimeidx],bins=pspan_bin_spec,density=plotdensity)

print sspan_bin
print sspanplot

tspan_bins_plot = np.zeros(histnbins)
sspan_bins_plot = np.zeros(histnbins)
precip_bins_plot = np.zeros(histnbins)


for ibin in range(0,histnbins):
        tspan_bins_plot[ibin] = 0.5 *(tspan_bin[ibin] + tspan_bin[ibin+1])
        sspan_bins_plot[ibin] = 0.5 *(sspan_bin[ibin] + sspan_bin[ibin+1])
        precip_bins_plot[ibin] = 0.5 *(precip_bin[ibin] + precip_bin[ibin+1])

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
res.xyYStyle = "Log"
res.xyXStyle = "Log"

res.tiMainString = 'LogLogplots_' + str(histnbins) + 'bins_from_' + str(mintsteps) + '-' + str(maxtsteps) + '_' + str(plotdensity) + '.nc' 

plot = []

if (plotdensity):
	res.tiYAxisString = "frequency"
else:
	res.tiYAxisString = "number"

res.tiXAxisString = "timespan (in 3 hourly increments)"
plot.append(Ngl.xy(wks,tspan_bins_plot,tspanplot,res))

res.tiXAxisString = "spatial span (in TRMM gridboxes)"
plot.append(Ngl.xy(wks,sspan_bins_plot,sspanplot,res))

res.tiXAxisString = "precip total (in mm)"
plot.append(Ngl.xy(wks,precip_bins_plot,precipplot,res))



panelres = Ngl.Resources()
panelres.nglPanelLabelBar = True
panelres.nglPanelYWhiteSpacePercent = 5.
panelres.nglPanelXWhiteSpacePercent = 5.

panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPanelTop                      = 0.98
panelres.nglPanelBottom                      = 0.02
panelres.nglPaperOrientation = "Auto"

Ngl.panel(wks,plot,[3,1],panelres)

Ngl.end()



