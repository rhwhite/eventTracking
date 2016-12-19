# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 17:45:25 2015
Code to plot figures of precipitation events density and total rainfall using command line arguments
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
import math
from scipy import stats
from rhwhitepackages.readwrite import shiftlons
from rhwhitepackages.readwrite import xrayOpen
import argparse
import resource

rsrcV = resource.RLIMIT_AS
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit starts as  :', soft
print 'Hard limit starts as :', hard
resource.setrlimit(rsrcV, (50000000000, hard)) #limit memory usage
#                          137438953472
soft, hard = resource.getrlimit(rsrcV)
print 'Soft limit changed to :', soft

parser = argparse.ArgumentParser(description="map event data")
parser.add_argument('--splittype',metavar='splittype',type=str,nargs=1,help='the type of split you want, day, speed, or maxspeed')
parser.add_argument('--speedtspan',metavar='speedtspan',type=int,nargs='?',default=4,help='how many time spans does the speed average cover?')
parser.add_argument('--tbound1',metavar='tbound1',type=float,nargs='+',help='lower bounds')
parser.add_argument('--tbound2',metavar='tbound2',type=float,nargs="+",help='upper bounds')
parser.add_argument('--unit',type=str,nargs=1,help='units of split type')
parser.add_argument('--Data',type=str,nargs=1,help='type of Data, TRMM, ERAI, or CESM')
parser.add_argument('--Version',type=str,nargs=1,help='Version of Data, Standard, low, 6th_from6 etc')
parser.add_argument('--anstartyr',type=int,nargs=1,help='start year for analysis')
parser.add_argument('--anendyr',type=int,nargs=1,help='end year for analysis')
parser.add_argument('--minlat',type=int,nargs='?',default=-45,help='min lat')
parser.add_argument('--maxlat',type=int,nargs='?',default=45,help='max lat')
parser.add_argument('--minlon',type=int,nargs='?',default=0,help='min lon')
parser.add_argument('--maxlon',type=int,nargs='?',default=360,help='max lon')
parser.add_argument('--sumlats',type=int,nargs='?',default=-1,help='number of lats to sum over')
parser.add_argument('--sumlons',type=int,nargs='?',default=-1,help='number of lons to sum over')
parser.add_argument('--minGB',type=int,nargs='?',default=0,help='minimum number of gridboxes to count as an event')


args = parser.parse_args()

print args.Data

splittype = args.splittype[0]
speedtspan = args.speedtspan
tbound1 = args.tbound1
tbound2 = args.tbound2
unit = args.unit[0]
Data = args.Data[0]
Version = args.Version[0]
anstartyr = args.anstartyr[0]
anendyr = args.anendyr[0]
MinLonF = args.minlon
MaxLonF = args.maxlon
MinLatF = args.minlat
MaxLatF = args.maxlat
sumlats = args.sumlats
sumlons = args.sumlons
minGB = args.minGB

if splittype=='day':
    splitname = 'Sizes'
elif splittype =='maxspeed':
    splitname = 'MaxSpeeds_' + str(speedtspan) + 'ts'

nbounds = len(tbound1)

FillValue = -9999

# Time period for analysis
nyears = anendyr - anstartyr + 1
print nyears
mapping = "center"
FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

if Data == "TRMM":
    Fstartyr = 1998
    Fendyr = 2014
    PrecipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"

    if sumlats > 0:
        PrecipClimFile = "Regrid_" + str(sumlats) + "_" + str(sumlons) + "_SeasAnn_TRMM_1998-2014_3B42_3hrly_nonan.nc"
    else:
        PrecipClimFile = "TRMM_1998-2014_clim_ann_1998-2014.nc"


    if Version == '5th_nanto25':
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_nto25/Precip/'
    elif Version == '5th_nantozero':
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5thresh_n2z/Precip/'
    elif Version == '7thresh':
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/7thresh/Precip/'
    elif Version == '6th_from6':
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/6th_from6/Precip/'
    elif Version == '5th_from48':
        DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW/TRMM_output/5th_from48/Precip/'

elif Data == "TRMMERAIgd":
    Fstartyr = 1998
    Fendyr = 2014
    PrecipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"
    PrecipClimFile = "Regrid_" + str(sumlats) + "_" + str(sumlons) + "_regrid2ERAI_TRMM_3B42_1998-2014_annclim.nc"

elif Data == "ERAI":
    Fstartyr = 1980
    Fendyr = 2014
    PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERAI/Precip_3hrly/'
    if sumlats > 0:
        PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ncra_ERAI_Totalprecip_1980-2015_preprocess.nc' #'SeasAnn_ERAI_Totalprecip_1980-2015_preprocess.nc'
    else:
                PrecipClimFile = 'Ann_ERAI_Totalprecip_1980-2015_preprocess.nc' #'SeasAnn_ERAI_Totalprecip_1980-2015_preprocess.nc'


elif Data == "ERA20C":
    Fstartyr = 1980
    Fendyr = 2011
    PrecipClimDir = '/home/disk/eos4/rachel/Obs/ERA_20C/'
    if sumlats > 0:
        PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ERA_20C_Ann_Totalprecip_' + str(Fstartyr) + '-' + str(Fendyr) + '.nc'
    else:
        PrecipClimFile = 'ERA_20C_Ann_Totalprecip_' + str(Fstartyr) + '-' + str(Fendyr) + '.nc'


elif Data == "CESM":
    Fstartyr = 1990
    Fendyr = 2010
    DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/CESM_output/' + Version + str(Fstartyr) + '/proc/'
    
    PrecipClimDir = '/home/disk/eos4/rachel/EventTracking/Inputs/CESM/f.e13.FAMPIC5.ne120_ne120.1979_2012.001/'
    if sumlats > 0:
        PrecipClimFile = 'Regrid_' + str(sumlats) + '_' + str(sumlons) + '_ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'
    else:
        PrecipClimFile = 'ncra_f.e13.FAMIPC5.ne120_ne120_TotalPrecip_1979-2012.nc'

DirIn = '/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + Data + '_output/' + Version + str(Fstartyr) + '/proc/'

FileInPrecip = xrayOpen(PrecipClimDir + PrecipClimFile)

if Data in ['TRMM']:
    if sumlats > 0:
        latin = FileInPrecip['lat']
        invar = 'PrecipAnnClim'
    else:
        latin = FileInPrecip['latitude']
        invar = 'pcp'

elif Data in ['TRMMERAIgd']:
    latin = FileInPrecip['lat']
    invar = 'pcp'

elif Data in ["ERAI"]:
    latin = FileInPrecip['lat']
    invar = 'tpnew'

elif Data in ["ERA20C"]:
    latin = FileInPrecip['lat']
    invar = 'tpnew'

elif Data in ["CESM"]:
    latin = FileInPrecip['lat']
    invar = 'PRECT'

if latin[0] > latin[1]:
    PrecipIn = FileInPrecip[invar][:,::-1,:].sel(lat=slice(MinLatF,MaxLatF))
else:
    if sumlats > 0:
        PrecipIn = FileInPrecip[invar].sel(lat=slice(MinLatF,MaxLatF))
    else:
        PrecipIn = FileInPrecip[invar].sel(latitude=slice(MinLatF,MaxLatF))

try:
    print np.amax(PrecipIn)
    if PrecipIn.units == "mm/hr":
        print 'converting'
        PrecipIn = PrecipIn * 24.0 #convert to mm/day
    elif PrecipIn.units == "mm/day":
        pass    
    else:
        error("unexpected unit in Precip file")
except AttributeError:
    if np.amax(PrecipIn) < 1.0:
        print "guessing we need to convert precip units! You should check this!"
        PrecipIn = PrecipIn * 24.0 #convert to mm/day
    else:
        print "guessing that we don't need to convert precip units - you should check this!"

def regressmaps(m,c,r,p,stderr,A,linreg,nlats,nlons):
# Calculate linear trend
    for ilat in range(0,nlats):
        for ilon in range(0,nlons):
            m[ilat,ilon],c[ilat,ilon],r[ilat,ilon],p[ilat,ilon],stderr[ilat,ilon] = stats.linregress(A,linreg[:,ilat,ilon])

def plotmap(plotvars1,plotvars2,plotmin1,plotmax1,plotmin2,plotmax2,vartitle1,vartitle2,title,figtitle,lons,lats,minlon,maxlon,minlat,maxlat):
    print plotmin1
    nplots = plotvars1.shape[0]
    print nplots
    wkres = Ngl.Resources()
    wkres.wkColorMap = "precip_diff_12lev"
    wks_type = "eps"
    wks = Ngl.open_wks(wks_type,figtitle,wkres)

    res = Ngl.Resources()
    res.cnInfoLabelOn         = False    # Turn off informational
                              # label.
    res.pmLabelBarDisplayMode = "Always" # Turn on label bar.
    res.cnLinesOn             = False    # Turn off contour lines.
    res.nglDraw  = False
    res.nglFrame = False

    res.sfMissingValueV = FillValue

    res.cnFillOn = True
    res.cnMissingValFillColor = "white"
    res.cnLineLabelsOn       = False
    res.pmLabelBarDisplayMode = "Always"
    res.cnLinesOn =  False

    # if lons start negative, shift everything over so there isn't a line down the middle of the Pacific
    if lons[0] < 0:
        nlonhalf = nlons/2
        lonsnew = np.zeros(lons.shape,np.float)
        lonsnew[0:nlonhalf] = lons[nlonhalf:nlons]
        lonsnew[nlonhalf:nlons] = lons[0:nlonhalf] + 360.0
        lons = lonsnew

        for iplot in range(0,nplots):
            plotvars1[iplot] = shiftlons(plotvars1[iplot],lons)
            plotvars2[iplot] = shiftlons(plotvars2[iplot],lons)
    else:
        lonsnew = lons

    res.sfXCStartV = float(lonsnew[0])
    res.sfXCEndV = float(lonsnew[len(lons)-1])
    res.sfYCStartV = float(lats[0])
    res.sfYCEndV = float(lats[len(lats)-1])

    res.mpProjection = "CylindricalEquidistant" # Change the map projection.
    res.mpCenterLonF = 180.           # Rotate the projection.
    res.mpFillOn     = True           # Turn on map fill.

    res.lbOrientation   = "Vertical"
    res.mpLimitMode = "LatLon"    # Limit the map view.
    res.mpMinLonF = MinLonF
    res.mpMaxLonF = MaxLonF
    res.mpMinLatF = MinLatF
    res.mpMaxLatF = MaxLatF
    res.mpOutlineBoundarySets = "AllBoundaries"

    res.lbLabelFontHeightF = 0.0125
    res.lbTitleFontHeightF = 0.0125
    
    res.tiMainFontHeightF = 0.015

    res.cnLevelSelectionMode = "ManualLevels" # Define your own

    toplot = []


    for iplot in range(0,nplots):
        tempplot = plotvars1[iplot]
        tempplot[np.where(np.isnan(tempplot))] = FillValue
        res.cnMinLevelValF       = plotmin1[iplot]          # contour levels.
        res.cnMaxLevelValF       = plotmax1[iplot]
        res.cnLevelSpacingF      = ((plotmax1[iplot]-plotmin1[iplot])/10.0)
        if iplot == 0:
            res.tiMainString = vartitle1[iplot] +"; mean = {:2.3g}".format(np.nanmean(tempplot))
        else:
            res.tiMainString = vartitle1[iplot] + "; mean = {:2.3g}".format(np.nanmean(tempplot)) + "%"
        toplot.append(Ngl.contour_map(wks,tempplot,res))

        tempplot = plotvars2[iplot]
        tempplot[np.where(np.isnan(tempplot))] = FillValue
        res.cnMinLevelValF       = plotmin2[iplot]          # contour levels.
        res.cnMaxLevelValF       = plotmax2[iplot]
        res.cnLevelSpacingF      = ((plotmax2[iplot]-plotmin2[iplot])/10.0)
        
        if iplot == 0:
            res.tiMainString = vartitle2[iplot] + "; mean = {:2.3g}".format(SumPrecipPercent) + "%"
        else:
            res.tiMainString = vartitle2[iplot] + "; mean = {:2.3g}".format(np.nanmean(tempplot)) + "%"
            toplot.append(Ngl.contour_map(wks,tempplot,res))

    
    textres = Ngl.Resources()
    textres.txFontHeightF = 0.015
    Ngl.text_ndc(wks,title,0.5,0.87,textres)


    panelres = Ngl.Resources()
    panelres.nglPanelLabelBar = True
    #panelres.nglPanelYWhiteSpacePercent = 5.
    #panelres.nglPanelXWhiteSpacePercent = 5.

    panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
    if nplots > 5:
        panelres.nglPanelTop                      = 0.8
        panelres.nglPanelBottom                      = 0.15
    else:
        panelres.nglPanelTop                      = 0.95
        panelres.nglPanelBottom                      = 0.01


    #panelres.nglPanelFigureStrings            = ["a","b","c","d","e","f"]
    #panelres.nglPanelFigureStringsJust        = "BottomLeft"

    panelres.nglPaperOrientation = "Auto"

    plot = Ngl.panel(wks,toplot,[nplots,2],panelres)

# Get lats and lons
iday = 0
if minGB > 0:
    fileadd = '_min' + str(minGB) + 'GB'
else:
    fileadd = ''

if sumlats > 0:
    FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
else:
    FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd +  '.nc'

print DirIn + FileIn
FileIn = xrayOpen(DirIn + FileIn)

print MinLatF, MaxLatF
lats = FileIn['lat'].sel(lat=slice(MinLatF,MaxLatF)).values
lons = FileIn['lon'].values

years = FileIn['years'].sel(years=slice(anstartyr,anendyr))

print lats

nlats = len(lats)
nlons = len(lons)
nyears = len(years)

difflat = lats[1]-lats[0]
difflon = lons[1]-lons[0]

if splittype == 'day':
    arraynbounds = nbounds + 1
else:
    arraynbounds = nbounds

Dendays = np.zeros([arraynbounds,nyears,nlats,nlons])      # +1 for total sum as first plot
Precipdays = np.zeros([arraynbounds,nyears,nlats,nlons])   # + 1 for total precip as first plot

for iday in range(0,nbounds):
    if sumlats > 0:
        FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
    else:
        FileIn = 'DenDirSpd_Map_Ann_' + splitname + '_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + Data + "_" + str(Fstartyr) + '-' + str(Fendyr) + '_' + Version + fileadd + '.nc'

    #Get lons and lats
    print DirIn + FileIn
    FileIn = xray.open_dataset(DirIn + FileIn)

    Dendays[iday,:,:,:] = FileIn['TDensity'].sel(lat=slice(MinLatF,MaxLatF),years=slice(anstartyr,anendyr)) 
    print tbound1,tbound2

    Precipdays[iday,:,:,:] = FileIn['TPrecip'].sel(lat=slice(MinLatF,MaxLatF),years=slice(anstartyr,anendyr))


print Dendays.shape

DenAnnAvg = np.nanmean(Dendays,axis=1)
PrecipAnnAvg = np.nanmean(Precipdays,axis=1)

DenAll = np.nansum(DenAnnAvg,axis=0)
PrecipAll = np.nansum(PrecipAnnAvg,axis=0)

#PrecipAllPercent = PrecipAll/365.0 # Convert from mm/year to mm/day
PrecipAllPercent = 100.0 * np.divide(PrecipAll,PrecipIn*365.0)  # convert PrecipIn from mm/day to mm/yr before dividing

SumPrecipPercent = 100.0 * np.nansum(PrecipAll)/np.nansum(PrecipIn * 365.0)
DenPercent = np.zeros(DenAnnAvg.shape)
PrecipPercent = np.zeros(PrecipAnnAvg.shape)

titlesDen = []
titlesPrec = []

arrayindex = 0
if splittype == 'day':
    DenPercent[0,:,:] = DenAll / (difflat * difflon)
    PrecipPercent[0,:,:] = PrecipAllPercent # Percent of total TRMM precip falling in these events
    titlesDen.append("Total annual event density, events/(yr deg~S1~2 )")
    titlesPrec.append("Percentage of total precipitation captured in events")
    arrayindex += 1

for iday in range(0,nbounds):
    DenPercent[iday+arrayindex,:,:] = 100.0 * (np.where(DenAll > 0, np.divide(DenAnnAvg[iday,:,:],DenAll),0.0)) 
    PrecipPercent[iday+arrayindex,:,:] = 100.0 * (np.where(PrecipAll > 0, np.divide(PrecipAnnAvg[iday,:,:],PrecipAll),0.0))

    if tbound2[iday] < tbound1[iday]*10 or tbound1[iday] == 0 or tbound2[iday] == 0:
        titlesDen.append('Event density: ' + str(tbound1[iday]) + ' to ' + str(tbound2[iday]) + unit + ' events')
        titlesPrec.append('Event precip: ' + str(tbound1[iday]) + ' to ' + str(tbound2[iday]) + unit + ' events')
    else:
        titlesDen.append("Event density: >" + str(tbound1[iday]) + unit + ' events')
        titlesPrec.append("Event precip: >" + str(tbound1[iday]) + unit + ' events')
    print np.nanmean(DenPercent[iday+arrayindex,:,:])
    print np.nanmean(PrecipPercent[iday+arrayindex,:,:])


# And now plot 

# Plot 1: Average density, and percentage easterly

figtitlein = FigDir + 'Paper_DenPrecipClim_' + Data + '_' + Version + '_Ana' + str(anstartyr) + '-' + str(anendyr) + '_' + str(tbound1[0]) + '_to_' + str(tbound2[nbounds-1]) + unit + fileadd
titlein = 'Events in ' + Data + " " + Version + ' years ' + str(anstartyr) + '-' + str(anendyr)

if splittype == 'day':
    if Data in ["TRMM"]:
        clims1,clims2,clims3,clims4 = [0.0,98.0,0.0,0.0,0.0,0.0],[300,100.0,1.0,0.2,0.05,1.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,90.0,40.0,40.0,40.0,40.0]  
    elif Data in ["TRMMERAIgd"]:
        if minGB in [4,9]:
            clims1,clims2,clims3,clims4 = [0.0,70.0,0.0,0.0,0.0,0.0],[6,100.0,20.0,7.0,4.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,50.0,50.0,50.0,50.0]

        else:
            clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[10,100.0,15.0,5.0,2.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,40.0,50.0,50.0,50.0]
    elif Data in ["ERAI"]:
        if minGB in [4,9]:
            clims1,clims2,clims3,clims4 = [0.0,70.0,0.0,0.0,0.0,0.0],[6,100.0,20.0,7.0,4.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,50.0,50.0,50.0,50.0]
        else:
            clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[10,100.0,15.0,5.0,2.0,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,80.0,40.0,50.0,50.0,50.0]

#       clims1,clims2,clims3,clims4 = [0.0,95.0,0.0,0.0,0.0,0.0],[30,100.0,2.0,1.5,0.5,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]
    elif Data in ["ERA20C"]:
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[15,100.0,4.0,3.0,2.0,4.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]        
    elif Data in ["CESM"]:
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0],[50,100.0,5.0,3.0,0.5],[0.0,0.0,0.0,0.0,0.0],[100,90.0,40.0,40.0,40.0]  

elif splittype == 'maxspeed':
    if Data == "TRMM":
        clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,95.0,0.0,0.0,0.0],[0.5,0.75,5.0,100.0,5.0,0.75,0.5],[0.0,0.0,0.0,20.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]
    elif Data in ["TRMMERAIgd"]:
        clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,90.0,0.0,0.0,0.0],[0.5,1.0,5.0,100.0,5.0,1.0,0.5],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]
    elif Data in ["ERAI"]:
        clims1,clims2,clims3,clims4 = [0.0,0.0,0.0,90.0,0.0,0.0,0.0],[0.5,1.0,5.0,100.0,5.0,1.0,0.5],[0.0,0.0,0.0,0.0,0.0,0.0,0.0],[30,30.0,30.0,100.0,30.0,30.0,30.0]
#               clims1,clims2,clims3,clims4 = [0.0,95.0,0.0,0.0,0.0,0.0],[30,100.0,2.0,1.5,0.5,3.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]

    elif Data == "CESM":
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0],[50,100.0,5.0,3.0,0.5],[0.0,0.0,0.0,0.0,0.0],[20,90.0,40.0,40.0,40.0]

    elif Data in ["ERA20C"]:
        clims1,clims2,clims3,clims4 = [0.0,90.0,0.0,0.0,0.0,0.0],[15,100.0,4.0,3.0,2.0,4.0],[0.0,0.0,0.0,0.0,0.0,0.0],[100,70.0,30.0,50.0,50.0,50.0]

plotmap(DenPercent,PrecipPercent,clims1,clims2,clims3,clims4 ,titlesDen,titlesPrec,titlein,figtitlein,lons,lats,MinLonF,MaxLonF,MinLatF,MaxLatF)

# 



Ngl.end()




