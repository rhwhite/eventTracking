
# coding: utf-8


# In[2]:

import os, errno
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray as xr
#import Ngl
#import math
from scipy import stats
from scipy.interpolate import griddata

from rhwhitepackages.readwrite import shiftlons
from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.stats import regressmaps

# plotting
# import matplotlib
import xray.plot as xplt
import matplotlib.pyplot as plt
import matplotlib.patches as patches

pthresh = 0.05

FigDir = '/home/disk/eos4/rachel/Figures/PrecipEvents/'

# In[3]:

#Select lats, lons, bounds, and other parameters
anstartyr = 1999
anendyr = 2014
splittype = 'day'
tbound1 = [0.0,1.0,2.0,5.0]
tbound2 = [1.0,2.0,5.0,100.0]

#tbound1 = [1.0, 1.125, 1.5, 1.875]
#tbound2 = [1.125, 1.25, 1.625, 2.0]

nbounds = len(tbound1)
if len(tbound2) != nbounds:
    exit("tbound arrays not equal lengths")
unit = 'day'

if splittype == 'day':
    diradd = '/Sizes/'
elif splittype == 'MaxSpeeds':
    diradd = '/MaxSpeeds/'

sumlats=0
sumlons=0

data = "TRMM"
version = "Standard"
minGB=0
mapping = 'center'

#minLat = [-45,-10,8,20,15]; maxLat = [-25,10,15,35,30]; minLon = [185,120,150,140,290]; maxLon = [220,160,220,220,340]  #225
minLat = [  8,  2, 22,-40]
maxLat = [  20,  8, 35,-25]
minLon = [155,200,150, 45]
maxLon = [200,280,180, 75]

# In[4]:

#Read in files
fstartyr = 1998
fendyr = 2014
precipClimDir = "/home/disk/eos4/rachel/Obs/TRMM/"
precipClimFile = "TRMM_3B42_1998-2014_annmean.nc"
precipRegFile = 'regress_TRMM_3B42_1999-2014_annmean.nc'
precipReg2File = 'regress_TRMM_3B42_1998-2014_annmean.nc'

olrDir = '/home/disk/eos4/rachel/Obs/OLR/'
olrRegFile = 'regress_olr_1999-2013.nc'

filePrecip = xrayOpen(precipClimDir + precipClimFile)
filePrecipReg =  xrayOpen(precipClimDir + precipRegFile)
filePrecipReg2 =  xrayOpen(precipClimDir + precipReg2File)

fileOlrReg = xrayOpen(olrDir + olrRegFile)

# Events file
dirIn = ('/home/disk/eos4/rachel/EventTracking/FiT_RW_ERA/' + data + '_output/'
            + version + str(fstartyr) + '/proc/' + diradd)

if minGB > 0:
    fileadd = '_min' + str(minGB) + 'GB'
else:
    fileadd = ''

iday = 1
if sumlats > 0:
    fileName = 'DenDirSpd_Map_Ann_Sizes_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + data + "_" + str(fstartyr) + '-' + str(fendyr) + '_' + version + fileadd + '_regrid_' + str(sumlons) + 'lons_' + str(sumlats) + 'lats.nc'
else:
    fileName = 'DenDirSpd_Map_Ann_Sizes_' + str(tbound1[iday]) + '-' + str(tbound2[iday]) + unit + '_' + mapping + '_' + data + "_" + str(fstartyr) + '-' + str(fendyr) + '_' + version + fileadd +  '.nc'

fileIn = xrayOpen(dirIn + fileName)
lats = fileIn['lat']
lons = fileIn['lon']
years = fileIn['years'].sel(years = slice(anstartyr,anendyr))

if lats[-1] < lats[0]:
    print("latitudes are north to south!")

nlats = len(lats)
nlons = len(lons)
nyears = len(years)


# In[ ]:

def ploteachpanel(precipclimin,origlonsin,origlatsin,titlein):
    if np.amax(origlonsin) < 190:
        newprecip2 = precipclimin.sel(longitude=slice(np.amin(origlonsin),0))
        newprecip2.coords['longitude'] = newprecip2.coords['longitude'] + 360
        precipclim = xr.concat((precipclimin.sel(longitude=slice(0,np.amax(origlonsin)+1)),newprecip2),dim='longitude')
    else:
       precipclim = precipclimin

    precipclim.plot(robust=True) # robust=True takes away outliers from colorbar selection process
    plt.title(titlein)
    ax = plt.gca()
    #if origlatsin[0] > 0:
    print 'inverting'
    ax.invert_yaxis()

    for ibox in range(0,len(minLat)):
        ax.add_patch(patches.Rectangle((minLon[ibox],minLat[ibox]),maxLon[ibox]-minLon[ibox],maxLat[ibox]-minLat[ibox],fill=False))
        
        centery = 0.5 *(minLat[ibox] + maxLat[ibox])
        centerx = 0.5 *(minLon[ibox] + maxLon[ibox]) # obviously use a different formula for different shapes

        plt.text(centerx, centery,str(ibox + 1),fontweight='bold',fontsize=14)

def precipplot(precipannin, origlonsin, origlatsin, titlein,nplots,nrows):
    plt.figure(figsize=(15,5*nplots))

    if nplots > 1:
        for iplot in range(0,nplots):
            plt.subplot(nrows,nplots/nrows,iplot+1)
            print precipannin[iplot].shape
            ploteachpanel(precipannin[iplot],origlonsin,origlatsin,titlein[iplot])

    else:
        ploteachpanel(precipannin,origlonsin,origlatsin,titlein)
        


# In[ ]:

#precipplot(filePrecip['pcp'].sel(time=slice(str(anstartyr),str(anendyr))).mean(dim='time'),filePrecip['pcp'].coords['longitude'], 'Annual mean rainfall from TRMM 1998-2014',1,1)


# In[ ]:

preciplons = filePrecipReg['pcp'].coords['longitude']
preciplats = filePrecipReg['pcp'].coords['latitude']

## Now plot regression.
#precipplot([filePrecipReg['pcp'],filePrecipReg2['pcp']],filePrecipReg['pcp'].coords['longitude'],filePrecipReg['pcp'].coords['latitude'], ['Regression of annual mean rainfall from TRMM 1998-2014','Regression of annual mean rainfall from TRMM 1999-2014'],2,2)
#plt.suptitle('Regression of total TRMM precip',size=14)
#plt.tight_layout()
#plt.subplots_adjust(top=0.9)
#plt.savefig(FigDir + 'PrecipRegression_withboxes.eps',format='eps',dpi=120.,facecolor='w',edgecolor='w')

#plt.show()

dataOLR = fileOlrReg['olr']
lonsOLR = fileOlrReg['olr'].coords['lon']
latsOLR = fileOlrReg['olr'].coords['lat']

# need to make these arrays

if np.amax(preciplons) < 190:
    olrtemp1 = preciplons.sel(longitude=slice(np.amin(preciplons),0)) + 360.0
    olrtemp2 = preciplons.sel(longitude=slice(0,np.amax(preciplons)))
    olrnewlons = xr.concat((olrtemp2,olrtemp1),dim='longitude')
else:
    olrnewlons = preciplons

print olrnewlons

OLRnew = dataOLR.sel(lat=preciplats,method='nearest')
OLRnew2 = OLRnew.sel(lon=olrnewlons,method='nearest')

## now plot OLR regression
#precipplot(OLRnew2,OLRnew2.coords['lon'],OLRnew2.coords['lat'],'Regression of annual mean NOAA interpolated OLR from 1999-2013',1,1)
#plt.savefig(FigDir + 'RegridOlrRegression_withboxes.eps',format='eps',dpi=120.,facecolor='w',edgecolor='w')


# Now mask out filePrecipReg['pcp'] based on where the sign is not opposite to
# OLRnew2

precipclimin = filePrecipReg['pcp']
origlonsin = filePrecipReg['pcp'].coords['longitude']
if np.amax(origlonsin) < 190:

    newprecip2 = precipclimin.sel(longitude=slice(np.amin(origlonsin),0))
    newprecip2.coords['longitude'] = newprecip2.coords['longitude'] + 360
    precipclim = xr.concat((
            precipclimin.sel(longitude=slice(0,np.amax(origlonsin)+1)),
                            newprecip2),
                            dim='longitude')

signs = OLRnew2[0,:,:].values * precipclim[0,:,:].values

OLRprecip = precipclim.copy(deep=True)
for ilon in range(0,len(OLRprecip.longitude)):
    for ilat in range(0,len(OLRprecip.latitude)):
        if signs[ilat,ilon] > 0: 
            OLRprecip[0,ilat,ilon] = 0.0

OLRprecip = OLRprecip.rename('mm/yr/yr')

if OLRprecip.long_name == 'precipitation (mm/hr)':
    OLRprecip.values = OLRprecip.values * 24.0 * 365.0
    OLRprecip.long_name = 'precipitation (mm/yr/yr)'
    OLRprecip.units = 'mm/yr/yr'
else:
    exit('i don\'t know the units, and I\'m not going to guess, sorry!')

print OLRprecip 
precipplot(OLRprecip,OLRprecip.coords['longitude'],
           OLRprecip.coords['latitude'],' ',1,1)

#precipplot(OLRprecip,OLRprecip.coords['longitude'],OLRprecip.coords['latitude'],'Regression' +
#    ' of precip masked by NOAA interpolated OLR from 1999-2013',1,1)
plt.savefig(FigDir + 'paper_OLR_Precip_Reg_withboxes.eps',format='eps',dpi=120.,facecolor='w',edgecolor='w')

quit()

#plt.show() 

# ##### x = np.linspace(0, 10)
# line, = plt.plot(x, np.sin(x), '--', linewidth=2)
# 
# dashes = [10, 5, 100, 5]  # 10 points on, 5 off, 100 on, 5 off
# line.set_dashes(dashes)
# 
# plt.show()

# In[ ]:

if splittype == "maxspeed":
    fileadd = "MaxSpeeds_" + str(speedtspan) + "ts_"
elif splittype == "speed":
    fileadd = "Speeds_"
elif splittype == "day":
    fileadd = "Sizes_"
else:
    exit("unexpected splittype")


# In[ ]:

def fileindata(fileIn,varname,fminLat,fmaxLat,fminLon,fmaxLon,anstartyr,anendyr):
    if fminLon <= 180 and fmaxLon > 180:
        temp1 = fileIn[varname].sel(lat=slice(fminLat,fmaxLat),lon=slice(fminLon,180),years=slice(anstartyr,anendyr)).sum(dim=("lon","lat"))
        temp2 = fileIn[varname].sel(lat=slice(fminLat,fmaxLat),lon=slice(-180,fmaxLon-360),years=slice(anstartyr,anendyr)).sum(dim=("lon","lat"))
        temp = temp1 + temp2
    else:
        if fminLon >= 180:
            fminLon = fminLon - 360
            fmaxLon = fmaxLon - 360
        temp = fileIn[varname].sel(lat=slice(fminLat,fmaxLat),lon=slice(fminLon,fmaxLon),years=slice(anstartyr,anendyr)).sum(dim=("lon","lat"))
    return temp

def plotlin_areas2(tbound1,tbound2,plottingin,fminLon,fmaxLon,fminLat,fmaxLat,intitle):
    numcols = 2
    # create numpy array of correct size
    tboundlist = []
    countx = 0
    county = 0
    plt.figure(figsize=(10,5))
    
    for ibound in range(0,nbounds):
    
        if tbound1[ibound] < 0:
            tboundtitle = str(tbound1[ibound]) + '-' + str(tbound2[ibound]) + unit
        else:
            tboundtitle = str(tbound1[ibound]) + '-' + str(tbound2[ibound]) + unit
        
        fileName = 'DenDirSpd_Map_Ann_' + fileadd + tboundtitle + '_' + mapping + '_' + data + "_" + str(fstartyr) + '-' + str(fendyr) + '_' + version +'.nc'
        fileIn = xrayOpen(dirIn + fileName)        
        if plottingin == "PperE":
            Ptemp = fileindata(fileIn,"TPrecip",fminLat,fmaxLat,fminLon,fmaxLon,anstartyr,anendyr)
            Etemp = fileindata(fileIn,"TDensity",fminLat,fmaxLat,fminLon,fmaxLon,anstartyr,anendyr)
            toplot = Ptemp/Etemp
        else:
            toplot = fileindata(fileIn,plottingin,fminLat,fmaxLat,fminLon,fmaxLon,anstartyr,anendyr)
        plt.subplot(2,2,ibound+1)
        toplot.plot.line() #(ax=axes[countx,county])
        plt.title(tboundtitle)

    plt.suptitle(intitle,size=14)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)


def plotlin_accum(tbound1,tbound2,plottingin,minLon,maxLon,minLat,maxLat,intitle,figtitle,ncols,plottype):
    nplots = len(minLon)
    nrows = np.ceil(float(nplots)/float(ncols))
    plt.figure(figsize=(10,5 * nrows))

    contribs = np.zeros([len(tbound1),anendyr-anstartyr+1],float)
    titles = []
    colors = ['r','b','g','k']
    positions = np.arange(0,anendyr - anstartyr,(anendyr-anstartyr)/len(tbound1))

    for ibound in range(0,nbounds):
		if tbound1[ibound] < 0:
			tboundtitle = str(tbound1[ibound]) + '-' + str(tbound2[ibound]) + unit
		else:
			tboundtitle = str(tbound1[ibound]) + '-' + str(tbound2[ibound]) + unit

		titles.append(tboundtitle)


    for iplot in range(0,nplots):
        plt.subplot(nrows,ncols,iplot+1)

        for ibound in range(0,nbounds):
            fileName = 'DenDirSpd_Map_Ann_' + fileadd + titles[ibound] + '_' + mapping + '_' + data + "_" + str(fstartyr) + '-' + str(fendyr) + '_' + version +'.nc'
            fileIn = xrayOpen(dirIn + fileName) 
            
            if plottingin == "PperE":
                Ptemp = fileindata(fileIn,"TPrecip",minLat[iplot],maxLat[iplot],minLon[iplot],maxLon[iplot],anstartyr,anendyr)
                Etemp = fileindata(fileIn,"TDensity",minLat[iplot],maxLat[iplot],minLon[iplot],maxLon[iplot],anstartyr,anendyr)
                contribs[ibound,:] = Ptemp/Etemp
            else:
                contribs[ibound,:] = fileindata(fileIn,plottingin,minLat[iplot],maxLat[iplot],minLon[iplot],maxLon[iplot],anstartyr,anendyr)

        if plottype == 'norm':
            xrind = xr.DataArray(contribs,coords=[('lowbound',tbound1),('years',range(anstartyr, anendyr+1))])
            xrtotal = xrind.sum(dim='lowbound')
            xrtotalnorm  = (xrtotal - xrtotal.mean(dim='years'))/ (xrtotal.max(dim='years') - xrtotal.min(dim='years'))
            xrcontribs = (xrind - xrind.mean(dim='years'))/ (xrind.max(dim='years')-xrind.min(dim='years'))
            miny = -0.8
            maxy = 0.8
            
        elif plottype == 'fraction':
            xrind = xr.DataArray(contribs,coords=[('lowbound',tbound1),('years',range(anstartyr, anendyr+1))])
            total = xrind.sum(dim='lowbound')

            xrcontribs = xrind/total
            miny = 0
            maxy = 0.8
    
        elif plottype == "accum":

            xrind = xr.DataArray(contribs,coords=[('lowbound',tbound1),('years',range(anstartyr, anendyr+1))])
            xrcontribs = xrind.copy(deep=True)
            for ibound in range(0,nbounds):
                xrcontribs[ibound,:] = xrind[0:ibound+1,:].sum(dim='lowbound')
            miny = 0 
            maxy = xrcontribs.max(dim=('years','lowbound'))

    	for ibound in range(0,nbounds):
            regresstemp = stats.linregress(xrcontribs.coords['years'],xrcontribs[ibound,:])
            if regresstemp[3] < pthresh:
                plt.plot(xrcontribs.coords['years'],xrcontribs[ibound,:],colors[ibound],label=titles[ibound])
            else:
                plt.plot(xrcontribs.coords['years'],xrcontribs[ibound,:],colors[ibound],label=titles[ibound],linestyle=':')

            #axes = plt.gca()
            #axes.set_ylim([miny,maxy])
            #plt.ylim((miny,maxy))
            if regresstemp[3] < pthresh:
                plt.text(xrcontribs.coords['years'][positions[ibound]],xrcontribs[ibound,positions[ibound]]+0.05, 'r^2 = ' + '{:5.3f}'.format(regresstemp[2]) + '; p= ' + '{:5.3f}'.format(regresstemp[3]),color=colors[ibound])

        if plottype == 'norm':
            regresstemp = stats.linregress(xrcontribs.coords['years'],xrtotalnorm)
            if regresstemp[3] < pthresh:
                plt.plot(xrcontribs.coords['years'],xrtotalnorm,'m',label='all')
            else:
                plt.plot(xrcontribs.coords['years'],xrtotalnorm,'m',label='all',linestyle=':')

        axes = plt.gca()
        #axes.set_ylim([miny,maxy])
        plt.ylim((miny,maxy))

    	plt.title(intitle + ' Area ' + str(iplot+1) ) #+ ': ' + str(minLon[area]) + '-' + str(maxLon[area]) + 'E; ' + str(minLat[area]) + '-' + str(maxLat[area]) + 'N')

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    lgd = ax.legend(handles, labels, loc=2, bbox_to_anchor=(1.0,1))

    plt.tight_layout()
    plt.savefig(FigDir + figtitle,format='eps',dpi=120.,facecolor='w',edgecolor='w',bbox_extra_artists=(lgd,), bbox_inches='tight')
    #plt.show()

    
    
def plotarea(area):
    loctitle = str(minLon[area]) + '-' + str(maxLon[area]) + 'E; ' + str(minLat[area]) + '-' + str(maxLat[area]) + 'N'    
    plotlin_areas2(tbound1,tbound2,'TDensity',minLon[area],maxLon[area],minLat[area],maxLat[area],'Number of events, ' + loctitle)
    plotlin_areas2(tbound1,tbound2,'TPrecip',minLon[area],maxLon[area],minLat[area],maxLat[area],'Precip in events, ' + loctitle)
    plotlin_areas2(tbound1,tbound2,'PperE',minLon[area],maxLon[area],minLat[area],maxLat[area],'Precip per events, ' + loctitle)


# 
# ***
# 
# 
# ## Area 1:

# In[ ]:

#plotarea(0)


# ***
# 
# ## Area 2

# In[ ]:

#plotarea(1)


# ***
# 
# ## Area 3

# In[ ]:

#plotarea(2)


# ***
# ## Area 4

# In[ ]:

#plotarea(3)


# ***
# ## Area 5

# In[ ]:

#plotarea(4)


# ***
# ## Area 1, single timespan

# In[ ]:

#plotarea(1)


# ***
# ## Area 3 single timespan
# 

# In[ ]:

#plotarea(2)


# 1. Plot precip in events against number of events?
# --- Where there is little trend in precip, there is no trend in precip per event
# 
# 2. For each region, plot line graphs showing the contribution to the total change in precip from the different timespan of events
# 

# In[ ]:

# 
for area in range(0,len(minLat)):
	plotlin_accum(tbound1,tbound2,'TPrecip',minLon,maxLon,minLat,maxLat,'Precip frac, ','FractionTime_Areas.eps',2,'fraction')


# Plot normalized changes in precip from each event
for area in range(0,len(minLat)):
    plotlin_accum(tbound1,tbound2,'TPrecip',minLon,maxLon,minLat,maxLat,'Norm. precip, ','NormalizedPrecip_Areas.eps',2,'norm')


for area in range(0,len(minLat)):
    plotlin_accum(tbound1,tbound2,'TDensity',minLon,maxLon,minLat,maxLat,'Norm. # events, ','NormalizedEvents_Areas.eps',2,'norm')

for area in range(0,len(minLat)):
    plotlin_accum(tbound1,tbound2,'PperE',minLon,maxLon,minLat,maxLat,'Norm. Precip per event, ','NormalizedPperE_Areas.eps',2,'norm')


for area in range(0,len(minLat)):
    plotlin_accum(tbound1,tbound2,'TPrecip',minLon,maxLon,minLat,maxLat,'Accum precip, ','AccumPrecip_Areas.eps',2,'accum')

