# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 2016

@author: rachel
"""

import os, errno
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import Ngl
import xray
from rhwhitepackages.readwrite import XrayOpen
from rhwhitepackages.plots import initcontourplot,setplotrange,initCScontourplot
import argparse
import copy

parser = argparse.ArgumentParser(description="generic ocean plots")
parser.add_argument('--var',metavar='var',type=str,help='variable to be plotted')
parser.add_argument('--lev',metavar='lev',type=int,nargs='?',default=-1,help='level to be plotted')
parser.add_argument('--expN',metavar='exp',type=str,nargs='?',default="CAM4POP_f19g16C_noTopo",help='experiment WITHOUT orography')
parser.add_argument('--res',metavar='res',type=str,nargs='?',default="f19g16",help='resolution')


args = parser.parse_args()
print "here's what I have as arguments: "
print args

var = args.var
ilev = args.lev
expN = args.expN
res = args.res

plotstartlon = 0
plotendlon = 360.0

plotstartlat = -25
plotendlat = 25

plotstartdepth = 0
plotenddepth = 5000

FigDir = "/home/disk/eos4/rachel/Figures/AMOC/Clim/"

rhow = 1000.0	# density of freshwater in kg/m3


if expN == "CAM4POP_f19g16C_noTopo":
	startyr = 211
	endyr = 240
elif expN == "CAM4POP_f19g16C_noR":
	startyr = 411
	endyr = 440
elif expN in ["CAM4POP_f45g37_10All_75S-0S_0", "CAM4POP_f45g37_10All_75S-45S_0", "CAM4POP_f45g37_10At25_55N_1", "CAM4POP_f45g37_10At25_55N_2", "CAM4POP_f45g37_10At25_55N_3", "CAM4POP_f45g37_10At25_55N_4"]:
	startyr = 50
	endyr = 99

nyears = endyr-startyr + 1

if res == "f19g16":
	expM = "CAM4POP_B1850_NCAR"
	FilenameM = "/home/disk/rachel/CESM_outfiles/CAM4POP_B1850_NCAR/ocn/hist/ClimAnn_b40.1850.track1.2deg.003.pop.h.500-529.nc"
	FilenameMA = "/home/disk/rachel/CESM_outfiles/CAM4POP_B1850_NCAR/atm/hist/ClimAnn_b40.1850.track1.2deg.003.cam.h0.500-529.nc"

	FilenameN = "/home/disk/rachel/CESM_outfiles/" + expN + "/ocn/hist/ClimAnn_" + expN + ".pop.h." + str(startyr) + "-" + str(endyr) + ".nc"
	FilenameNA = "/home/disk/rachel/CESM_outfiles/" + expN + "/atm/hist/ClimAnn_" + expN + ".cam.h0." + str(startyr) + "-" + str(endyr) + ".nc"

elif res == "f45g37":
	expM = expN
	expN = "CAM4POP_f45g37"
	FilenameN = "/home/disk/rachel/CESM_outfiles/FluxEXP/" + expN + "/ocn/ClimAvg_" + expN + ".pop.h.0005-0105.nc"
	FilenameNA = "/home/disk/rachel/CESM_outfiles/FluxEXP/" + expN + "/atm/ClimAvg_" + expN + ".cam2.h0.0005-0105.nc"

	FilenameM = "/home/disk/rachel/CESM_outfiles/FluxEXP/" + expM + "/ocn/ClimAvg_" + expM + ".pop.h." + '{:04d}'.format(startyr) + "-" + '{:04d}'.format(endyr) + ".nc"
	FilenameMA = "/home/disk/rachel/CESM_outfiles/FluxEXP/" + expM + "/atm/ClimAvg_" + expM + ".cam2.h0." + '{:04d}'.format(startyr) + "-" + '{:04d}'.format(endyr) + ".nc"
else:
	exit("resolution not defined" + res)

FileM = XrayOpen(FilenameM,False)
FileMA = XrayOpen(FilenameMA,False)

FileN = XrayOpen(FilenameN,False)
FileNA = XrayOpen(FilenameNA,False)

Tlats = FileM['TLAT']
Tlons = FileM['TLONG']
HT = FileM['HT']/100.0 	# convert to m

if var == "AMOC":
	varMGlob = np.squeeze(FileM['MOC'][:,0,:,:,:].sum(dim='moc_comp'))
	varNGlob = np.squeeze(FileN['MOC'][:,0,:,:,:].sum(dim='moc_comp'))

	varMAtl = np.squeeze(FileM['MOC'][:,1,:,:,:].sum(dim='moc_comp'))
        varNAtl = np.squeeze(FileN['MOC'][:,1,:,:,:].sum(dim='moc_comp'))

	varMPac = varMGlob-varMAtl
        varNPac = varNGlob-varNAtl
elif var in ["SS","TS","RHOS"]:
	if var == "SS":
		a1,a2,adiff = 32,36,0.3
		b1,b2,bdiff = -2,2,0.4
		varM = np.squeeze(FileM['SALT'].sel(z_t=500))
		varN = np.squeeze(FileN['SALT'].sel(z_t=500))
		
	elif var == "TS":
		a1,a2,adiff = 0,30,3.0
		b1,b2,bdiff = -4,4,0.8
                varM = np.squeeze(FileM['TEMP'].sel(z_t = 500))
                varN = np.squeeze(FileN['TEMP'].sel(z_t = 500))

	elif var == "RHOS":
		a1,a2,adiff = 1.02,1.03,0.001
		b1,b2,bdiff = -0.002,0.002,0.0004
                varM = np.squeeze(FileM['RHO'].sel(z_t = 500))
                varN = np.squeeze(FileN['RHO'].sel(z_t = 500))	
elif var in ["P-E","PRECT","EVAP","SHFLX"]:
	if var == "P-E":
		a1,a2,adiff = -4.0,4.0,0.8
                b1,b2,bdiff = -2.0,2.0,0.4

		varM = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileMA['PRECT'] - FileMA['QFLX']/rhow)	# convert from m/s to mm/day
                varN = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileNA['PRECT'] - FileNA['QFLX']/rhow)

	elif var in ["PRECT","EVAP"]:
		a1,a2,adiff = 0.0,10.0,1.0
		b1,b2,bdiff = -1.0,1.0,0.2
		if var == "PRECT":
			varM = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileMA['PRECT'])  # convert from m/s to mm/day
			varN = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileNA['PRECT'])
		elif var == "EVAP":
			varM = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileMA['QFLX']/rhow)  # convert from m/s to mm/day
			varN = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileNA['QFLX']/rhow)
	elif var in ["SHFLX"]:
                a1,a2,adiff = -30,30.0,6.0
                b1,b2,bdiff = -10.0,10.0,2.0
		varM = np.squeeze(FileMA[var])  # convert from m/s to mm/day
		varN = np.squeeze(FileNA[var])

	lons = FileMA['lon'].values
	lats = FileMA['lat'].values


	varM = np.where(np.isnan(varM),-9999,varM)
	varN = np.where(np.isnan(varN),-9999,varN)

	varM,lons = Ngl.add_cyclic(varM, lons)
        varN = Ngl.add_cyclic(varN)

elif var in ["PRECTZM"]:
	varM = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileMA['PRECT'].mean(dim='lon'))  # convert from m/s to mm/day
	varN = 1000.0 * 60.0 * 60.0 * 24.0 * np.squeeze(FileNA['PRECT'].mean(dim='lon'))
        lats = FileMA['lat'].values


elif var == "FSNSC":
	a1,a2,adiff = 50,400,35
	b1,b2,bdiff = 0.0,10.0,5.0
	varM = np.squeeze(FileMA['FSNSC'])

        landfrac = np.squeeze(FileMA['LANDFRAC'])
	print landfrac
        lons = FileMA['lon']
        lats = FileMA['lat'].values

	landfrac = Ngl.add_cyclic(landfrac)
        varM,lons = Ngl.add_cyclic(varM, lons)
	print varM.shape
	print lons.shape

	nlats = lats.shape[0]
	nlons = lons.shape[0]
	print nlats,nlons
	varMSO = np.zeros(varM.shape)
	varMSH = np.zeros(varM.shape)
	varMNA = np.zeros(varM.shape)
	varMSO[:,:] = -9999
        varMSH[:,:] = -9999
	varMNA[:,:] = -9999
	for ilat in range(0,nlats):
		for ilon in range(0,nlons):
			if landfrac[ilat,ilon] < 0.5:
				if lats[ilat] < 0 and lats[ilat] > -70:
					varMSH[ilat,ilon] = 10.0
				if lats[ilat] < -45 and lats[ilat] > -70:
					varMSO[ilat,ilon] = 10.0
				if lats[ilat] < 55 and lats[ilat] > 25 and lons[ilon] > 270:
					varMNA[ilat,ilon] = 10.0

        print np.sum(varMSH - varM)
        print np.sum(varMSO - varM)
	print np.sum(varMNA - varM)
else:
	exit(var + " not defined!")

FigTitle = var + "_climmean_" + expM + "-" + expN + '_' + str(nyears)
wkres = Ngl.Resources()
if var == "PRECT":
	wkres.wkColorMap = "precip_diff_12lev"
elif var == "FSNSC":
        wkres.wkColorMap = "WhiteYellowOrangeRed"
else:
	wkres.wkColorMap = "sunshine_diff_12lev"
wks_type = "eps"

wks = Ngl.open_wks(wks_type,FigDir + FigTitle,wkres)


gsres = Ngl.Resources()
gsres.gsLineColor       = "Black"


plot = []

nrows = 0
if var == "AMOC":

	lats = FileM['lat_aux_grid']
	depths = FileM['moc_z']/100.0
	print lats.shape
	print depths.shape
	print varMGlob.shape
	resCS = Ngl.Resources()
	resCS = initCScontourplot(resCS,plotstartlat,plotstartdepth,plotendlat,plotenddepth,lats.values,depths.values)	
	resCS.nglYAxisType = "LinearAxis"
        resCS.nglXAxisType = "LinearAxis"
	resCS.vpWidthF = 0.8
	resCS.vpHeightF = 0.4


        setplotrange(resCS,-20,20,4)
        resCS.tiYAxisString = "depth, m"

	resCS.tiMainFontHeightF = 0.025
        resCS.tiMainString = "Pacific and Indian meridional overturning"
	plot.append(Ngl.contour(wks,varMPac,resCS))
        resCS.tiMainString = "Atlantic meridional overturning"
        plot.append(Ngl.contour(wks,varMAtl,resCS))

        #resCS.tiMainString = "Pacific and Indian meridional overturning"
        #plot.append(Ngl.contour(wks,varNPac,resCS))
        #resCS.tiMainString = "Atlantic meridional overturning"
        #plot.append(Ngl.contour(wks,varNAtl,resCS))

	resCS.tiMainString = "Change in Pacific and Indian meridional overturning"
	resCS.tiYAxisString = "depth, m"

        setplotrange(resCS,-5,5,1)

	temp = Ngl.contour(wks,varMPac-varNPac,resCS)
        Ngl.add_polyline(wks,temp,np.array([0,0]),np.array([0,5000]),gsres)
	plot.append(temp)

        resCS.tiMainString = "Change in Atlantic meridional overturning"
        setplotrange(resCS,-5,5,1)
        
	temp = (Ngl.contour(wks,varMAtl-varNAtl,resCS))
        Ngl.add_polyline(wks,temp,np.array([0,0]),np.array([0,5000]),gsres)
        plot.append(temp)
	nrows = 4

elif var in ["P-E","PRECT","EVAP","SHFLX"]:
	resMP = Ngl.Resources()
        resMP = initcontourplot(resMP,plotstartlat,plotstartlon,plotendlat,plotendlon,lats.values,lons)
	resMP.mpCenterLonF = 0. 	
	resMP.mpGridAndLimbOn = False
	setplotrange(resMP,a1,a2,adiff)
        plot.append(Ngl.contour_map(wks,varM,resMP))
        plot.append(Ngl.contour_map(wks,varN,resMP))
        setplotrange(resMP,b1,b2,bdiff)
        plot.append(Ngl.contour_map(wks,varM-varN,resMP))
	nrows = 3
elif var in ["PRECTZM"]:
	resLP = Ngl.Resources()
	resLP.nglDraw = False
	resLP.nglFrame = False
	resLP.xyMarkLineMode = "Lines"
        resLP.tiMainString = expM
	resLP.xyLineColor = "blue"
	resLP.tiMainFontHeightF = 0.012
	resLP.tiYAxisString = "Zonal mean precip, mm/day"
        resLP.trXMinF = plotstartlat
        resLP.trXMaxF = plotendlat
	temp = Ngl.xy(wks,lats,varM.values,resLP)
        plot.append(temp)
        resLP.tiMainString = expN
        temp = Ngl.xy(wks,lats,varN.values,resLP)
	plot.append(Ngl.xy(wks,lats,varN.values,resLP))

        a1 = 0.5
        resLP.trYMinF = -a1
        resLP.trYMaxF = a1
        resLP.tiMainString = expM + " - " + expN

        temp = Ngl.xy(wks,lats,varM.values-varN.values,resLP)
        Ngl.add_polyline(wks,temp,np.array([0,0]),np.array([-a1,a1]),gsres)
        plot.append(temp)

	nrows = 3

elif var in ["FSNSC"]:
        resMP = Ngl.Resources()
        resMP = initcontourplot(resMP,plotstartlat,plotstartlon,plotendlat,plotendlon,lats,lons)
        #resMP.cnMissingValFillColor = "white"
        #resMP.sfMissingValueV = -9999
        resMP.cnFillMode = "AreaFill"
        resMP.mpCenterLonF = 0.
        resMP.mpGridAndLimbOn = False
        resMP.cnConstFEnableFill = True
	setplotrange(resMP,a1,a2,adiff)
        plot.append(Ngl.contour_map(wks,varM,resMP))
        plot.append(Ngl.contour_map(wks,varMSH,resMP))
        setplotrange(resMP,b1,b2,bdiff)
        plot.append(Ngl.contour_map(wks,varMSO,resMP))
	plot.append(Ngl.contour_map(wks,varMSH,resMP))
        plot.append(Ngl.contour_map(wks,varMNA,resMP))

#---Attach blank plot with special labels to map plot, and return
	resBMP = Ngl.Resources()
        resBMP.mpProjection = "CylindricalEquidistant" # Change the map projection.
        resBMP.mpCenterLonF = 0.           # Rotate the projection.
        resBMP.mpGridAndLimbOn = False
        resBMP.pmLabelBarDisplayMode = "Always"
        resBMP.lbPerimOn = False
	resBMP.nglDraw  = False
        resBMP.nglFrame = False
        resBMP.mpOutlineBoundarySets = "AllBoundaries"
        resBMP.mpLimitMode = "LatLon"    # Limit the map view.
        resBMP.mpMinLonF = plotstartlon
        resBMP.mpMaxLonF = plotendlon
        resBMP.mpMinLatF = plotstartlat
        resBMP.mpMaxLatF = plotendlat
	plot.append(Ngl.map(wks,resBMP))


else:
	resMP = Ngl.Resources()
	#resMP = initcontourplot(resMP,plotstartlat,plotstartlon,plotendlat,plotendlon,Tlats.values,Tlons.values)
        #resMP.nglYAxisType = "LinearAxis"
        #resMP.nglXAxisType = "LinearAxis"
	resMP.cnFillOn = True
	resMP.cnFillMode = "AreaFill"
	resMP.gsnAddCyclic = True
	resMP.MissingValFillColor = "white"
	resMP.sfMissingValueV = -9999
	resMP.mpGridMaskMode        = "MaskLand"
        resMP.nglDraw  = False
        resMP.nglFrame = False

	resMP.sfXArray = Tlons.values
	resMP.sfYArray = Tlats.values

        resMP.cnLineLabelsOn = False
        resMP.cnLinesOn = False
        resMP.cnLevelSelectionMode = "ManualLevels"

	contour = Ngl.contour(wks,varM,resMP)

        resMP.mpProjection = "CylindricalEquidistant" # Change the map projection.
        resMP.mpCenterLonF = 0.           # Rotate the projection.
        #resMP.mpFillOn     = True           # Turn on map fill.
        resMP.mpLimitMode = "LatLon"    # Limit the map view.
        resMP.mpMinLonF = Ngl.get_float(contour.sffield,"sfXCActualStartF")
	resMP.mpMaxLonF = Ngl.get_float(contour.sffield,"sfXCActualEndF")
	resMP.mpMinLatF = Ngl.get_float(contour.sffield,"sfYCActualStartF")
	resMP.mpMaxLatF = Ngl.get_float(contour.sffield,"sfYCActualEndF")
        #resMP.mpOutlineBoundarySets = "AllBoundaries"

        #resMP.sfYArray = Tlats.values
        #resMP.sfXArray = Tlons.values




	setplotrange(resMP,a1,a2,adiff)
        plot.append(Ngl.contour_map(wks,varM,resMP))
        plot.append(Ngl.contour_map(wks,varN,resMP))
        setplotrange(resMP,b1,b2,bdiff)
        plot.append(Ngl.contour_map(wks,varM-varN,resMP))


textres = Ngl.Resources()
textres.txFontHeightF = 0.012
Ngl.text_ndc(wks,var + '; ' +  expM  + ' - ' + expN ,0.5,0.9,textres)

panelres = Ngl.Resources()
panelres.nglPanelLabelBar                 = False     # Turn on panel labelbar
panelres.nglPaperOrientation = "Auto"
panelres.nglPanelTop = 0.9
panelres.nglPanelBottom = 0.1
panelres.nglPanelYWhiteSpacePercent = 10.
#panelres.nglPanelXWhiteSpacePercent = 5.

Ngl.panel(wks,plot,[nrows,1],panelres)

Ngl.end()

