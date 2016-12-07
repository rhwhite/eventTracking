"""
Created on Jun 4 2016

@author: rhwhite
"""

import os, errno
from netCDF4 import Dataset
import netCDF4
import numpy as np
import datetime as dt
import pandas
import xray
import Ngl
from scipy import stats

dayinmon = [31,28,31,30,31,30,31,31,30,31,30,31]
dayinyear = 365


def calcann(tstart,tend,var,annvar,years)
	# First check if the data is monthly or daily
	yeartest1 = years[0]
	yeartest2 = years[11]
	yeartest3 = years[26]
	if yeartest1 != yeartest2 and yeartest2 != yeartest3:
		res = "month"
	else:
		yeartest4 = years[366]
		if yeartest1 != yeartest4:
			res = "day"
		else:
			print "resolution less than 1 day - need to sort this!"
			res = "somehour"

	# Check whether first year is a full year:
        yearstart = years[tstart]
	if res == "month":
		yearstart2 = years[tstart+11]
	elif res == "day":
		yearstart2 = years[tstart + 364]
	else:
		print "don't have a check for sub-daily data yet, just going to assume 6 hourly"
		yearstart2 = years[tstart + 365*4]	
	if yearstart == yearstart2:
		tstart = tstart
	else:
		# Find beginning of new year
		while year == yearstart:
			tstart = tstart + 1
			year = years[tstart]
			
	startyear = year[tstart]
	yearnow = -999
	inyearcount = 0
	for icount in range(tstart,tend):
		if years[icount] == yearnow:
			if res == month:
			        # If monthly need to multiply by dayinmon
				temp = temp + var[icount,:,:] * dayinmon[inyearcount]
			else:
				temp = temp + var[icount,:,:]

			inyearcount = inyearcount + 1
		else:

			if res == month:
				if inyearcount != 12:
					print inyearcount
					print "not 12 months in year"
					break
				else:
					varann[years[yearnow - startyear,:,:] = temp / dayinyear

			elif res == day:
				if inyearcount < 364 or inyearcount > 365:
					print inyearcount
					print "not 365 days per year!"
					break
				else:
					varann[years[yearnow - startyear,:,:] = temp / dayinyear
			else:
				varann[years[yearnow - startyear,:,:] = temp / inyearcount

			yearnow = years[icount]
			inyearcount = 1
			if res == month:
				temp = var[icount,:,:] * dayinmon[inyearcount]
			else:
				temp = var[icount,:,:]




def calcseas(tstart,tend,var,seasvar,years,seasons)
        # First check if the data is monthly or daily
        yeartest1 = years[0]
        yeartest2 = years[11]
        yeartest3 = years[26]
        if yeartest1 != yeartest2 and yeartest2 != yeartest3:
                res = "month"
        else:
                yeartest4 = years[366]
                if yeartest1 != yeartest4:
                        res = "day"
                else:
                        print "resolution less than 1 day - need to sort this!"
                        res = "somehour"

        	elif months[icount] == 3:
        	        seasprecip[iyear,1,:,:] = (Precip[icount,:,:] * 31.0 + Precip[icount + 1,:,:] * 30.0 + Precip[icount+2,:,:] * 31.0)/92.0
        	elif months[icount] == 6:
        	        seasprecip[iyear,2,:,:] = (Precip[icount,:,:] * 30.0 + Precip[icount + 1,:,:] * 31.0 + Precip[icount+2,:,:] * 31.0)/92.0
        	elif months[icount] == 9:
        	        seasprecip[iyear,3,:,:] = (Precip[icount,:,:] * 30.0 + Precip[icount + 1,:,:] * 31.0 + Precip[icount+2,:,:] * 30.0)/91.0
        	elif months[icount] == 12:
        	        seasprecip[iyear,0,:,:] = (Precip[icount,:,:] * 31.0 + Precip[icount + 1,:,:] * 31.0 + Precip[icount+2,:,:] * 28.0)/90.0
	



