
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
#from scipy import stats
from rhwhitepackages.readwrite import shiftlons
from rhwhitepackages.readwrite import xrayOpen
from rhwhitepackages.stats import regressmaps

# plotting
# import matplotlib
import xray.plot as xplt
import matplotlib.pyplot as plt
import matplotlib.patches as patches


# In[3]:

x = np.linspace(0, 10)
line, = plt.plot(x, np.sin(x), '--', linewidth=2)

dashes = [10, 5, 100, 5]  # 10 points on, 5 off, 100 on, 5 off
line.set_dashes(dashes)

plt.show()

