#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 11:47:39 2019
Start-up file for Spyder
@author: sthalabard
"""

import numpy as np
from math import *
import scipy as scp
from scipy import fftpack as fft
from scipy import io, integrate

from matplotlib import pyplot as plt
from matplotlib.pyplot import *
from matplotlib import animation,cm,rc,colors

from mpl_toolkits.mplot3d import Axes3D

from time import time as ctime
import os,glob,subprocess
import warnings; warnings.simplefilter('ignore')


##########
from IPython.display import display,HTML,Image
display(HTML("<style>.container { width:95% !important; }</style>"))

import sys   
#print(sys.executable)
#%%
def newfig(a=1,b=1,figheight=7,aspectratio=0.9,**kwargs):
    rc('legend', frameon=False,fancybox=False,shadow=False,fontsize=14,loc='best')
    rc('lines', linewidth=1)
    font = {'family':'serif','size':26}
    rc('font',**font)
    rc('text', usetex=True)
    rc('xtick',labelsize=32)
    rc('ytick',labelsize=32)
    rc('savefig',format='pdf')
    return plt.subplots(a,b,figsize=(b*figheight,a*figheight*aspectratio),**kwargs)

class dic2struc:
    def __init__(self, **entries):
        self.__dict__.update(entries)
