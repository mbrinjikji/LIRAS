#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:33:52 2023

@author: marah
"""

import os
import numpy as np
import astropy.io.fits as pf
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore',category=AstropyWarning)
import matplotlib.pyplot as plt
import glob
from scipy import ndimage

def dualsplit(pathname):
    
    pathname = pathname

    os.chdir(pathname)
    sciencelist= glob.glob('*fits')
    sciencelist.sort()
    print(sciencelist)
    
    if not os.path.exists('SX'): os.mkdir('SX')
    if not os.path.exists('DX'): os.mkdir('DX')
    for fits in sciencelist:
        scihdr = pf.getheader(fits)
        science = pf.getdata(fits)
        # sx = science[:,0:1024]
        # dx = science[:,1024:2048]
        sx = science[:1024,:]
        dx = science[1024:,:]
        pf.writeto(pathname + 'SX/' + fits,sx,scihdr,overwrite=True,output_verify='ignore')
        pf.writeto(pathname + 'DX/' + fits,dx,scihdr,overwrite=True,output_verify='ignore')
