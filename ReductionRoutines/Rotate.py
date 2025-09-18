#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 10:06:13 2023

@author: marah
"""

import numpy as np
import glob
import os
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from scipy.ndimage import shift
from scipy.ndimage import rotate


def rotate(lm_offset,pathname,savepath):

    pathname = pathname
    savepath = savepath
    
    
    os.chdir(pathname)
    sciencelist= glob.glob('*fits')
    sciencelist.sort()
    print(sciencelist)
    
    '''this line will make a new folder called "rotated" if one doesnt already exist'''
    if not os.path.exists('rotated'): os.mkdir('rotated')
    '''Loop through each file'''
    for fits in sciencelist:
        scidat = pf.getdata(fits)
        scihdr = pf.getheader(fits)
        parang = scihdr['LBT_PARA']
        angle = float(parang)
        '''We need to set an NaN values to 1 or rotate wont work'''
        scidat[np.isnan(scidat)] = 1
        '''Rotate the image'''
        finalangle = angle*-1 + lm_offset
        image = rotate(scidat,(finalangle),reshape=False)
        final = image
        '''This will save the file in the new folder'''
        pf.writeto(savepath + fits,final,scihdr,overwrite=True)
