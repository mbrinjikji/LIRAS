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
from astropy.stats import SigmaClip,sigma_clip
from photutils.background import Background2D, MedianBackground
from scipy import ndimage


def skysub(pathname,skypathname,filemidpoint):

    pathname = pathname
    
    skypathname = skypathname
    
    filemidpoint=filemidpoint

    os.chdir(pathname)
    sciencelist= glob.glob('*fits')
    sciencelist.sort()
    print(sciencelist)
    
    scihdr = [pf.getheader(fits) for fits in sciencelist]
    
    alist = []
    blist = []
    
    # for i in range(0,len(sciencelist)):
    #     if scihdr[i]['FLAG'] == 'NOD_A':
    #         alist.append(sciencelist[i])
    #     if scihdr[i]['FLAG'] == 'NOD_B':
    #         blist.append(sciencelist[i])
    
    
    alist = sciencelist[0:filemidpoint]
    blist = sciencelist[filemidpoint:]

    
    '''Sky BG creation'''
    
    a_science = [pf.getdata(fits) for fits in alist]    
    b_science = [pf.getdata(fits) for fits in blist]  
    
    a = np.array(a_science)
    b = np.array(b_science)
    
    a_median = np.median(a,axis=0)
    b_median = np.median(b,axis=0)
    
    a_data = a_median
    b_data = b_median
    '''The background for each nod position is going to be the part of the image where there is no observations'''
    
    imsize1=int(scihdr[0]['NAXIS1']/2)
    imsize2=int(scihdr[0]['NAXIS2']/2)
    
    
    a = a_data[:,imsize1:]
    b = b_data[:,:imsize1]
    
    '''combine the images, make a directory to save, and save the file'''
    final = np.concatenate((b,a),axis=1)
    if not os.path.exists(pathname+'/sky/'): os.mkdir(pathname+'/sky/')
    pf.writeto(skypathname,final,overwrite=True,output_verify='ignore')
    
    '''this is the sky subtraction step -- only run after you made the sky bg file'''
    '''Sky Subtraction'''
    
    skydata = pf.getdata(skypathname)
    
    if not os.path.exists('SkySub'): os.mkdir('SkySub')
    for fits in sciencelist:
        scihdr = pf.getheader(fits)
        science = pf.getdata(fits)
        flag = scihdr['FLAG']    
        if 'NOD_A' in flag:
            skysub = science - skydata
            pf.writeto(pathname + 'SkySub/' + fits[:-5] + '_skysub.fits',skysub,scihdr,overwrite=True,output_verify='ignore')
        else:
            skysub = science - skydata
            pf.writeto(pathname + 'SkySub/' + fits[:-5] + '_skysub.fits',skysub,scihdr,overwrite=True,output_verify='ignore')
        
    
