#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 23:03:30 2023

@author: marah
"""

import os
import numpy as np
import astropy.io.fits as pf
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore',category=AstropyWarning)
from scipy.ndimage import shift
import itertools
import glob

    
def reduce_noflats(darkpath,bpath,sciencepath):
    
    darkpath = darkpath ##path to median dark
    
    bpath = bpath ##path to badpixelmap
    
    sciencepath = sciencepath ##path to raw science files
    
    '''load in bpmask'''
    bpmask = pf.getdata(bpath)
    bpmask = bpmask[:,:]
    bp_ind=np.where(bpmask > 0)
    	
    '''load in median dark from Dark.py'''
    dark = pf.getdata(darkpath)
    dark = dark[:,:]
    
    def fixBadPix(im,NOneSidedShifts,Thres,badpixelmap):
        ShiftList = range(-NOneSidedShifts, NOneSidedShifts+1)
        ShiftTable = list(itertools.product(ShiftList,ShiftList))
        tmp = []
        for i in range(len(ShiftTable)):
            tmp.append(shift(im,ShiftTable[i],order=3))
        std_im,med_im = np.std(tmp,axis=0), np.median(tmp,axis=0)
        test = abs(im-med_im) - Thres*std_im
        bad = np.where(test > 0.)
    
        fim = np.copy(im)
        fim[bad] = med_im[bad]
        bad2 = np.where(badpixelmap==1)
        
        fim[bad2] = med_im[bad2]
        
        print('Done matching badpixelmap.')
        return fim
    
    
    os.chdir(sciencepath)
    sciencelist= glob.glob('*fits')
    sciencelist.sort()
    print(sciencelist)
    scihdr = [pf.getheader(fits) for fits in sciencelist]
    if not os.path.exists('reduced'): os.mkdir('reduced')
    for fits in sciencelist:
        science = pf.getdata(fits)[0]
        scihdr = pf.getheader(fits)
        print(fits)
        tims = (science-dark)
        fims = fixBadPix(tims,3,5,bpmask)
        fims[fims>10000] = 0.0
        pf.writeto(sciencepath + '/reduced/' + fits,fims,header=scihdr,overwrite=True)
    

def reduce_withflats(darkpath,flatpath,bpath,sciencepath):
    
    darkpath = darkpath ##path to median dark
    
    flatpath = flatpath ###path to median flat
    
    bpath = bpath ##path to badpixelmap
    
    sciencepath = sciencepath ##path to raw science files



    '''load in median flat, badpixel corrected and scaled, from Flat.py'''
    flat = pf.getdata(flatpath)
    flat = flat[:,:]
    flat[flat<1E-1] = 1.
    	# --------------------------------------- #
    	# ------ Step 1b: Flat BPmask -------- #
    '''load in bpmask'''
    bpmask = pf.getdata(bpath)
    bpmask = bpmask[:,:]
    bp_ind=np.where(bpmask > 0)
    	
    '''load in median dark from Dark.py'''
    dark = pf.getdata(darkpath)
    dark = dark[:,:]
    
    def fixBadPix(im,NOneSidedShifts,Thres,badpixelmap):
        ShiftList = range(-NOneSidedShifts, NOneSidedShifts+1)
        ShiftTable = list(itertools.product(ShiftList,ShiftList))
        tmp = []
        for i in range(len(ShiftTable)):
            tmp.append(shift(im,ShiftTable[i],order=3))
        std_im,med_im = np.std(tmp,axis=0), np.median(tmp,axis=0)
        test = abs(im-med_im) - Thres*std_im
        bad = np.where(test > 0.)
    
        fim = np.copy(im)
        fim[bad] = med_im[bad]
        bad2 = np.where(badpixelmap==1)
        
        fim[bad2] = med_im[bad2]
        
        print('Done matching badpixelmap.')
        return fim
    
    
    os.chdir(sciencepath)
    sciencelist= glob.glob('*fits')
    sciencelist.sort()
    print(sciencelist)
    scihdr = [pf.getheader(fits) for fits in sciencelist]
    if not os.path.exists('reduced'): os.mkdir('reduced')
    for fits in sciencelist:
        science = pf.getdata(fits)[0]
        scihdr = pf.getheader(fits)
        print(fits)
        tims = (science-dark)/flat
        fims = fixBadPix(tims,3,5,bpmask)
        fims[fims>10000] = 0.0
        pf.writeto(sciencepath + '/reduced/' + fits,fims,header=scihdr,overwrite=True)
    
