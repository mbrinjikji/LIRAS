# -*- coding: utf-8 -*-
import numpy as np
import scipy.ndimage
import scipy.signal
import astropy.io.fits as pf
import os
import glob


def flat(flatpath,flatdarkpath):
    
    flatpath = flatpath ##path to flats
    
    flatdarkpath = flatdarkpath ##path to flats for darks

    # ------ Step 2: Subtract median Dark & median stack the Flats -------- #
    os.chdir(flatpath)
    flatlist = glob.glob('*fits')
    flatlist.sort()
    print(flatlist)
    
    n=len(flatlist)
    
    flatheader = pf.getheader(flatlist[0])
    
    
    # # ------ Step 1: Load median dark for flats------- #
    # '''read in the median dark file created by Dark.py for the flats'''
    med_darkname = flatdarkpath
    med_dark_fits = pf.open(med_darkname)
    med_dark_im = med_dark_fits[0].data
    # # --------------------------------------- #
    
    imsize1=flatheader['NAXIS1']
    imsize2=flatheader['NAXIS2']
    
    flatarray=np.zeros((imsize1,imsize2,n),dtype=np.float32)
    for i in range(0,n-1):
        print(i)
        fits=pf.open(flatpath+flatlist[i])
        im=fits[0].data
        im = im[1,:,:]
        im = im - med_dark_im
        flatarray[:,:,i]=im
        
    med_flat=np.median(flatarray,axis=2)
    
    
    if not os.path.exists(flatpath+'/reduced/'): os.mkdir(flatpath+'/reduced/')
    flatname = flatpath+'/reduced/median_flat.fits'
    pf.writeto(flatname,med_flat,flatheader,overwrite=True,output_verify='ignore')
    
    # ------ step 2b: Bad Pixel mask from Median dark -------- #
    ## http://www.ster.kuleuven.be/~pieterd/python/html/core/numpy_scipy.html sigma clipping
    flat_sm = scipy.signal.medfilt(med_flat, 15)
    std = np.std(med_flat[100:924,100:1000])
    bad = np.abs(med_flat - flat_sm) / std > 3.0
    bpmask = bad.astype(int)
    med_flat_cln = med_flat.copy()
    med_flat_cln[bad] = flat_sm[bad]
    
    
    bpmaskname = flatpath+'/reduced/bpmask_flat.fits'
    pf.writeto(bpmaskname,bpmask,flatheader,overwrite=True,output_verify='ignore')
    med_flatname = flatpath+'/reduced/median_flat_cln.fits'
    pf.writeto(med_flatname,med_flat_cln,flatheader,overwrite=True,output_verify='ignore')
    
    # # --------------------------------------- #
    
    # ------ Step 3: Normalize Flat -------- #
    ind=np.where(med_flat_cln == 0)
    if np.sum(ind) > 0:	
     	med_flat_cln[ind] = 0.001
     
    flatscalefactor = np.median(med_flat_cln)
    scl_flat = med_flat_cln/flatscalefactor
    # scl_flat[scl_flat>1.1] = 1.
    # scl_flat[scl_flat<0.9] = 1.
    
    
    scl_flatname = flatpath+'/reduced/median_flat_cln_scl.fits'
    pf.writeto(scl_flatname,scl_flat,flatheader,overwrite=True,output_verify='ignore')



