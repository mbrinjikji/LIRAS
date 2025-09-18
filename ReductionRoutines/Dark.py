# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import glob
import scipy.ndimage
import scipy.signal
import astropy.io.fits as pf
import os



def dark(path,cube,median,frame):
    
    path = path ###path to dark data

    cube = cube ###True if the image is a cube, False if not

    median = median ###True if you want to median the cube, False if not?

    frame = frame ###frame # in cube you want to use

    
    os.chdir(path)
    darklist = glob.glob('*fits')
    darklist.sort()
    print(darklist)
    
    
    n=len(darklist)
    
    darkheader = pf.getheader(darklist[0])
    
    
    
    '''Read in darks and Median Combine'''
    imsize1=darkheader['NAXIS1']
    imsize2=darkheader['NAXIS2']
    darkarray=np.zeros((imsize2,imsize2,n),dtype=np.float32)
    for i in range(0,n-1):
        print(i)
        fits=pf.open(path + darklist[i])
        im=fits[0].data
        if cube == True:
            if median == True:
                if len(im.shape) > 1: #check for data cubes
                      assert not np.any(np.isnan(im))
                      im = np.median(im,axis=0) #if data cube, then median frames
            if median == False:
                im = im[frame,:,:]
        if cube == False:
            im = im[0,:,:]
        darkarray[:,:,i]=im
    
    med_dark=np.median(darkarray,axis=2)
    
    if not os.path.exists(path+'/reduced/'): os.mkdir(path+'/reduced/')
    darkname = path+'/reduced/median_dark.fits'
    pf.writeto(darkname,med_dark,darkheader,overwrite=True,output_verify='ignore')
    
    '''Bad Pixel Map Creation'''
    
    dark_sm = scipy.signal.medfilt(med_dark, 15)
    std = np.std(med_dark[100:924,100:1000])
    bad = (np.abs(med_dark - dark_sm) / std) > 3.0
    bpmask = bad.astype(int)
    med_dark_cln = med_dark.copy()
    med_dark_cln[bad] = dark_sm[bad]
    
    # badflat = np.ones(med_dark.shape)
    # badflat[bad]=0
    # med_dark_cln[badflat==0] = med_dark[badflat==0]
    
    bpmaskname = path+'/reduced/bpmask_dark.fits'
    pf.writeto(bpmaskname,bpmask,darkheader,overwrite=True,output_verify='ignore')
    med_darkname = path+'/reduced/median_dark_cln.fits'
    pf.writeto(med_darkname,med_dark_cln,darkheader,overwrite=True,output_verify='ignore')


