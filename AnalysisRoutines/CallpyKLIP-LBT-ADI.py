#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 18:55:55 2024

@author: marah
"""
import os
import numpy as np
import astropy.io.fits as fits
import scipy.ndimage
import pyklip.klip
import pyklip.parallelized
import pyklip.instruments.Instrument as Instrument
import pyklip.fakes
import pyklip.kpp.metrics.crossCorr
import pyklip.kpp.stat.statPerPix_utils
import multiprocessing as mp
import glob
import matplotlib.pyplot as plt

import os

def pyklip_ADI(pathname, outputdir, prefix, ann, subs, numbs, maxnumbs, mode, star_centx, star_centy):
    mp.set_start_method("fork",force=True)
    pathname = pathname
    outputdir = outputdir
    prefix = prefix
    ann = ann
    subs = subs
    numbs = numbs
    maxnumbs = maxnumbs
    mode = "ADI"
    
    star_centx = star_centx
    star_centy = star_centy
    
    os.chdir(pathname)
    fitslst = glob.glob('*.fits')
    fitslst.sort()
    print(fitslst)
    
    # load in the calibration frame to calibrate the brightness of any sources with respect to the star
    # with fits.open(file_unsat) as hdulist:
    #     calib_frame = hdulist[0].data # image of the unsaturated star for photometric calibration
    #     # crop it down because we don't need it so big
    #     calib_orig_cent = (starcentx, starcenty)
    #     stamprad = 15
    #     calib_frame = calib_frame[calib_orig_cent[0]-stamprad:calib_orig_cent[0]+stamprad+1, calib_orig_cent[1]-stamprad:calib_orig_cent[1]+stamprad+1]
    #     calib_exptime = hdulist[0].header['ITIME']
    
    # the location of the star in the calibration frame
    calib_centx = star_centx
    calib_centy = star_centy
    
    
    # plt.figure()
    # plt.imshow(img_cube[0], cmap="inferno")
    # plt.xlim([star_centx-150, star_centx+150])
    # plt.ylim([star_centy-150, star_centy+150])
    # plt.title("Science Frame")
    
    # plt.figure()
    # plt.imshow(calib_frame, cmap="inferno")
    # plt.title("Calibration Frame")
    
    centers = np.array([[star_centx, star_centy] for _ in range(0,len(fitslst))])
    
    data = [fits.getdata(f) for f in fitslst]
    
    print('defining dataset')
    dataset = Instrument.GenericData(data, centers, filenames=fitslst)
    dataset.IWA = 5
    dataset.OWA = 200

    print('starting klip')
    pyklip.parallelized.klip_dataset(dataset, outputdir=outputdir, fileprefix=prefix, annuli=ann,
                                 subsections=subs, numbasis=numbs, maxnumbasis=maxnumbs, mode=mode,
                                 movement=0,numthreads=4)
    
    
    


