#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 20:49:29 2024

@author: marah
"""

from matplotlib import pyplot as plt
from scipy import optimize as op
from astropy.io import fits
from scipy import ndimage
import numpy as np
from scipy.ndimage import gaussian_filter
import emcee
from multiprocessing import Pool
import multiprocessing as mp
from corner import corner
import pandas as pd


def photometry_psffit(filepathname,psf_xy,im_xy,pad,n_walkers,n_burnin, n_iter):
    
    filepathname = filepathname
    
    psf_xy = psf_xy

    im_xy = im_xy

    pad = pad

    n_walkers = n_walkers
    
    n_burnin = n_burnin
    
    n_iter = n_iter
    
    mp.set_start_method("fork", force=True)
    
    def lnprior(theta):
        xs, ys, f = theta
        if 0<=xs<=5 and 0<=ys<=5 and -5<=f<=-3:
            return 0
        
        return -1e300
    
    def lnlike(theta, data):
        xs, ys, f = theta
        
        psf, im, pad = data
        
        s = im.shape
        x, y = np.meshgrid(np.arange(s[0], dtype=np.float64), np.arange(s[1], dtype=np.float64))
        d = np.sqrt((x-pad)**2 + (y-pad)**2)
        indx = np.where(d <= 5)
        mask = np.full(im.shape, np.nan)
        mask[indx] = 1
        
        data_ = (psf, im, mask)
        
        my_log_like = lnChiSq(theta, data_)
        
        return my_log_like
    
    def lnChiSq(theta, data):
        xs, ys, f = theta
        
        psf, im, mask = data
        
        s = psf.shape
        x, y = np.meshgrid(np.arange(s[0], dtype=np.float64), np.arange(s[1], dtype=np.float64))
        psf = ndimage.map_coordinates(psf, (y-ys, x-xs), cval=np.nan) * (10**f)
        resid = (im - psf) * mask
        
        #n = np.count_nonzero(~np.isnan(resid))
        
        return -0.5 * np.nansum(resid**2.0)
    
    def lnprob(theta, data):
    
        ll = lnlike(theta, data)
        
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -1e300, theta
    
        if not np.isfinite(ll):
            return -1e300, theta
    
        my_log_prob = lp + ll
    
        return my_log_prob, theta
    
    psf_path = filepathname #oct 2024 K
    
    im_path = filepathname 
    
    psf = fits.getdata(psf_path)
    im = fits.getdata(im_path)
    
    psf_hdr = fits.getheader(psf_path)
    psf /= (psf_hdr['ITIME']*psf_hdr['NCOADDS'])
    
    im_hdr = fits.getheader(im_path)
    im /= (im_hdr['ITIME']*im_hdr['NCOADDS'])
    
    im_mask = gaussian_filter((im),10)
    im -= im_mask
    
    
    # Extract stamp from full image
    
    psf_stamp = psf[psf_xy[1]-pad:psf_xy[1]+pad, psf_xy[0]-pad:psf_xy[0]+pad]
    im_stamp = im[im_xy[1]-pad:im_xy[1]+pad, im_xy[0]-pad:im_xy[0]+pad]
    
    data = (psf_stamp, im_stamp, pad)
    
    args = (data,)
    
    f_guess = np.log10(np.nanmax(im_stamp)/np.nanmax(psf_stamp))
    
    lower_guess = (0, 0, f_guess-0.5)
    upper_guess = (2, 2, f_guess+0.5)
    
    pos0 = np.array([
        np.random.uniform(lg, ug, n_walkers)
        for lg, ug in zip(lower_guess, upper_guess)
    ]).T
    
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(
            n_walkers,
            len(lower_guess),
            lnprob,
            args=(*args,),
            pool=pool
        )
        
        # Run burn-in stage
        if n_burnin > 0:
            pos, prob, state, blobs = sampler.run_mcmc(pos0, n_burnin, progress=True)
            sampler.reset()
        else:
            pos = pos0
    
        # Run sampling stage
        pos, prob, state, blobs = sampler.run_mcmc(pos, n_iter, progress=True)
    
    #%%
    
    chains = pd.DataFrame(sampler.get_blobs(flat=True), columns=['xs','ys','f'])
    chains = chains.dropna()
    
    chains['delta_mag'] = -2.5*np.log10(10**chains['f'])
    
    fig = corner(chains[['delta_mag','xs','ys']], 
        show_titles=True, levels=(1-np.exp(-2), 1-np.exp(-0.5)))
    
    plt.show()




