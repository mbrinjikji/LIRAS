#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 11:04:19 2024

@author: marah
"""
import glob, os
import astropy.io.fits as pf
import numpy as np
from photutils.datasets import load_star_image
from photutils.aperture import aperture_photometry, CircularAperture, CircularAnnulus
from photutils.detection import DAOStarFinder
from astropy.stats import mad_std
import matplotlib.pyplot as plt
from photutils.aperture import ApertureStats
from astropy.stats import SigmaClip
from scipy.ndimage import gaussian_filter


def aperturephotometry(file1,star_position,comp_position,platescale,dstarpc,star_inner_ann_radius,star_outer_ann_radius,comp_inner_ann_radius,comp_outer_ann_radius):
    
    file1 = file1
    
    star_position = star_position
    
    comp_position = comp_position
    
    platescale = platescale
    
    dstarpc = dstarpc
    
    star_inner_ann_radius = star_inner_ann_radius
    
    star_outer_ann_radius = star_outer_ann_radius
    
    comp_inner_ann_radius = comp_inner_ann_radius
    
    comp_outer_ann_radius = comp_outer_ann_radius
    

    '''Import Data'''
    
    
    file = pf.getdata(file1)
    file_hdr = pf.getheader(file1)
    
    
    '''Sky Annulus Photometry'''
    
    data = file
    
    
    '''star'''
    
    aperture = CircularAperture(star_position, r=star_inner_ann_radius)
    # 
    
    annulus_aperture = CircularAnnulus(star_position, r_in=star_inner_ann_radius, r_out=star_outer_ann_radius)
    # 
    aperstats = ApertureStats(data, annulus_aperture)
    
    bkg_mean = aperstats.mean
    
    plt.figure()
    plt.imshow(data, origin='lower',vmin=0,vmax=max(data))
    plt.show()
    
    annulus_aperture.plot(color='blue', lw=1.5, alpha=0.5)
    
    # print(bkg_mean)  
    
    phot_table = aperture_photometry(data, aperture)
    
    for col in phot_table.colnames:
    
        phot_table[col].info.format = '%.8g'  # for consistent table output
    
    # print(phot_table)
    
    aperture_area = aperture.area_overlap(data)
    
    # print(aperture_area)  
    
    total_bkg = bkg_mean * aperture_area
    
    # print(total_bkg)  
    
    phot_bkgsub = phot_table['aperture_sum'] - total_bkg
    
    phot_table['total_bkg'] = total_bkg
    
    phot_table['aperture_sum_bkgsub'] = phot_bkgsub
    
    for col in phot_table.colnames:
    
        phot_table[col].info.format = '%.8g'  # for consistent table output
    
    # print(phot_table)
    
    '''Median sky background'''
    
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    
    aper_stats = ApertureStats(data, aperture, sigma_clip=None)
    
    bkg_stats = ApertureStats(data, annulus_aperture, sigma_clip=sigclip)
    
    # print(bkg_stats.median)  
    
    total_bkg = bkg_stats.median * aper_stats.sum_aper_area.value
    
    # print('annulus sum (initial) = ', aper_stats.sum)
    # print('annulus size = ', aper_stats.sum_aper_area.value)
    # print('annulus std = ', aper_stats.std)
    
    # print('bkg sum = ', total_bkg)
    # print('bkg annulus size = ',bkg_stats.sum_aper_area.value)  
    # print('bkg std = ', bkg_stats.std)
    
    apersum_bkgsub = aper_stats.sum - total_bkg
    
    # print('final annulus sum = ', apersum_bkgsub)  
    
    '''companion'''
    
    comp_aperture = CircularAperture(comp_position, r=comp_inner_ann_radius)
    # 
    
    comp_annulus_aperture = CircularAnnulus(comp_position, r_in=comp_inner_ann_radius, r_out=comp_outer_ann_radius)
    # 
    comp_aperstats = ApertureStats(data, comp_annulus_aperture)
    
    comp_bkg_mean = comp_aperstats.mean
    
    # plt.figure()
    # plt.imshow(data, origin='lower',vmin=0,vmax=100)
    # plt.show()
    
    comp_annulus_aperture.plot(color='blue', lw=1.5, alpha=0.5)
    
    # print(comp_bkg_mean)  
    
    comp_phot_table = aperture_photometry(data, comp_aperture)
    
    for col in comp_phot_table.colnames:
    
        comp_phot_table[col].info.format = '%.8g'  # for consistent table output
    
    # print(comp_phot_table)
    
    comp_aperture_area = comp_aperture.area_overlap(data)
    
    # print(comp_aperture_area)  
    
    comp_total_bkg = comp_bkg_mean * comp_aperture_area
    
    # print(comp_total_bkg)  
    
    comp_phot_bkgsub = comp_phot_table['aperture_sum'] - comp_total_bkg
    
    comp_phot_table['total_bkg'] = comp_total_bkg
    
    comp_phot_table['aperture_sum_bkgsub'] = comp_phot_bkgsub
    
    for col in comp_phot_table.colnames:
    
        comp_phot_table[col].info.format = '%.8g'  # for consistent table output
    
    # print(comp_phot_table)
    
    '''Median sky background'''
    
    comp_aper_stats = ApertureStats(data, comp_aperture, sigma_clip=None)
    
    comp_bkg_stats = ApertureStats(data, comp_annulus_aperture, sigma_clip=sigclip)
    
    # print(comp_bkg_stats.median)  
    
    comp_total_bkg = comp_bkg_stats.median * comp_aper_stats.sum_aper_area.value
    
    # print('annulus sum (initial) = ', comp_aper_stats.sum)
    # print('annulus size = ', comp_aper_stats.sum_aper_area.value)
    # print('annulus std = ', comp_aper_stats.std)
    
    # print('bkg sum = ', comp_total_bkg)
    # print('bkg annulus size = ',comp_bkg_stats.sum_aper_area.value)  
    # print('bkg std = ', comp_bkg_stats.std)
    
    comp_apersum_bkgsub = comp_aper_stats.sum - comp_total_bkg
    
    # print('final annulus sum = ', comp_apersum_bkgsub)  
    
    starflux = apersum_bkgsub
    
    compflux = comp_apersum_bkgsub
    
    deltamag = -2.5*np.log10(np.abs(compflux)/starflux)
    
    dpix = np.sqrt((comp_position[0][0]-star_position[0][0])**2 + (comp_position[0][1]-star_position[0][1])**2)
    
    darcsec = dpix * platescale
    
    dau = darcsec*dstarpc
    
    print('deltamag = ', deltamag)
    
    print('distance (pix) = ', dpix)
    
    print('distance (") = ', darcsec)
    
    print('distance (AU) = ', dau)