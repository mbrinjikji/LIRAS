#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 15:32:57 2024

@author: marah
"""


import numpy as np
import math as m
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore',category=AstropyWarning)
import astropy.io.fits as pf
from lmfit import Model

def moffat(coords,a,xcen,ycen,sig,beta):
    term1 = (coords[0] - xcen)**2
    term2 = (coords[1] - ycen)**2
    term3 = sig**2
    return a*(1+(term1+term2)/term3)**(-1.*beta)


def moffat_model(file,comp_center,star_center):
    science = pf.getdata(file,header=False,ignore_missing_end=True)
    y,x = np.meshgrid(np.arange(len(science[0])),np.arange(len(science[0])))
    r = 10.
    sdist = np.where(((x-star_center[0])**2 + (y-star_center[1])**2) < r**2) #comp is short for companion
    #get x and y positions for what's in circle on companion
    comp_x = sdist[0]
    print(comp_x)
    comp_y = sdist[1]
    print(comp_y)
    counts = comp_x*0.
    #get the actual counts at each location
    for i in range(0,len(comp_x)):
        counts[i] = science[comp_x[i],comp_y[i]]
    all_comp = np.array((comp_x,comp_y)) #make array of positions
    guessa = counts.max()
    sallcounts = np.sum(counts)
    print('max counts =', guessa)
    print('all counts =', sallcounts)
    gmod = Model(moffat)
    # result = gmod.fit(counts,coords=all_comp,a=allcounts,xcen=comp_center[0],ycen=comp_center[1],sig=10,beta=2)
    sresult = gmod.fit(counts,coords=all_comp,a=sallcounts,xcen=star_center[0],ycen=star_center[1],sig=10,beta=2)
    srexcen = sresult.params['xcen'].value
    sreycen = sresult.params['ycen'].value
    srea = sresult.params['a'].value
    sresig = sresult.params['sig'].value
    srebeta = sresult.params['beta'].value
    
    
    cdist = np.where(((x-comp_center[0])**2 + (y-comp_center[1])**2) < r**2) #comp is short for companion
    #get x and y positions for what's in circle on companion
    comp_x = cdist[0]
    print(comp_x)
    comp_y = cdist[1]
    print(comp_y)
    counts = comp_x*0.
    #get the actual counts at each location
    for i in range(0,len(comp_x)):
        counts[i] = science[comp_x[i],comp_y[i]]
    all_comp = np.array((comp_x,comp_y)) #make array of positions
    guessa = counts.max()
    callcounts = np.sum(counts)
    print('max counts =', guessa)
    print('all counts =', callcounts)
    gmod = Model(moffat)
    cresult = gmod.fit(counts,coords=all_comp,a=callcounts,xcen=comp_center[0],ycen=comp_center[1],sig=10,beta=2)
    crexcen = cresult.params['xcen'].value
    creycen = cresult.params['ycen'].value
    crea = cresult.params['a'].value
    cresig = cresult.params['sig'].value
    crebeta = cresult.params['beta'].value

    
    
    return srexcen,sreycen,srea,sallcounts,sresig,srebeta,crexcen,creycen,crea,callcounts,cresig,crebeta


# file = '/Users/marah/Documents/gradschool/Research/Thesis/LBT/companions/epoch1/J04335252_rot2.fits'
# comp_center = [crexcen, creycen]
# star_center = [srexcen, sreycen]

# print(moffat_model(file,comp_center,star_center))



def seps(comp_center,star_center,platescale,dstarpc):
    dpix = np.sqrt((comp_center[0]-star_center[0])**2 + (comp_center[1]-star_center[1])**2) #calculation distance from candidate to star in pixel
    #print('Separation in Pixels = ', dpix)
    #print('parallax = ', parallax)

    #convert pixel distance into arcseconds
    darcsec = dpix * platescale #arcsec/pixel scale from LBTI corrections
    #print('Separation in arcsec = ', darcsec)
    #print('Distance to star (pc) = ', dstarpc)
    dau = darcsec*dstarpc #projected distance from candidate to star
    
    position_angle = np.arctan2((-1.0) * (comp_center[1]-star_center[1]), (comp_center[0]-star_center[0]))/ (2 * np.pi) * 360.

    
    return dpix, darcsec, dau, position_angle
    

# comp_center = [crexcen, creycen]
# star_center = [srexcen, sreycen]
# platescale = 0.010537
# dstarpc = 140

# print(seps(comp_center,star_center,platescale,dstarpc))
