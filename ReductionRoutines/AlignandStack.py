#!/usr/local/bin/python
#------------------------------------------------------------------------------------------------
import numpy as np
import glob
import os
import astropy.io.fits as pf
import sep
from scipy.ndimage.interpolation import shift
#------------------------------------------------------------------------------------------------
'''
This part aligns all the images and then combines them into one. 
'''

#------------------------------------------------------------------------------------------------
def alignandstack(pathforalign,pathforstack,finalimage,thresholdvalue):
    
    def getCoo(thresholdvalue,scidat):
        '''this step finds the star coordinates. it can be finicky so i try and get the thresh/minarea as close to the value of the star as possible'''
        data = data = scidat.byteswap(True).view(scidat.dtype.newbyteorder())
        # data[data>1500] = 1.
        bkg = sep.Background(data)
        print(bkg.globalrms)
        thresh = thresholdvalue * bkg.globalrms
        # thresh = 700.
        objects = sep.extract(data, thresh, minarea=50)

        print(objects['x'], objects['y'], objects['npix'])

        return ( objects['x'], objects['y'] )
    
    def alignIm(pathforalign, pathforstack, fits):
        
        pathforalign = pathforalign
        
        pathforstack = pathforstack
        
        scidat = pf.getdata(fits)
        scihdr = pf.getheader(fits)
        Im = np.copy(scidat)
        xcoo, ycoo = getCoo(thresholdvalue,scidat)

        if len(xcoo) > 0:
            shIm = []
            if len(xcoo) == 1: 
                Im = np.pad( Im, (512,512), 'constant', constant_values=(0,0) )
                ysh, xsh = 512-xcoo[0], 512-ycoo[0]
                tmpIm = shift( Im, (xsh, ysh), order=1)
                shIm.append( tmpIm )
                inimage = tmpIm 

            if len(xcoo) > 1: 
                # maskS2 = Im - makeCircularMask(Im.shape,30,xcoo[1],ycoo[1])*Im
                # maskS1 = Im - makeCircularMask(Im.shape,30,xcoo[0],ycoo[0])*Im

                # maskS2 = np.pad( maskS2, (512,512), 'constant', constant_values=(0,0) )
                # maskS1 = np.pad( maskS1, (512,512), 'constant', constant_values=(0,0) )

                # ysh, xsh = 512-xcoo[0], 512-ycoo[0]
                # shIm.append( shift( maskS2, (xsh, ysh), order=1) )

                # ysh, xsh = 512-xcoo[1], 512-ycoo[1]
                # shIm.append( shift( maskS1, (xsh, ysh), order=1) )

                # inimage = np.mean(shIm,axis=0)
                
                Im = np.pad( Im, (512,512), 'constant', constant_values=(0,0) )
                ysh, xsh = 512-xcoo[0], 512-ycoo[0]
                tmpIm = shift( Im, (xsh, ysh), order=1)
                shIm.append( tmpIm )
                inimage = tmpIm 

            inmax = inimage.max()
            inmin = inimage.min()
            newimage = (inimage - inmin)/(inmax - inmin)

            outimage = newimage

            # Reormalize the array 
            finIm = outimage*(inmax - inmin) + inmin 
            pf.writeto(pathforstack+fits,finIm,scihdr,overwrite=True,output_verify='ignore')
            return ( finIm )

        return ( [] )
    
    pathforalign = pathforalign
    
    pathforstack = pathforstack
    
    finalimage = finalimage
    
    thresholdvalue = thresholdvalue
    
    os.chdir(pathforalign)
    fitslst = glob.glob('*skysub.fits')
    fitslst.sort()
    print(fitslst)

    '''this line runs the align step. keep for aligning, comment when combining all aligned files into one'''
    if not os.path.exists('align'): os.mkdir('align')
    results = [alignIm(pathforalign,pathforstack,fits) for fits in fitslst]
    
    os.chdir(pathforstack)
    fitslist = glob.glob('*skysub.fits')
    fitslist.sort()
    print(fitslist)


    '''the next six lines are the combine step. uncomment and run after commenting the line above'''
    science = [pf.getdata(fits) for fits in fitslist]
    scihdr = [pf.getheader(fits) for fits in fitslist]
    scihdr = scihdr[0]
    science = np.array(science)
    sci = np.median(science, axis=0)

    pf.writeto(finalimage, sci, scihdr, overwrite=True)


# #------------------------------------------------------------------------------------------------
# def alignIm(pathforstack, fits):
    
#     pathforstack = pathforstack
    
#     scidat = pf.getdata(fits)
#     scihdr = pf.getheader(fits)
#     Im = np.copy(scidat)
#     xcoo, ycoo = getCoo(scidat)

#     if len(xcoo) > 0:
#         shIm = []
#         if len(xcoo) == 1: 
#             Im = np.pad( Im, (512,512), 'constant', constant_values=(0,0) )
#             ysh, xsh = 512-xcoo[0], 512-ycoo[0]
#             tmpIm = shift( Im, (xsh, ysh), order=1)
#             shIm.append( tmpIm )
#             inimage = tmpIm 

#         if len(xcoo) > 1: 
#             # maskS2 = Im - makeCircularMask(Im.shape,30,xcoo[1],ycoo[1])*Im
#             # maskS1 = Im - makeCircularMask(Im.shape,30,xcoo[0],ycoo[0])*Im

#             # maskS2 = np.pad( maskS2, (512,512), 'constant', constant_values=(0,0) )
#             # maskS1 = np.pad( maskS1, (512,512), 'constant', constant_values=(0,0) )

#             # ysh, xsh = 512-xcoo[0], 512-ycoo[0]
#             # shIm.append( shift( maskS2, (xsh, ysh), order=1) )

#             # ysh, xsh = 512-xcoo[1], 512-ycoo[1]
#             # shIm.append( shift( maskS1, (xsh, ysh), order=1) )

#             # inimage = np.mean(shIm,axis=0)
            
#             Im = np.pad( Im, (512,512), 'constant', constant_values=(0,0) )
#             ysh, xsh = 512-xcoo[0], 512-ycoo[0]
#             tmpIm = shift( Im, (xsh, ysh), order=1)
#             shIm.append( tmpIm )
#             inimage = tmpIm 

#         inmax = inimage.max()
#         inmin = inimage.min()
#         newimage = (inimage - inmin)/(inmax - inmin)

#         outimage = newimage

#         # Reormalize the array 
#         finIm = outimage*(inmax - inmin) + inmin 
#         pf.writeto(pathforstack+fits,finIm,scihdr,overwrite=True,output_verify='ignore')
#         return ( finIm )

#     return ( [] )

# #------------------------------------------------------------------------------------------------
# def getCoo(thresholdvalue,scidat):
#     '''this step finds the star coordinates. it can be finicky so i try and get the thresh/minarea as close to the value of the star as possible'''
#     data = scidat.byteswap(True).newbyteorder()
#     # data[data>1500] = 1.
#     bkg = sep.Background(data)
#     print(bkg.globalrms)
#     thresh = thresholdvalue * bkg.globalrms
#     # thresh = 700.
#     objects = sep.extract(data, thresh, minarea=50)

#     print(objects['x'], objects['y'], objects['npix'])

#     return ( objects['x'], objects['y'] )

# #------------------------------------------------------------------------------------------------
