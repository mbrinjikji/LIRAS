import os
import numpy as np
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore',category=AstropyWarning)
import astropy.io.fits as pf
import pyklip.instruments.Instrument as Instrument
import pyklip.parallelized as parallelized
import pyklip.klip as pyklip
import pyklip.rdi as rdi
import time

def findfiles(path):
    for (path,dirs,files) in os.walk(path):
        return files

def remove(string):
    return string.replace(" ","")


#get all the files that have been 
def getobjectlist(datadir):
    print(datadir)
    objectlist = []
    files = findfiles(datadir)
    #print(files)
    files.sort()
    filesc = []
    for i in files:
        if '.fits' in i:
            filesc.append(i)
    #print(filesc)
    for f in filesc:
        #print(f)
        science = pf.open(datadir+f,mode='update',ignore_missing_end=True)
        OBJECT = science[0].header['OBJNAME']
        ITIME = science[0].header['ITIME']
        # if ITIME == 60 or ITIME == 30 or ITIME == 10 or ITIME == 3 or ITIME == 0.1:
        if ITIME > 0:
            objectlist.append(OBJECT)
        science.close()
    return objectlist,files

def preppsflib(star_center,mode,prefix,ann,subs,numbs,make_corr,targetname,outpath,corr_path,datapsf,objectlist,files): 
    filesc = []
    for i in files:
        if '.fits' in i:
            filesc.append(i)
    

 	#get all the files used for the PSF library and the target star into a list
    objectlistu = np.unique(objectlist)
    print(len(objectlistu))

    starlist = []
    filelist = []
    targetstar = targetname
    print('targetname = ', targetname)
    for n in filesc:
        #padname = n.replace('.fits','_padded.fits')
        science = pf.open(datapsf+ '/' + n,mode='update',ignore_missing_end=True)
        OBJECT = science[0].header['OBJNAME']
        ITIME = science[0].header['ITIME']
        # if ITIME == 60 or ITIME == 30 or ITIME == 10 or ITIME == 3 or ITIME == 0.1:
        if ITIME > 0:
            if OBJECT in objectlistu:
                filelist.append(datapsf+n)
            if targetstar == OBJECT:
                starlist.append(datapsf+n)
                #starlist.append(datatarget+padname)
        science.close()
    #print('filelist = ', filelist)
    #print(len(filelist))
    print('starlist = ', starlist)
    print(len(starlist))


    #get the center of all the stars from one file header    
    sciencehead = pf.open(filelist[0],ignore_missing_end=True)
    klipxcen = star_center[1]
    klipycen = star_center[0]
    print('padded center = ',klipxcen,klipycen)
    time.sleep(5)
    aligned_center = [klipxcen,klipycen]
    print('aligned_center = ',aligned_center)
    sciencehead.close()


    #make the dataset for the correlation matrix and psf library
    dataset = Instrument.GenericData(filelist,guess_star=[klipxcen,klipycen],highpass=False,find_star=False)
    # for i in range(0,len(dataset.input)):
    #     dataset[i] = pyklip.align_and_scale(dataset.input[i],[klipycen,klipxcen],old_center=dataset.centers[i])
        # print['aligned dataset centers = ', dataset.centers[i]]
   # print('dataset filenames =', dataset.filenames)

    psflib_imgs = dataset.input
    psflib_filenames = dataset.filenames
    # print('psflib_filenames = ', psflib_filenames)
    print('dataset center = ', dataset.centers[0],dataset.centers[1],dataset.centers[2])
    print(len(dataset.centers))
    for c in range(0,len(dataset.centers)):
        dataset.centers[c] = [klipycen,klipxcen]
    print('dataset center = ', dataset.centers[0],dataset.centers[1],dataset.centers[2])

    

    #make the correlation matrix
    print('Will you make the correlation matrix?')
    if make_corr == True:
        print('You have decided to make the correlation matrix from scratch.')
        psflib = rdi.PSFLibrary(psflib_imgs,aligned_center,psflib_filenames,compute_correlation=True)
        psflib.save_correlation(corr_path,overwrite=True)
        corr_matrix_hdulist = pf.open(corr_path)
        corr_matrix = corr_matrix_hdulist[0].data
        psflib = rdi.PSFLibrary(psflib_imgs,aligned_center,psflib_filenames,correlation_matrix=corr_matrix)

    #if you already have the correlation matrix, just open it
    if make_corr == False:
        print('You have decided to read in an existing correlation matrix.')
            
    stardataset = Instrument.GenericData(starlist,guess_star=[klipxcen,klipycen],highpass=False,find_star=False)
    #print(stardataset)
    for c in range(0,len(stardataset.centers)):
        stardataset.centers[c] = [klipycen,klipxcen]


    psflib.prepare_library(stardataset) #this needs to be just the star you are running the reduction on, take out this star from the library


    objectname = remove(targetstar)
    print(objectname)

    parallelized.klip_dataset(stardataset,mode=mode,outputdir=outpath,fileprefix=objectname+prefix,annuli=ann,subsections=subs,movement=0,numbasis=numbs,numthreads=8,minrot=0,calibrate_flux=False,aligned_center=[klipxcen,klipycen],annuli_spacing="log",maxnumbasis=100,corr_smooth=1,spectrum=None,psf_library=psflib,highpass=False,lite=False,save_aligned=False,restored_aligned=None,dtype=None,algo="klip",time_collapse="median",wv_collapse="median")


def runeachobject(outpath,datapsf,corr_path,star_center,mode,prefix,ann,subs,numbs):
    #path where you want the output images from klip to go
    outpath = outpath
    #path where the padded fits files are stored
    datapsf = datapsf
    
    corr_path = corr_path

    star_center = star_center

    mode = mode

    prefix = prefix

    ann = ann

    subs = subs

    numbs = numbs

    objectlist,files = getobjectlist(datapsf)
    objectlistu = np.unique(objectlist)
    print(objectlistu)
    numobjects = len(objectlistu)
    print(numobjects)
    count=0
    desiredlist = ['Name of Object(s) to run KLIP on']
    for i in objectlistu:
        print(i)
        if i in desiredlist:
            print('starting')
            print(i)
            time.sleep(5)
            #print(count+ '/' +numobjects)
            preppsflib(star_center=star_center,mode=mode,prefix=prefix,ann=ann,subs=subs,numbs=numbs,filemake_corr=False,targetname=i,outpath=outpath,corr_path=corr_path,datapsf=datapsf,objectlist=objectlist,files=files)
            count += 1

# if __name__ == '__main__':
#     runeachobject()