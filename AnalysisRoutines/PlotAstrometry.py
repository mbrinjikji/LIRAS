from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric
from astropy.time import Time
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

plt.rcParams['xtick.labelsize']=15
plt.rcParams['ytick.labelsize']=15


def plot_astrometry(mjd_ref,dstar,ra_deg,dec_deg,starpmra,starpmde,sep_epoch1,pa_epoch1,sep_epoch2,pa_epoch2):

    mjd_ref = mjd_ref #obs date in mjd
    
    dstar = dstar ##distance to star in mas

    ra_deg = ra_deg ###RA of star in degrees

    dec_deg = dec_deg ###DEC of star in degrees

    starpmra = starpmra ###RA proper motion of star [mas/yr]

    starpmde = starpmde ###DEC proper motion of star [mas/yr]

    sep_epoch1 = sep_epoch1 #separation of candidate in first observing epoch in arcsec
    pa_epoch1 = pa_epoch1 #position angle of candidate in first observing epoch in degrees
    
    sep_epoch2 = sep_epoch2 #separation of candidate in second observing epoch in arcsec
    pa_epoch2 = pa_epoch1 #position angle of candidate in second observing epoch in degrees
    
    
    mjd_list = np.arange(5000.0, dtype = np.float64) -2500.0 + mjd_ref 
    t = Time(mjd_list, format = 'mjd')
    
    earth_position = get_body_barycentric('earth', t) #Get earth position in AU for each t
    
    earth_x = np.asarray([earth_position[i].x.value for i in range(0, len(earth_position))])
    earth_y = np.asarray([earth_position[i].y.value for i in range(0, len(earth_position))])
    earth_z = np.asarray([earth_position[i].z.value for i in range(0, len(earth_position))])
    
    mjd_ind = np.argmin(np.abs(mjd_list - mjd_ref))
    
    plx = dstar / 1000.0 
    
    ra = np.radians(ra_deg) #ra and dec coordinates in degrees
    dec = np.radians(dec_deg)
    
    
    # Now calculate the motion of the star in RA/Dec due to parallax only
    plx_ra = plx * (earth_x * np.sin(ra) - earth_y * np.cos(ra))
    plx_de = plx * ((earth_x * np.cos(ra) * np.sin(dec)) + (earth_y * np.sin(ra) * np.sin(dec)) - (earth_z * np.cos(dec)))
    
    star_pmra = starpmra/1000.0 #The numerator is in mas/year
    star_pmde = starpmde/1000.0
    
    pm_ra = ((mjd_list-mjd_ref)/365.25 * star_pmra)
    pm_de = ((mjd_list-mjd_ref)/365.25 * star_pmde)
    
    tot_ra = pm_ra + plx_ra
    tot_de = pm_de + plx_de
    
    dr = np.pi/180.0
    rd = 180.0/np.pi
        
    delta_ra = (pm_ra + (plx_ra - plx_ra[mjd_ind])) * (-1.0)
    delta_de = (pm_de + (plx_de - plx_de[mjd_ind])) * (-1.0)
    
    offset_ra = (sep_epoch1 * np.sin(pa_epoch1 * dr)) + delta_ra 
    offset_de = (sep_epoch1 * np.cos(pa_epoch1 * dr)) + delta_de
    
    offset_rho = np.sqrt((offset_ra**2.0) + (offset_de**2.0))
    offset_theta = (((np.arctan2(offset_de, -offset_ra) * rd) + 270.0) % 360.0)
    
    goi_ra_offset_epoch1 = sep_epoch1 * np.sin(pa_epoch1 * dr)
    goi_de_offset_epoch1 = sep_epoch1 * np.cos(pa_epoch1 * dr)
    
    goi_ra_offset_epoch2 = sep_epoch2 * np.sin(pa_epoch2 * dr)
    goi_de_offset_epoch2 = sep_epoch2 * np.cos(pa_epoch2 * dr)
    
    
    fig,ax = plt.subplots(1,1,figsize=(8,6),layout='tight')
    ax.plot(offset_ra,offset_de,label='Background Track')
    ax.plot(goi_ra_offset_epoch1,goi_de_offset_epoch1,marker='^',label='Epoch 1',linestyle='None')
    ax.plot(goi_ra_offset_epoch2,goi_de_offset_epoch2,marker='^',label='Epoch 2',linestyle='None')
    ax.set_xlabel(r'$\Delta$RA [arcsec]',fontsize=18)
    ax.set_ylabel(r'$\Delta$Dec [arcsec]',fontsize=18)
    ax.legend(loc='best',fontsize=16)
    plt.show()

