import thacherphot as tp
import numpy as np
import matplotlib.pyplot as plt
import glob,sys,time
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.visualization import scale_image
from photutils import CircularAperture, aperture_photometry, CircularAnnulus
import quick_image as qi
import sys
import hcongrid as h


#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

def set_path(user='Yousef',date='2016Jan12'):
    """
    Returns data path of user identified in the input of the function.

    example
    -------
    path = set_path(user='Yousef')

    """
    if user == 'Staniya':
        path = '/Users/staniya/Astronomy/MINERVA/'+date+'/'

    elif user == 'Yousef':
        path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/'+date+'/'

    elif path == 'Alden':
        path = '/Users/Alden/python/'+date+'/'

    elif path == 'Swift':
        path = '/Users/jonswift/Dropbox (Thacher)/Observatory/Data/MINERVA/'+date+'/'
        
    return path

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
def make_bias(path=None):

    """
    Make master bias from raw data files
    """
    if not path:
        path = set_path()
    
    biases,bct = tp.get_files(dir=path, tag='Bias')
    masterbias = tp.master_bias(biases,outdir='./')

    return masterbias


#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
def make_dark(path=None, bias=None):

    """
    Make master dark from raw data files
    """
    if not path:
        path = set_path()
        
    if not bias:
        bias = make_bias(path=path)

    darks,dct = tp.get_files(dir = path,tag='Dark')
    masterdark = tp.master_dark(darks,bias=bias,outdir='./')

    return masterdark


#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
def make_flat(path=None,band='gp',bias=None,dark=None):

    """
    Make master flat from raw data files in specified band
    """
    if not path:
        path = set_path()

    if not bias:
        bias = make_bias(path=path)

    if not dark:
        dark = make_dark(path=path, bias=bias)
        
    flats,ct   = tp.get_files(dir=path, tag='SkyFlat.'+band)
    masterflat = tp.master_flat(flats, bias=bias, dark=dark, outdir='./', suffix=band)

    return masterflat
    
    
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

def stack_ims(band='gp',source='NGC188',dark=None, bias=None, flat=None, path=None):

    
    if not path:
        path = set_path()

    if not bias:
        bias = make_bias(path=path)

    if not dark:
        dark = make_dark(path=path)

    if not flat:
        dark = make_flat(path=path,bias=bias,dark=dark,band=band)

        
    files,sz = tp.get_files(dir=path, tag=source+'.'+band)

    # Designate reference file 
    reffile = files[sz/2]
    image0, header0 = qi.readimage(reffile)
    ysz, xsz = np.shape(image0)
    refim = h.pyfits.open(reffile)
    refh = h.pyfits.getheader(reffile)
    stack = np.zeros((xsz,ysz,sz))
    
    for i in range(sz):
        im = h.pyfits.open(files[i])
        newim = h.hcongrid((im[0].data-dark-bias)/flat, im[0].header,refh)
        stack[:,:,i] = newim
    
    final = np.median(stack, axis=2)


    return final

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

def find_stars(image):

    star_data = image.data[0:400, 0:400]
    from astropy.stats import sigma_clipped_stats
    mean, median, std = sigma_clipped_stats(star_data, sigma=3.0)
    from photutils import daofind
    sources = daofind(star_data - median, fwhm=3.0, threshold=5.*std)
    
    return sources

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
'''
def plot_stars()



'''

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
'''
def phot_all():
    
    
    
'''