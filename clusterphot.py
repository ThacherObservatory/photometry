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
from length import length


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
        
    if length(bias) == 1:
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

    if length(bias) == 1:
        bias = make_bias(path=path)

    if length(dark) == 1:
        dark = make_dark(path=path, bias=bias)
        
    flats,ct   = tp.get_files(dir=path, tag='SkyFlat.'+band)
    masterflat = tp.master_flat(flats, bias=bias, dark=dark, outdir='./', suffix=band)

    return masterflat
    
    
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

def stack_ims(path=None,band='gp',source='NGC188',bias=None,dark=None,flat=None):

    
    if not path:
        path = set_path()

    if length(bias) == 1:
        bias = make_bias(path=path)

    if length(dark) == 1:
        dark = make_dark(path=path)

    if length(flat) == 1:
        flat = make_flat(path=path,band=band,bias=bias,dark=dark)

        
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


    return final,refh

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

def find_stars(image, plot = False, fwhm = 20.0, threshold=3.):

    from astropy.stats import sigma_clipped_stats
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    from photutils import daofind
    sources = daofind(image - median, fwhm=fwhm, threshold=threshold*std)
    
   
   #for loop for vetting identified stars here; use numpy.delete for deleting elements from corresponding miniarrays
   # for i in len(sources):
        
    
    if plot == True:
       # from astropy.visualization import SqrtStretch
       # from astropy.visualization.mpl_normalize import ImageNormalize
        positions = (sources['xcentroid'], sources['ycentroid'])
        apertures = CircularAperture(positions, r=4.)
        #norm = ImageNormalize(stretch=SqrtStretch())
        #plt.imshow(image, cmap='Greys', origin='lower', norm=norm)
        qi.display_image(image)
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        
    return sources

#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#


def do_phot(image,position,radius = 5, r_in=15., r_out=20.):
    
    aperture = CircularAperture(position,r=radius)

    bkg_aperture = CircularAnnulus(position,r_in=r_in,r_out=r_out)

    # perform the photometry; the default method is 'exact'
    phot = aperture_photometry(image, aperture)
    bkg = aperture_photometry(image, bkg_aperture)

    # calculate the mean background level (per pixel) in the annuli
    bkg_mean = bkg['aperture_sum'] / bkg_aperture.area()
    bkg_sum = bkg_mean * aperture.area()
   
    #this may need editing; 'phot' in second line below may need brackets with 'flux_sum' inside
    phot['bkg_sum'] = bkg_sum
    phot['flux_sum'] = phot - bkg_sum
    
    return phot
    
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
def vet_sources(xcen,ycen,fwhm=20.0):
    '''
    Input the xcentroid and ycentriod values of sources in an image and the
    distance less than which sources will be rejected (in pixels)

    Returns vetted list of x and y centroid values.

    '''
    if len(xcen) != len(ycen):
        print 'X centroid and Y centroid values not equal in length!'
        return None,None

    s = np.argsort(xcen)
    xsort = xcen[s]
    ysort = ycen[s]
    n = len(xsort)
    i = 0
    while i < len(xsort)-1:
        x = xsort[i]
        y = ysort[i]
        xinds, = np.where(np.abs(xsort-x) < fwhm)
        distvec = np.sqrt( (xsort[xinds]-x)**2 + (ysort[xinds]-y)**2)
        dinds, = np.where(distvec < fwhm)
        if len(dinds) > 1:
            xsort = np.delete(xsort,xinds[dinds[1:]])
            ysort = np.delete(ysort,xinds[dinds[1:]])
        i += 1

    xvet = np.copy(xsort)
    yvet = np.copy(ysort)

    return xvet,yvet
