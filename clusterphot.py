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
    
   # stars already found accurately, vet_sources will be implemented when working properly
   # vet_sources(10.0,10.0)
        
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
   
    #look at ipython notebook; his may need editing; 'phot' in second line below may need brackets with 'flux_sum' inside
    #phot['bkg_sum'] = bkg_sum
    #phot['flux_sum'] = phot['flux'] - bkg_sum
    
    #these two lines from ipython notebook
    flux_bkgsub = phot['aperture_sum'] - bkg_sum
    phot['aperture_sum_bkgsub'] = flux_bkgsub
    
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
    
    
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#
def phot_all(image, xcen, ycen):
    
    flux = []
    for i in range(len(xcen)):
        phot = do_phot(image, [(xcen[i],ycen[i])])
        flux = np.append(flux,phot['aperture_sum_bkgsub'])
        
    return flux
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

def cal_image(plot = True,path=None,band='gp',source='BD710031',bandindex = 3,bias=None,dark=None,flat=None):
    
    if not path:
        path = set_path()

    if length(bias) == 1:
        bias = make_bias(path=path)

    if length(dark) == 1:
        dark = make_dark(path=path)

    if length(flat) == 1:
        flat = make_flat(path=path,band=band,bias=bias,dark=dark)
        
        
    files,sz = tp.get_files(dir=path, tag=source+'.'+band)
    
    reffile = files[bandindex]
    image0, header0 = qi.readimage(reffile)
    refh = h.pyfits.getheader(reffile)
   
    im = h.pyfits.open(files[bandindex])
    newim = h.hcongrid((im[0].data-dark-bias)/flat, im[0].header,refh)
    
    if plot:
        qi.display_image(newim)
        
    return newim,header0
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#   

def headerplot(fluxdata,ykeyword,badimindices = [],band = 'gp',path = None, source = 'BD710031', bias = None, dark=None, flat=None):
    
    
    if not path:
        path = set_path()

    if length(bias) == 1:
        bias = make_bias(path=path)

    if length(dark) == 1:
        dark = make_dark(path=path)

    if length(flat) == 1:
        flat = make_flat(path=path,band=band,bias=bias,dark=dark)
    
    files,osz = tp.get_files(dir=path, tag=source+'.'+band)
    #remove faulty images
    sz = osz
    if len(badimindices) != 0:
        for g in range(osz):
            if g in badimindices:
                del files[g]
                sz=sz-1
    
    keyworddata = []
    headers = []
    
    #extract headers from each image
    for i in range(sz):
        files[i],header = cal_image(plot=False,path=path,band=band,source=source,bandindex=i,bias=bias,dark=dark,flat=flat)
        headers.append(header)
      
    #extract keyword data from each header 
    for j in range(sz):
        header = headers[j]
        keyworddata.append(header[ykeyword])
        
    #determine peak flux value, index and its header for annotation 
    peakflux = max(fluxdata)
    peakfluxindex = fluxdata.argmax()
    keyworddataofpeakflux = keyworddata[peakfluxindex]
        
    #normalize flux data
    nfluxdata = fluxdata
    for q in range(len(fluxdata)):
        nfluxdata[q] = nfluxdata[q]/peakflux
    
    #display plot; set y axis
    plt.scatter(keyworddata,nfluxdata)
    axes = plt.gca()
    axes.set_ylim([0,1])
    
    #labels
    plt.title(ykeyword.capitalize()+' v. Flux (Band: ' + band + ')')
    plt.xlabel(ykeyword.capitalize())
    plt.ylabel('Flux')
    xcoord = keyworddataofpeakflux
    ycoord = 1
    plt.annotate('Peak Flux = '+str(peakflux), xy=(xcoord,ycoord), xytext=(xcoord ,ycoord - .05))
    
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#   

def flux_all(path=None,band='gp',source='BD710031',badimindices = [],bias=None,dark=None,flat=None):

    if not path:
        path = set_path()

    if length(bias) == 1:
        bias = make_bias(path=path)

    if length(dark) == 1:
        dark = make_dark(path=path)

    if length(flat) == 1:
        flat = make_flat(path=path,band=band,bias=bias,dark=dark)
        
    files,osz = tp.get_files(dir=path, tag=source+'.'+band)
    sz = osz
    #remove faulty images
    if len(badimindices) != 0:
        for g in range(osz):
            if g in badimindices:
                del files[g]
                sz=sz-1
            
    flux = []
    
    for i in range(sz):
        image,header = cal_image(plot=False,path=path,band=band,source=source,bandindex=i,bias=bias,dark=dark,flat=flat)
        w = wcs.WCS(header)
        x,y = w.all_world2pix(float(header['TARGRA']), float(header['TARGDEC']), 1)
        position = (float(x),float(y))
        phot = do_phot(image, position)
        flux = np.append(flux,phot['aperture_sum_bkgsub'])
        
    return flux
    
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#  

def display_raws(path=None,band='gp',source='BD710031'):
    
    if not path:
        path = set_path()
    
    files,sz = tp.get_files(dir=path, tag=source+'.'+band)
    
    for i in range(len(files)):
        image = fits.getdata(files[i],0)
        qi.display_image(image,fignum=i)
    
    
    
    