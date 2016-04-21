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

#this function needs work; path is not being set properly
def set_path(user == 'Yousef'):
    if user == 'Staniya':
        path = '/Users/staniya/Astronomy/MINERVA/2016Jan12/'
        #path= '/Users/staniya/Astronomy/MINERVA/2016Jan13'
    elif user == 'Yousef':
        path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan12/'
        return path
        #path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan13/'
    #elif path == 'Alden':
    #Alden, input your path here
    return

set_path('Yousef')
#-------------------------------------------------------------------------#
#-------------------------------------------------------------------------#

#get raw image
gimage, gheader = fits.getdata(path + 'n20160112.T4.NGC188.gp.0096.fits', 0, header = True)
rimage, rheader = fits.getdata(path + 'n20160112.T4.NGC188.rp.0107.fits', 0, header = True)
iimage, iheader = fits.getdata(path + 'n20160112.T4.NGC188.ip.0121.fits', 0, header = True)
zimage, zheader = fits.getdata(path + 'n20160112.T4.NGC188.zp.0140.fits', 0, header = True)

#get files
files = glob.glob(path+"*fits")

#Image Calibration -------------------------------------------------------#
#get bias,dark,flats files
biases,bct = tp.get_files(dir = path, tag='Bias')
darks,dct = tp.get_files(dir = path,tag='Dark')
gflats = tp.get_files(dir = path, tag = 'SkyFlat.gp')
iflats = tp.get_files(dir = path, tag = 'SkyFlat.ip')
rflats = tp.get_files(dir = path, tag = 'SkyFlat.rp')
zflats = tp.get_files(dir = path, tag = 'SkyFlat.zp')

#create master bias,flat
bias = tp.master_bias(biases,outdir='./')
dark = tp.master_dark(darks,bias=bias,outdir='./')

#master flats for the 4 photometric bands
gflat = tp.master_flat(gflats[0],bias=bias,dark=dark,outdir='./', suffix = 'gp')
rflat = tp.master_flat(rflats[0],bias=bias,dark=dark,outdir='./', suffix = 'rp')
iflat = tp.master_flat(iflats[0],bias=bias,dark=dark,outdir='./', suffix = 'ip')
zflat = tp.master_flat(zflats[0],bias=bias,dark=dark,outdir='./', suffix = 'zp')

#calibrated images for the 4 photometric bands
gcal = (gimage - dark - bias)/gflat
rcal = (rimage - dark - bias)/rflat
ical = (iimage - dark - bias)/iflat
zcal = (zimage - dark - bias)/zflat

"""
#display images
qi.display_image(gcal,siglo = 2, sighi = 2,fignum = 2)
qi.display_image(rcal,siglo = 2, sighi = 2,fignum = 2)
qi.display_image(ical,siglo = 2, sighi = 2,fignum = 2)
qi.display_image(zcal,siglo = 2, sighi = 2,fignum = 2)
"""

#Photometry on a Single Star --------------------------------------------#
#potential candidate: x ~= 1830, y = 75
#potential eclipsing binary: x = 66, y = 1350

#display image and plot an x over the star at x = 1830, y = 75
qi.display_image(gcal,siglo = 2, sighi = 2,fignum = 1)
plt.plot(1830,75,'x')

#aperture photometry
x0 = 1830.0
y0 = 75.0
radius,ynew,xnew,fwhm,aspect,snrmax,totflux,totap,chisq = \
    tp.optimal_aperture(x0,y0,gcal,skyrad=[15,20])

#Combining all images in a photometric band -----------------------------#

#get all images for each photometric band
gfiles,gct = tp.get_files(dir = path, tag = 'NGC188.gp')
rfiles,rct = tp.get_files(dir = path, tag = 'NGC188.rp')
ifiles,ict = tp.get_files(dir = path, tag = 'NGC188.ip')
zfiles,zct = tp.get_files(dir = path, tag = 'NGC188.zp')

"""
This was copied from the imaging script of Fall 2015 (A. Wood, K. O'Neill, M. Wilcox)
files = glob.glob('Mantis*[0-9]'+band+'_cal.fit*')
zsz = len(files)
reffile = files[zsz/2]
image0,header0 = readimage(reffile)
ysz,xsz = np.shape(image0)
refim = h.pyfits.open(reffile)
refh = h.pyfits.getheader(reffile)
stack = np.zeros((xsz,ysz,zsz))
for i in range(zsz):
    im = h.pyfits.open(files[i])
    newim = h.hcongrid(im[0].data,im[0].header,refh)
    stack[:,:,i] = newim
    
final = np.median(stack,axis=2)
""" 
#stack all images and produce final calibrated image for each photometric band
gzsz = len(gfiles)
greffile = gfiles[gzsz/2]
image0, header0 = qi.readimage(greffile)
gysz, gxsz = np.shape(image0)
grefim = h.pyfits.open(greffile)
grefh = h.pyfits.getheader(greffile)
gstack = np.zeros((gxsz,gysz,gzsz))
for i in range(gzsz):
    gim = h.pyfits.open(gfiles[i])
    gnewim = h.hcongrid((gim[0].data-dark-bias)/gflat, gim[0].header,grefh)
    gstack[:,:,i] = gnewim
    
gfinal = np.median(gstack, axis = 2)

rzsz = len(rfiles)
rreffile = rfiles[rzsz/2]
image0, header0 = qi.readimage(rreffile)
rysz, rxsz = np.shape(image0)
rrefim = h.pyfits.open(rreffile)
rrefh = h.pyfits.getheader(rreffile)
rstack = np.zeros((rxsz,rysz,rzsz))
for i in range(rzsz):
    rim = h.pyfits.open(rfiles[i])
    rnewim = h.hcongrid((rim[0].data-dark-bias)/rflat, rim[0].header,rrefh)
    rstack[:,:,i] = rnewim
    
rfinal = np.median(rstack, axis = 2)

izsz = len(ifiles)
ireffile = ifiles[izsz/2]
image0, header0 = qi.readimage(ireffile)
iysz, ixsz = np.shape(image0)
irefim = h.pyfits.open(ireffile)
irefh = h.pyfits.getheader(ireffile)
istack = np.zeros((ixsz,iysz,izsz))
for i in range(izsz):
    iim = h.pyfits.open(ifiles[i])
    inewim = h.hcongrid((iim[0].data-dark-bias)/iflat, iim[0].header,irefh)
    istack[:,:,i] = inewim
    
ifinal = np.median(istack, axis = 2)

zzsz = len(zfiles)
zreffile = zfiles[zzsz/2]
image0, header0 = qi.readimage(zreffile)
zysz, zxsz = np.shape(image0)
zrefim = h.pyfits.open(zreffile)
zrefh = h.pyfits.getheader(zreffile)
zstack = np.zeros((zxsz,zysz,zzsz))
for i in range(zzsz):
    zim = h.pyfits.open(zfiles[i])
    znewim = h.hcongrid((zim[0].data-dark-bias)/zflat, zim[0].header,zrefh)
    zstack[:,:,i] = znewim
    
zfinal = np.median(zstack, axis = 2)

#framework; this was just copied from ipython notebook; later to be custom tailored
from photutils.datasets import load_star_image
hdu = load_star_image()    
star_data = hdu.data[0:400, 0:400]

from astropy.stats import sigma_clipped_stats
mean, median, std = sigma_clipped_stats(star_data, sigma=3.0)

from photutils import daofind
sources = daofind(star_data - median, fwhm=3.0, threshold=5.*std) 





