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
import hcongrid


#to change the path, just change the assignment below to your name (but don't forget to 
#actually type in your path in the if statements, as well as any potential import statements)
path = 'Yousef'

#def get_path
if path == 'Staniya':
    path= '/Users/staniya/Astronomy/MINERVA/2016Jan12/'
    #path= '/Users/staniya/Astronomy/MINERVA/2016Jan13'
if path == 'Yousef':
    path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan12/'
    #path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan13/'
#if path == 'Alden':
    #Alden, input your path here

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

#display images
qi.display_image(gcal,siglo = 2, sighi = 2,fignum = 2)
qi.display_image(rcal,siglo = 2, sighi = 2,fignum = 2)
qi.display_image(ical,siglo = 2, sighi = 2,fignum = 2)
qi.display_image(zcal,siglo = 2, sighi = 2,fignum = 2)

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
gimages = tp.get_files(dir = path, tag = 'NGC188.gp')
rimages = tp.get_files(dir = path, tag = 'NGC188.rp')
iimages = tp.get_files(dir = path, tag = 'NGC188.ip')
zimages = tp.get_files(dir = path, tag = 'NGC188.zp')

#Note: not necessary; combine all images into one master image for each photometric band
#mastergcal = hcongrid.hcongrid(gimages, )
#masterical = hcongrid.hcongrid(gimages, )
#masterrcal = hcongrid.hcongrid(gimages, )
#masterzcal = hcongrid.hcongrid(gimages, )

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


#framework; this was just copied from ipython notebook; later to be custom tailored
from photutils.datasets import load_star_image
hdu = load_star_image()    
star_data = hdu.data[0:400, 0:400]

from photutils import daofind
sources = daofind(star_data - median, fwhm=3.0, threshold=5.*std) 






