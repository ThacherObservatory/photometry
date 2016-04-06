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
#import quickimage as qi

#to change the path, just change the assignment below to your name (but don't forget to 
#actually type in your path in the if statements, as well as any potential import statements)
path = 'Yousef'

#if path == 'Shin':
    path= '/Users/shinnosuke/Documents/Astronomy/__MACOSX/2016Jan12/'
    #path= 'Users/shinnosuke/Documents/Astronomy/__MACOSX/2016Jan13/'
    import swiftimage as qi
if path == 'Yousef':
    path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan12/'
    #path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan13/'
    import swiftimage as qi
#if path == 'Alden':
    #Alden, input your path here

#########################################################################

#test invidiual image
image, header = fits.getdata(path + 'n20160112.T4.NGC188.gp.0096.fits', 0, header = True)

files = glob.glob(path+"*fits")

#image calibration using biases
biases,bct = tp.get_files(dir = path, tag='Bias')
darks,dct = tp.get_files(dir = path,tag='Dark')
#flats,fct = tp.get_files(dir = path, tag='SkyFlat')
gflats = tp.get_files(dir = path, tag = 'SkyFlat.gp')
iflats = tp.get_files(dir = path, tag = 'SkyFlat.ip')
rflats = tp.get_files(dir = path, tag = 'SkyFlat.rp')
zflats = tp.get_files(dir = path, tag = 'SkyFlat.zp')

bias = tp.master_bias(biases,outdir='./')

dark = tp.master_dark(darks,bias=bias,outdir='./')

#master flats for the 4 photometric bands
gflat = tp.master_flat(gflats[0],bias=bias,dark=dark,outdir='./')
iflat = tp.master_flat(iflats[0],bias=bias,dark=dark,outdir='./')
rflat = tp.master_flat(rflats[0],bias=bias,dark=dark,outdir='./')
zflat = tp.master_flat(zflats[0],bias=bias,dark=dark,outdir='./')

#calibrated image equations for the 4 photometric bands
#note: image will be defined later (when a specific image is chosen)
gcal = (image - dark - bias)/gflat
ical = (image - dark - bias)/iflat
rcal = (image - dark - bias)/rflat
zcal = (image - dark - bias)/zflat

qi.readimage(image)



