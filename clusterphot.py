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
#actually type in your path in the if statements
path = 'Yousef'

#if path == 'Shin':
    #Shin, input your path here
if path == 'Yousef':
    path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan12/'
    #path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/2016Jan13/'
    import swiftimage as qi
#if path == 'Alden':
    #Alden, input your path here

#########################################################################

#test invidiual image
qi.readimage('n20160112.T4.NGC188.gp.0096.fits')

files = glob.glob(path+"*fits")

#image calibration using biases
biases,bct = tp.get_files(dir = path, tag='Bias')
darks,dct = tp.get_files(dir = path,tag='Dark')
#flats,fct = tp.get_files(dir = path, tag='SkyFlat')
gflats = tp.get_files(dir = path, tag = 'SkyFlat.gp')
iflats = tp.get_files(dir = path, tag = 'SkyFlat.ip')
rflats = tp.get_files(dir = path, tag = 'SkyFlat.rp')
zflats = tp.get_files(dir = path, tag = 'SkyFlat.zp')

bias = tp.master_bias(biases,outdir='./Photometry/')

dark = tp.master_dark(darks,bias=bias,outdir='./Photometry/')

#master flats for the 4 photometric bands
gflat = tp.master_flat(gflats,bias=bias,dark=dark,outdir='./Photometry/')
iflat = tp.master_flat(iflats,bias=bias,dark=dark,outdir='./Photometry/')
rflat = tp.master_flat(rflats,bias=bias,dark=dark,outdir='./Photometry/')
zflat = tp.master_flat(zflats,bias=bias,dark=dark,outdir='./Photometry/')

#calibrated image equations for the 4 photometric bands
gcal = (image - dark - bias)/gflat
ical = (image - dark - bias)/iflat
rcal = (image - dark - bias)/rflat
zcal = (image - dark - bias)/zflat




