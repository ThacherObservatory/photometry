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
import quickimage as qi

#Yousef notes: for me at least there is a problem with import statements,
#particularly with thacherphot not being able to import the module fitgaussian

#if path == 'Shin':
    #Shin's path
if path == 'Yousef':
    path = '/Users/Yousef/Desktop/Astronomy/MINERVA Star Clusters/'
#if path == 'Alden':
    #Alden's path


files = glob.glob(path+"*fits")
#test image in files
qi.readimage(files[10],plot=True)

image,header = qi.readimage(files[50],plot=True,siglo=2,sighi=2)
w = wcs.WCS(header)


#image calibration (bias, dark,flat)
darks,dct = tp.get_files(tag='Dark')
biases,bct = tp.get_files(tag='Bias')
flats,fct = tp.get_files(tag='SkyFlat')

bias = tp.master_bias(biases,outdir='./Photometry/',readnoise=False)

dark = tp.master_dark(darks,bias=bias,outdir='./Photometry/')

flat = tp.master_flat(flats,bias=bias,dark=dark,outdir='./Photometry/')
