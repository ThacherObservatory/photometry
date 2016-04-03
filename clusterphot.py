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


files = glob.glob(path+"*fits")

#image calibration (bias, dark,flat)
darks,dct = tp.get_files(tag='Dark')
biases,bct = tp.get_files(tag='Bias')
flats,fct = tp.get_files(tag='SkyFlat')

bias = tp.master_bias(biases,outdir='./Photometry/',readnoise=False)

dark = tp.master_dark(darks,bias=bias,outdir='./Photometry/')

flat = tp.master_flat(flats,bias=bias,dark=dark,outdir='./Photometry/')
