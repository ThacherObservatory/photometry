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

if path == 'Shin'
    //Shin's path
if path == 'Yousef'
    //Yousef's path
if path == 'Alden'
    //Alden's path


files = glob.glob(path+"*fits")

image,header = qi.readimage(files[50],plot=True,siglo=2,sighi=2)
w = wcs.WCS(header)
