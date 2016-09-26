# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 14:04:05 2016

@author: sara
"""
<<<<<<< HEAD
import numpy as np
import matplotlib.pyplot as plt
=======
import matplotlib.animation as animation
import numpy as np
import matplotlib.pyplot as plt
import pyfits as pf
>>>>>>> db16ccd1d282e41ad5622a64865ae5e7191d279c
from fitgaussian import *
import robust as rb
import sys, os, time, glob
#import scipy as sp
#import matplotlib.patheffects as PathEffects
import pdb
#import djs_phot_mb as djs
#from select import select
#from astropysics.coords import AngularCoordinate as angcor
#import astropy.io.fits
#from astropy import wcs
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import robust as rb
import glob
import img_scale
<<<<<<< HEAD
#import quick_image
#import allsky
=======
import quick_image.py
import allsky.py
>>>>>>> db16ccd1d282e41ad5622a64865ae5e7191d279c

"""
Header Info:
Keyword	Description
SIMPLE	Set to T if the file conforms to the FITS standard and to F if not.
BITPIX	Datatype of the data in the HDU.
XTENSION	Type of extension: IMAGE, TABLE, BINTABLE.
NAXIS	Number of dimensions of a data array. For tables this is always 2.
NAXISn	Length of each dimension. Here 1 <= n <= NAXIS. For tables NAXIS1 is the number of table rows.
TFIELDS	Number of fields, i.e., table columns, in the table.
TFORMn	Format of data in the nth table column. See table 7.3 in page 50 of the FITS Standard.
BSCALE	Factor by which data in an image has been scaled.
BZERO	Zero point of the data in an image array.
TSCALEn	Same as BSCALE, but for table column n.
TZEROn	Same as BZERO, but for table column n.
"""
<<<<<<< HEAD
#----------------------------------------------------------------------#
# get_files:                                                           #
#----------------------------------------------------------------------#

def get_files(prefix='IMG',dir="/Users/jonswift/Dropbox (Thacher)/Observatory/AllSkyCam/Data/",
              suffix='.FIT'):
    
    """
    Overview:
    ---------
    Returns list of files with a user defined prefix and suffix withing a
    specified directory
    
    
    Calling sequence:
    -----------------
    files = get_files('HATp33b',dir='/home/users/bob/stuff/')
    
    """

    files = glob.glob(dir+prefix+"*"+suffix)
    
    fct = len(files)

    return files,fct

'''def listsum(List):
    theSum = 0
    for i in List:
        theSum = theSum + i
    return theSum'''

def pixel_median(xp,yp,prefix='IMG',dir="/Users/jonswift/Dropbox (Thacher)/Observatory/AllSkyCam/Data/",
              suffix='.FIT'):
    files, files_len = get_files(prefix,dir,suffix)
    master_bias = fits.PrimaryHDU()

    pixle_median  = []     
    for image in files:
         hdu = fits.open(image)[0]
         xd = hdu.header['NAXIS1']
         yd = hdu.header['NAXIS2'] 
         datarow = hdu.data[xp]
         datavalue = datarow[yp]
         pixle_median.append(datavalue)
    
    Median = np.median(pixle_median)
    return Median
    
=======
def master_bias(prefix='IMG',dir="/Users/jonswift/Dropbox (Thacher)/Observatory/AllSkyCam/Data/",
              suffix='.FIT'):
    files, files_len = get_files(prefix,dir,suffix)
    
    for image in files:
        


>>>>>>> db16ccd1d282e41ad5622a64865ae5e7191d279c


























    
    