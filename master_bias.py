"""

Created on Sun Sep 25 20:55:01 2016

@author: Astro_Lub

"""
from astropy.io import fits
#import matplotlib.pyplot as plt
import numpy as np
#import robust as rb
import glob
#import img_scale

def master_bias(prefix="IMG",suffix=".FIT",dir="/Users/Julien/Dropbox/Test/"):

    """

    Args:
        prefix: constant between all images (non-numerical, decided by camera)
        suffix: file type
        dir: directory

    Returns:

    """
    imnames = glob.glob(dir+prefix+"*"+suffix)
    numimg = len(imnames)
    hdu= fits.open(imnames[0])[0]
    xd = hdu.header['NAXIS1']
    yd = hdu.header['NAXIS2']
    cube = np.zeros((yd,xd,numimg))
    for i in range(numimg):
        cube[:,:,i] = fits.getdata(imnames[i])
    master = np.zeros((yd,xd))
    for i in range(yd):
        for j in range(xd):
            master[i,j] = np.median(cube[i,j,:])
    return master


#We have now the name of every image\
