"""

Created on Sun Sep 25 20:55:01 2016

@author: Astro_Lub,syao,george

"""
from astropy.io import fits
#import matplotlib.pyplot as plt
import numpy as np
#import robust as rb
import glob
#import img_scale
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
