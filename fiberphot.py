import thacherphot as tp
import numpy as np
import matplotlib.pyplot as plt
from length import length
import quick_image as qi
import glob
from fitgaussian import *
from astropy.stats import sigma_clipped_stats
from photutils import daofind


path = './'
files = glob.glob(path+"*fit")

image,header = qi.readimage(files[0],plot=True,siglo=2,sighi=2)

fwhm=10
threshold=10
mean, median, std = sigma_clipped_stats(image, sigma=3.0)
source = daofind(image - median, fwhm=fwhm, threshold=threshold*std)

ycen = np.int(np.round(source['xcentroid']))
xcen = np.int(np.round(source['ycentroid']))

subim = image[xcen-100:xcen+100,ycen-100:ycen+100]

qi.display_image(subim)

#biases,bct = tp.get_files(dir=path, tag='Bias')
#masterbias = tp.master_bias(biases,outdir='./')

#darks,dct = tp.get_files(dir = path,tag='Dark')
#masterdark = tp.master_dark(darks,bias=bias,outdir='./')

#cal = image-dark-bias

cal = subim

params = fitgaussian(cal)

y,x = params[1],params[2]

plt.plot([x],[y],'rx',markersize=20)

apdict = tp.optimal_aperture(x,y,cal,[20,25])

