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
name='fiber_raw4*fit'
files = glob.glob(path+name)

def get_fiber_flux(file,fwhm=20,threshold=100):

    image,header = qi.readimage(file,plot=False)
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)

    source = daofind(image - median, fwhm=fwhm, threshold=threshold*std)
    xcen = np.int(np.round(source['xcentroid'][0]))
    ycen = np.int(np.round(source['ycentroid'][0]))

    plt.plot([xcen],[ycen],'rx',markersize=20)

    params = fitgaussian(image)
    
    y,x = params[1],params[2]

    apdict = tp.optimal_aperture(x,y,image,[150,160])
    return apdict


goodfiles = []
for file in files:
    phot = get_fiber_flux(file)
    if phot['totflux'] < 8e7:
        goodfiles = np.append(goodfiles,file)
        
totflux = []
cog = np.zeros((100,len(goodfiles)))
for i in range(len(goodfiles)):
    phot = get_fiber_flux(goodfiles[i])
    totflux = np.append(totflux,phot['totflux'])
    cog[:,i] = phot['curve_of_growth'][1]

for i in range(np.shape(cog)[1]):
    plt.plot(cog[:,i],label=str(totflux[i]))
plt.legend()
    
