import thacherphot as tp
import numpy as np
import matplotlib.pyplot as plt
from length import length
import quick_image as qi
import glob
from fitgaussian import *
from astropy.stats import sigma_clipped_stats
from photutils import daofind
from plot_params import *

plot_params()

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

    
plt.figure(78)
plt.clf()
mflux = np.mean(totflux)
for i in range(np.shape(cog)[1]):
    plt.plot(cog[:,i]*totflux[i]/mflux)

out = 'Normalization = '+str(mflux)+' counts'
plt.annotate(out,[0.87,0.6],horizontalalignment='right', \
             xycoords='figure fraction',fontsize='large')
plt.axhline(y=1.0,color='k',ls='--',lw=2)
plt.xlabel('Aperture radius (pixels)',fontsize=18)
plt.ylabel('Normalized Flux',fontsize=18)
plt.xlim(0,50)
plt.title('curves of growth for fiber E output',fontsize=20)
plt.savefig('fiber_raw4_curves.png',dpi=300)

