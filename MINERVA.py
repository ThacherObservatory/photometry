import quick_image as qi
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


path = '/Users/jonswift/Astronomy/EBs/Data/MINERVA/2015Oct15/'

glob.glob(path+'*fits')

files = glob.glob(path+"*fits")
#qi.readimage(files[10],plot=True)

image,header = qi.readimage(files[50],plot=True,siglo=2,sighi=2)
w = wcs.WCS(header)

target = SkyCoord('19 51 39.82 +48 19 55.4', unit=(u.hourangle, u.deg))
world0 = np.array([[target.ra.degree, target.dec.degree]])
pix0 = w.wcs_world2pix(world0,1) 
x0 = pix0[0,0]
y0 = pix0[0,1]

plt.figure(1)
plt.plot([x0],[y0],'y+',markersize=30)
plt.xlim(0,2048)
plt.ylim(0,2048)


darks,dct = tp.get_files(tag='Dark')
biases,bct = tp.get_files(tag='Bias')
flats,fct = tp.get_files(tag='SkyFlat')
targetfiles,tct = tp.get_files(tag='KIC10935310')

bias = tp.master_bias(biases,outdir='./Photometry/',readnoise=False)

dark = tp.master_dark(darks,bias=bias,outdir='./Photometry/')

flat = tp.master_flat(flats,bias=bias,dark=dark,outdir='./Photometry/')

cal = (image - dark - bias)/flat
qi.display_image(cal,siglo=2,sighi=2)
plt.plot([x0],[y0],'y+',markersize=30,linewidth=1e6)
plt.xlim(0,2048)
plt.ylim(0,2048)

# define the aperture
radius,ynew,xnew,fwhm,aspect,snrmax,totflux,totap,chisq = \
    tp.optimal_aperture(x0,y0,image,skyrad=[15,20])

pref = np.array([[x0,y0],[836,1154],[645,1152], [1234,1002], [1131,1338], [377,426]])
coords0 = w.wcs_pix2world(pref,1) 
ras = coords0[:,0]
decs = coords0[:,1]

info = tp.batch_phot(targetfiles[:-1],ras,decs,bias=bias,dark=dark,flat=flat,
                     outdir='./Photometry/',skyrad=np.array([15,20]))

plt.figure(56)
plt.plot(info['jd']-np.median(info['jd']),info['flux'][:,0],'o-')

sys.exit()



tp.lightcurve(info)

position = (xnew,ynew)



aperture = CircularAperture(position, r=radius)

bkg_aperture = CircularAnnulus(position, r_in=15., r_out=20.)

# perform the photometry; the default method is 'exact'
phot = aperture_photometry(image, aperture)
bkg = aperture_photometry(image, bkg_aperture)

# calculate the mean background level (per pixel) in the annuli
bkg_mean = bkg['aperture_sum'] / bkg_aperture.area()
bkg_mean
bkg_sum = bkg_mean * aperture.area()

# plot the apertures
plt.imshow(scale_image(image, scale='sqrt', percent=98.), 
           origin='lower',cmap='gray')
aperture.plot(color='blue')
bkg_aperture.plot(color='cyan', hatch='//', alpha=0.8)
plt.xlim(xnew-100,xnew+100)
plt.ylim(ynew-100,ynew+100)

phot['bkg_sum'] = bkg_sum
