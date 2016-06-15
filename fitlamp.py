# -*- coding: utf-8 -*-
"""
Created on Tue Jun 14 16:25:34 2016

@author: ONeill
"""
# Purpose: 1) Extrapolate baseline flux of MINERVA light source
#          2) Compare to known flux of light source through fiber optic cable
#          3) Post cable repair, determine if repair improves fluxcable:fluxbase
# KO/JS 6/14/16: Drew from fiberphot.py and thacherphot.py, modified for purpose

import thacherphot as tp
import numpy as np
import matplotlib.pyplot as plt
from length import length
import quick_image as qi
import glob
from fitgaussian import *
from astropy.stats import sigma_clipped_stats

path = '/Users/ONeill/Dropbox/Observatory/Data/minerva/FiberTests/'
files = glob.glob(path+"lamp2*fit")

image,header = qi.readimage(files[0],plot=True,siglo=2,sighi=2)

params = fitgaussian(image)

#using params make model: amplitude, x center, y center, x width, y width, angle rotated
#using thacherphot to overplot contours of gaussian
#for every x y pixel have data value and model value, so subtract from each other and find residuals - should just look like noise

fit = gaussian(*params)
level = (np.max(image) - np.median(image))*np.array([0.95,0.5,0.1])
plt.figure(1)
plt.contour(fit(*indices(image.shape)),level,colors='blue')

def find_flux(image):
    model = fit(*indices(image.shape))
    qi.display_image(model)

    resid = model - image
    qi.display_image(resid)

    #plt.figure(2)
    #plt.hist(resid, bins=50)

    plt.figure(3)
    qi.display_image(resid)

    A = params[0]
    Xsig = params[3]
    Ysig = params[4]
    #Check if sigmas x and y = xcen and ycen

    flux = (A * 2 * pi * Xsig * Ysig)/10.
    
    num = (np.sum(image))/10.
    
    return num, flux

fluxes = []
totals = []

for file in files:
    image,header = qi.readimage(file)
    tot, flux = find_flux(image)
    fluxes = np.append(fluxes,flux)
    totals = np.append(totals,tot)