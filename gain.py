# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 21:15:56 2016

@author: sara
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
#import robust as rb

def gain(flat1, flat2, masterdark):

    Sum = fits.getdata(flat1)+fits.getdata(flat2)-(fits.getdatat(masterdark)*2)
    Dif = fits.getdata(flat1)-fits.getdata(flat2)-(fits.getdata(masterdark)*2)

    Sumimg = Sum/2.0

    mean = np.mean(Sumimg)
    variance = np.std(Dif)**2/2.0
    gain = 1/np.polyfit(mean,variance,1)[0]

    return gain

flat1 = np.random.normal(20.7,7,1000) #mean, std, size
flat2 = np.random.normal(20.7,7,1000)
masterdark = np.random.normal
