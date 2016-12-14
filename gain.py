# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 21:15:56 2016

@author: sara
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import robust as rb
import glob

def gain(flat1, flat2, center=[50,50], sidelen=10):#(masterdark?)

    xcent = center[0]
    ycent = center[1]
    
    Sum = fits.getdata(flat1, header=False)+fits.getdata(flat2,header=False)#-(fits.getdatat(masterdark)*2)
    Dif = fits.getdata(flat1,header=False)-fits.getdata(flat2,header=False)#-(fits.getdata(masterdark)*2)

    Sumimg = Sum/2.0
    mean = np.mean(Sumimg[ycent-(sidelen/2):ycent+(sidelen/2),xcent-(sidelen/2):xcent+(sidelen/2)])
    variance = np.std(Dif[ycent-(sidelen/2):ycent+(sidelen/2),xcent-(sidelen/2):xcent+(sidelen/2)])**2/2.0
    #gain = 1/np.polyfit(mean,variance,1)[0]

    return mean,variance

#flat1 = np.random.normal(20.7,7,1000) #mean, std, size
#flat2 = np.random.normal(20.7,7,1000)
#masterdark = np.random.normal
files = glob.glob('/Users/georgelawrence/python/Astronomy/gaintests/g*')
mean = []
variance = []
for ind,f in enumerate(files):
    try:
        mean.append(gain(f,files[ind+1])[0])
        variance.append(gain(f,files[ind+1])[1])
    except:
        mean.append(gain(f,files[0])[0])
        variance.append(gain(f,files[0])[1])
        
plt.figure(1)
plt.clf()
plt.ion()
plt.plot(mean,variance,'ro')
plt.show()
print 1/np.polyfit(mean,variance,1)[0]