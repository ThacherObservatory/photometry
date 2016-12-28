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
from plot_params import *

def gain(file1, file2, center=[1300,1300], sidelen=100):#(masterdark?)
    xcent = center[0]
    ycent = center[1]
    
    flat1 = np.float64(fits.getdata(file1, header=False))
    flat2 = np.float64(fits.getdata(file2,header=False))
    Sum = (flat1+flat2)/2
    Dif = flat1-flat2

    mean = np.mean(Sum[ycent-(sidelen/2):ycent+(sidelen/2),xcent-(sidelen/2):xcent+(sidelen/2)])
    variance = np.std(Dif[ycent-(sidelen/2):ycent+(sidelen/2),xcent-(sidelen/2):xcent+(sidelen/2)])**2/2.0

    return mean,variance


inttimes = [2,4,6,8,10,12,14,16,18]

mean = []
meanerr = []
variance = []
varerr = []

for int in inttimes:
    print '... starting gain calculations for integration time of '+str(int)+ 'seconds'
    print '--------------------------------------------------'
    files = glob.glob('/Users/jonswift/Dropbox (Thacher)/Observatory/Data/Science/16Dec2016/gain*_'+str(int)+'.fit')

    mntmp = []
    vartmp = []
    numfiles = len(files)
    for i in 2*np.arange(numfiles/2):
        print files[i].split('/')[-1],files[i+1].split('/')[-1]
        print gain(files[i],files[i+1])
        print ''
        mntmp = np.append(mntmp,gain(files[i],files[i+1])[0])
        vartmp = np.append(vartmp,gain(files[i],files[i+1])[1])
    mean = np.append(mean,np.mean(mntmp))
    meanerr = np.append(meanerr,np.std(mntmp))
    variance = np.append(variance,np.mean(vartmp))
    varerr = np.append(varerr,np.std(vartmp))

plot_params()
plt.figure(1)
plt.clf()
plt.ion()
plt.errorbar(mean,variance,xerr=meanerr,yerr=varerr,fmt='o',color='k',capthick=1.5)
plt.xlabel('Mean Value (counts)',fontsize=18)
plt.ylabel('Variance',fontsize=18)
fit,cov = np.polyfit(mean,variance,1,cov=True)
gain  = 1/fit[0]
gainerr = np.sqrt(cov[0,0])
plt.annotate('gain = %.3f $\pm$ %.3f e$^-$/ADU'%(gain,gainerr),[0.15,0.8],
             horizontalalignment='left',xycoords='figure fraction',fontsize=14)
x = np.linspace(np.min(mean),np.max(mean),100)
y = np.polyval(fit,x)
plt.plot(x,y,'r--',lw=1.5)
#plt.savefig('gain.png',dpi=300)
plt.show()


#gain(f,files[ind+1])[0]<18000 and gain(f,files[ind+1])[0]>12000 and
#gain(f,files[0])[0]<18000 and gain(f,files[0])[0]>12000 and
