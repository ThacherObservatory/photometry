# -*- coding: utf-8 -*-
"""
Need to convert all the pyfits procedures to astropy procedures.

"""

# 001V image has an airplane
# Some of the brighter stars are saturated


from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import robust as rb
import glob
import img_scale

def readimage(imfile,plot=False,siglo=3,sighi=7):
    image,header = fits.getdata(imfile,0,header=True)
    med = np.median(image)
    sig = rb.std(image)
    if plot:
        plt.ion()
        plt.figure(1)
        vmin = med - siglo*sig
        vmax = med + sighi*sig
        plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gray')

    return image,header

def display_image(image,siglo=3,sighi=7,fignum=2):
    med = np.median(image)
    sig = rb.std(image)
    plt.ion()
    plt.figure(fignum)
    vmin = med - siglo*sig
    vmax = med + sighi*sig
    plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gray')
    return


def makeband(band='V'):

    files = glob.glob('Mantis*[0-9]'+band+'.fit')
    zsz = len(files)
    image0,header0 = readimage(files[0])
    ysz,xsz = np.shape(image)
    
    stack = np.zeros((xsz,ysz,zsz))
    for i in range(zsz):
        image,header = readimage(files[i])
        
        stack[:,:,i] = image
        
    final = np.median(stack,axis=2)

    if band == 'V':
        tag = 'Blue'
        
    if band == 'R':
        tag = 'Green'
    
    if band == 'ip':
        tag = 'Red'
        
    fits.writeto(tag+'.fit', final, header)
    
    return
    
        
    
    
def make_RGB(sigmax=3,sigmin=1,write=False):
    
    Blue,header = fits.getdata('Blue.fit',0,header=True)
    Green,header = fits.getdata('Green.fit',0,header=True)    
    Red,header = fits.getdata('Red.fit',0,header=True)
    
    bmed = np.median(Blue)
    gmed = np.median(Green)
    rmed = np.median(Red)
    
    bsig = rb.std(Blue)
    gsig = rb.std(Green)                
    rsig = rb.std(Red)
        
    final = np.zeros((Blue.shape[0],Blue.shape[1],3),dtype=float)
    
    final[:,:,0] = img_scale.sqrt(Red,scale_min=rmed+sigmin*rsig,scale_max=rmed+sigmax*rsig)
    final[:,:,1] = img_scale.sqrt(Green,scale_min=gmed+sigmin*gsig,scale_max=gmed+sigmax*gsig)
    final[:,:,2] = img_scale.sqrt(Blue,scale_min=bmed+sigmin*bsig,scale_max=bmed+sigmax*bsig)

    plt.ion()
    plt.figure(99)
    plt.imshow(final,aspect='equal')

    if write:
        plt.savefig('RGB.png',dpi=300)
    
    return

