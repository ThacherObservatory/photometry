#from __future__ import division
from length import length
import numpy as np
import ipdb
#------------------------------------------------------------------------------
#+
# NAME:
#   djs_photfrac
#
# PURPOSE:
#   Create a list of pixel numbers and their fractional contribution to
#   an annular region.
#
# CALLING SEQUENCE:
#   djs_photfrac, xcen, ycen, Rvec, xdimen=, ydimen=, $
#    [ xPixNum=, yPixNum=, pixnum=, fracs=, fillfrac= ]
#
# INPUTS:
#   xcen:       X center(s)
#   ycen:       Y center(s)
#   Rvec:       Either a 2-element array with two radii to define an annulus,
#               or a scalar to define a circular aperature.
#
# OPTIONAL INPUTS:
#   xdimen:     Number of X pixels.
#   ydimen:     Number of Y pixels.
#
# OUTPUTS:
#   pixnum:     Pixel number, 0-indexed, for referencing array using one index.
#   xPixNum:    Pixel number in X, 0-indexed.
#   yPixNum:    Pixel number in Y, 0-indexed.
#   fracs:      Return value of covering fraction of the annulus
#               over the pixel number.
#   fillfrac:   Ratio of returned pixel areas to the annulus area;
#               this ratio should equal 1.0 if the aperature falls completely
#               within the image boundaries
#
# COMMENTS:
#   The total counts within this region is given by
#     totcounts = total( pData(pixnum) * fracs )
#   The area within this region is given by
#     area = total(fracs)
#   The average counts is given by
#     totcounts = total( pData(pixnum) * fracs ) / total(fracs)
#   To test for bad pixels, e.g. values greater than vmax, within
#   the aperature,
#     if (where(pData(pixnum) GT vmax) EQ -1) then <no bad values> $
#     else <bad values exist>
#
#   If no pixels within the given annulus are found, then return pixnum=-1.
#
# BUGS:
#   - can wrap around on edge of you use PixNum.  XPixNum,YPixNum do
#     not exhibit this problem
#
# PROCEDURES CALLED:
#   djs_ceil()
#   djs_floor()
#
# REVISION HISTORY:
#   Written D. Schlegel, 27 November 1996, Durham
#   Bug identified - 2 Nov 2000, D. Finkbeiner
#-
#------------------------------------------------------------------------------
# INTERNAL SUPPORT PROCEDURES:
#
# djs_photfrac_intcirc
#------------------------------------------------------------------------------
# Return the integral of the circle with radius Radius within the pixel
# defined by x=[xA,xB] and y=[yA,yB].  Add to this the area beneath the
# pixel, i.e. the value (xB-xA)*yA.
#
# This function assumes that xB > xA >= 0 and yB > yA >= 0, and that
# the circle passes through the pixel:
#   xA^2 + yA^2 <= Radius^2 <= xB^2 + yB^2

def djs_photfrac_intcirc(xa, xb, ya, yb, radius):
   xanorm = xa / radius
   xbnorm = xb / radius
   yanorm = ya / radius
   ybnorm = yb / radius
   
   
   gg = 0.0 * xanorm
   pix1 = np.nonzero((ybnorm >= 1.0))
   count = np.size(pix1)
   if (count != 0):   
      gg[pix1] = xanorm[pix1]

   pix2 = np.nonzero(ybnorm < 1.0)
   count = np.size(pix2)

   if (count != 0):   
      gg[pix2] = np.maximum(xanorm[pix2], np.sqrt(1.0 - ybnorm[pix2] * ybnorm[pix2]))
   
   hh = np.minimum(xbnorm, np.sqrt(1.0 - yanorm * yanorm))
   
   result = radius * radius * ((gg - xanorm) * ybnorm + (xbnorm - hh) * yanorm + 0.5 * (hh * np.sqrt(1.0 - hh * hh) + np.arcsin(hh) - gg * np.sqrt(1.0 - gg * gg) - np.arcsin(gg)))
   return result

#------------------------------------------------------------------------------
def djs_photfrac(xcen, ycen, rvec, xdimen=None, ydimen=None, xpixnum=None, ypixnum=None, pixnum=None, fracs=None, fillfrac=None):
   """Syntax - djs_photfrac, xcen, ycen, Rvec, $
    xdimen=, ydimen=, xPixNum=, yPixNum=, pixnum=, $
    fracs=, fillfrac="""
  
   outputs = {'pixnum':None,
              'xpixnum':None, 'ypixnum':None,
              'fracs':None, 'fillfrac':None}
   
   pixnum = -1
   fracs = 0
   
   # If Rvec contains one element, then use the annulus [0,Rvec],
   # otherwise use the annulus [Rvec[0],Rvec[1]].
   if (length(rvec)==1):
      radius1 = 0.0
      radius2 = abs(rvec)
   else:   
      radius1 = abs(rvec[0])
      radius2 = abs(rvec[1])
   
   sqradius1 = radius1 * radius1
   sqradius2 = radius2 * radius2
   
   # Limit computations to pixels in a square centered on (xCen, yCen)
   # which completely bounds the outer circle described by Radius2.
   istart = np.maximum(0, 
                np.floor(xcen + 0.5 - radius2).astype('int'))
   iend = np.minimum((xdimen - 1), 
                np.ceil(xcen - 0.5 + radius2).astype('int'))
   ilen = iend - istart + 1
   if ((istart > xdimen) |
       (iend < 0) |
       (ilen < 1) ):
      print 'Error - No pixels in X range'
      pass
   
   jstart = np.maximum(0, 
                np.floor(ycen + 0.5 - radius2).astype('int'))
   jend = np.minimum((ydimen - 1), 
                np.ceil(ycen - 0.5 + radius2).astype('int'))
   jlen = jend - jstart + 1
   if ((jstart > ydimen) |
       (jend < 0) |
       (jlen < 1)):
      print 'Error - No pixels in Y range' 
      pass

   # Compute pixel corner locations relative to (xCen,yCen)
   xa = (istart + np.arange(ilen) - 0.5) - xcen
   xb = xa + 1
   ya = (jstart + np.arange(jlen) - 0.5) - ycen
   yb = ya + 1
   
   # Split the pixel across the axis y=yCen by adding one more X pixel
   # e.g., by adding one more element to xA and xB.
   isplit = np.floor(xcen - istart + 0.5).astype('int')# Integer X pixel number to split
   

   q_splitx = (isplit >= 0) & (isplit < ilen)
   
   if (q_splitx):   
      xa = np.append(xa,0.0)
      xb = np.append(xb,xb[isplit])
      xb[isplit] = 0.0
   
   # Split the pixel across the axis x=xCen by adding one more Y pixel
   # e.g., by adding one more element to yA and yB.
   # Integer Y pixel number to split
   jsplit = np.floor(ycen - jstart + 0.5).astype('int')
   q_splity = (jsplit >= 0) &( jsplit < jlen)
   if (q_splity):   
      ya = np.append(ya,0.0)
      yb = np.append(yb, yb[jsplit])
      yb[jsplit] = 0.0
   
   # Force all relative coordinates to be non-negative and reorder
   # values such that xB>xA and yB>yA.
   xa = np.abs(xa)
   xb = np.abs(xb)
   ya = np.abs(ya)
   yb = np.abs(yb)
   xap = np.minimum(xa, xb)
   xbp = np.maximum(xa, xb)
   yap = np.minimum(ya, yb)
   ybp = np.maximum(ya, yb)
   
   # Compute distances to lower left corner of pixel, RadiusA,
   # and upper right corner, RadiusB.
   npixelx = length(xap)
   npixely = length(yap)
   npixels = npixelx * npixely
   
   #iindx   = (np.arange(npixelx*npixely, dtype=np.uint32).reshape(npixelx, npixely))%npixelx
   #jindx   = (np.arange(npixelx*npixely, dtype=np.uint32).reshape(npixelx, npixely))/npixelx
   iindx, jindx = np.meshgrid(
                    np.arange(npixelx),
                    np.arange(npixely))
   
   sqradiusa = (xap[iindx.flatten()].reshape(np.shape(iindx))**2+
   yap[jindx.flatten()].reshape(np.shape(jindx))**2) 
   sqradiusb = (xbp[iindx.flatten()].reshape(np.shape(iindx))**2+
   ybp[jindx.flatten()].reshape(np.shape(jindx))**2 )
   
   # Integrate within the annulus defined by the circles [Radius1,Radius2]
   # within each pixel.
   integral = np.zeros((npixely, npixelx))
   
   qpix0 = (sqradiusb > sqradius1) & (sqradiusa < sqradius2)
   pix1 = np.nonzero( (sqradiusb < sqradius2) & qpix0==True)
   count=np.size(pix1)
   if (count != 0):   
      integral[pix1] = (xbp[iindx[pix1]] - xap[iindx[pix1]]) * ybp[jindx[pix1]]
   
   pix2 = np.nonzero((sqradiusb > sqradius2) & qpix0==True)
   count=np.size(pix2)
   if (count != 0):   
      integral[pix2] = djs_photfrac_intcirc(xap[iindx[pix2]], xbp[iindx[pix2]], yap[jindx[pix2]], ybp[jindx[pix2]], radius2)
   
   pix1 = np.nonzero((sqradiusa >= sqradius1) & qpix0==True)
   count= np.size(pix1)
   if (count != 0):   
      integral[pix1] = integral[pix1] - (xbp[iindx[pix1]] - xap[iindx[pix1]]) * yap[jindx[pix1]]
   
   pix2 = np.nonzero((sqradiusa < sqradius1) & qpix0==True)
   count = np.size(pix2)
   if (count != 0):   
      integral[pix2] = integral[pix2] - djs_photfrac_intcirc(xap[iindx[pix2]], xbp[iindx[pix2]], yap[jindx[pix2]], ybp[jindx[pix2]], radius1)
   
   # Collapse the split pixels back into the original pixels
   if (q_splity):   
      integral[jsplit,0:(npixelx - 1)+1] = integral[jsplit,0:(npixelx - 1)+1] + integral[npixely - 1,0:(npixelx - 1)+1]
   if (q_splitx):   
      integral[0:(npixely - 1)+1,isplit] = integral[0:(npixely - 1)+1,isplit] + integral[0:(npixely - 1)+1,npixelx - 1]
   
   # Set the return values
   xpixnum = iindx[0:(npixely - 1 - q_splity)+1,0:(npixelx - 1 - q_splitx)+1] + istart
   ypixnum = jindx[0:(npixely - 1 - q_splity)+1,0:(npixelx - 1 - q_splitx)+1] + jstart
   fracs = integral[0:(npixely - 1 - q_splity)+1,0:(npixelx - 1 - q_splitx)+1]
   # Limit the return values to only those pixels with non-zero contributions
   pix1 = np.nonzero(fracs != 0.0)
   if (length(pix1)==0):   
      return  # ??? Should never happen!
   xpixnum = xpixnum[pix1]
   ypixnum = ypixnum[pix1]
   pixnum = 0 + xpixnum + xdimen * ypixnum
   fracs = fracs[pix1]
   
   # Test to see if aperature exceeds image boundary by computing the
   # ratio of the filled pixels to the area of the annulus.  If all
   # pixels are within the image boundary, then fillfrac=1.0.
   fillfrac = np.sum(fracs) / (np.pi * (sqradius2 - sqradius1))
   
   outputs = {'pixnum':pixnum,
              'xpixnum':xpixnum, 'ypixnum':ypixnum,
              'fracs':fracs, 'fillfrac':fillfrac}
   return outputs

#------------------------------------------------------------------------------
