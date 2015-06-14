import numpy as np
from length import length
import ipdb
#------------------------------------------------------------------------------
#+
# NAME:
#   djs_photcen
#
# PURPOSE:
#   Recenter an object position.
#
# CALLING SEQUENCE:
#   djs_photcen, xcen, ycen, image, $
#    [ calg=, cbox=, cmaxiter=, cmaxshift=, fwhm=, /fixfw, ceps=, qmaxshift= ]
#
# INPUTS:
#   xcen:       X center(s)
#   ycen:       Y center(s)
#   image:      2-dimensional image
#
# OPTIONAL INPUTS:
#   calg:       Centering algorithm.  Choose from iweight, gauss1, gauss2, none.
#                 iweight = intensity-weighted center, computed independently
#                           in both X and Y
#                 gauss1  = gaussian fit, including a constant term, computed
#                           independently in both X and Y
#                 gauss2  = 2D gaussian fit, including a constant term
#                 none    = no centering
#               Default to iweight.
#   cbox:       Centering box width.  Default to 7.
#   cmaxiter:   Maximum number of iterations for centering algorithm.
#               Default to 10.
#   cmaxshift:  Maximum center shift.  If this shift is exceeded in either
#               X or Y, then return the center position passed in XCEN, YCEN.
#               A value of 0 imposes no maximum shift.  Default to 3.
#   fwhm:       FWHM for gauss1 or gauss2 centering algorithms.  If this is
#               a scalar, then the same value is used for both X and Y.
#               If this is a vector, then the first two elements are used
#               for X and Y respectively.
#   fixfw:      If set and nonzero, then fix the FWHM for gauss1 or gauss2 fits.
#   ceps:       Stop iterating when relative shifts in both X and Y are less
#               than this value.  Default to 0.
#
# OUTPUTS:
#   xcen:       Re-centered X position (modified).
#   ycen:       Re-centered X position (modified).
#
# OPTIONAL OUTPUS:
#   qmaxshift:  Return 1 if maximum shift was reached in either X or Y.
#               Return 0 otherwise.
#
# PROCEDURES CALLED:
#   curvefit()
#   np.ceil()
#   np.floor()
#   gauss2dfit()
#
# REVISION HISTORY:
#   01-Dec-1996  Written by D. Schlegel, Durham.
#   10-Aug-1998  Added option for calg='gauss2' (DJS).
#   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
#                of FITS headers; make IDL 5 complient (DJS).
#-
#------------------------------------------------------------------------------
# INTERNAL SUPPORT PROCEDURES:
#
# djs_photcen_gfunc1
# djs_photcen_gfit1()
#------------------------------------------------------------------------------
# Function of the form
#   F(x) = a * exp(-e^2) + c
# where e=(x-x0)/sigma, and the fitting parameters avec=[a,x0,c,sigma].
# This function can be passed to curvefit().

def djs_photcen_gfunc1(xvec, avec, fval, pderiv):
   n_params = 4
   def _ret():  return (xvec, avec, fval, pderiv)
   
   global fwhcom
   
   if (array(avec, copy=0).nelements() <= 3):   
      fixfw = 1
   else:   
      fixfw = 0
   
   if (fixfw == 0):   
      fwhcom = avec[3]
   
   ee = (xvec - avec[1]) / fwhcom
   ee2 = ee * ee
   bx = exp(-0.5 * ee2)
   fval = avec[0] * bx + avec[2]
   
   if (n_params >= 4):   
      # Calculate partial derivatives
      dfda = bx
      dfdx0 = bx * ee / fwhcom
      dfdc = 0 * bx + 1.0
      pderiv = concatenate([concatenate([bx]), concatenate([dfdx0]), concatenate([dfdc])])
      if (fixfw == 0):    #  FWHM is not fixed, so include it in fit
         dfdsig = dfdx0 * ee
         pderiv = concatenate([concatenate([pderiv]), concatenate([dfdsig])])
   
   
   return _ret()


#------------------------------------------------------------------------------
def djs_photcen_gfit1(xvec, yvec, wvec, fwhm, xcen, fixfw):

   n_params = 6
   
   global fwhcom
   
   # If too few data points, then return with no change to xcen
   if (array(xvec, copy=0).nelements() < 4):   
      print 'Error - Too few data points for gaussian fit'
      return xcen
   
   # Determine initial guesses
   c0 = array(yvec, copy=0).min()
   sig0 = fwhm / 2.354820
   a0 = max(yvec, imax) - c0
   x0 = xvec[imax]
   avec = concatenate([a0, x0, c0])
   
   if (fixfw != 0):   
      fwhcom = fwhm
   else:   
      avec = concatenate([avec, sig0])
   
   yfit = curvefit(xvec, yvec, wvec, avec, function_name='djs_photcen_gfunc1')
   
   newx = avec[1]
   return newx


#------------------------------------------------------------------------------
def djs_photcen(xcen, ycen, image, calg=None, cbox=None, cmaxiter=None, cmaxshift=None, fwhm=None, fixfw=None, ceps=None, qmaxshift=None):
   """Syntax - result = djs_photcen( xcen, ycen, image, $'
   [ calg=, cbox=, cmaxiter=, cmaxshift=, $'
   fwhm=, /fixfw, ceps= ] )"""

   outputs={'xcen':None, 'ycen':None, 'qmaxshift':None}

# Need 3 parameters
#   _opt = (calg, cbox, cmaxiter, cmaxshift, fwhm, fixfw, ceps, qmaxshift)
#   def _ret():
#      _optrv = zip(_opt, [calg, cbox, cmaxiter, cmaxshift, fwhm, fixfw, ceps, qmaxshift])
#      _rv = [xcen, ycen, image]
#      _rv += [_o[1] for _o in _optrv if _o[0] is not None]
#      return tuple(_rv)
#   
   
   #----------
   # If XCEN and YCEN are arrays, then call this routine recursively
   
   nx = length(xcen)
   ny = length(ycen)
   if (nx != ny):   
      print 'Dimensions of NX and NY do not agree'
      pass
   
   if (nx >1):   
      #outputs['xcen']=np.zeros(length(xcen))
      #outputs['ycen']=np.zeros(length(ycen))
      #outputs['qmaxshift']=np.zeros(length(qmaxshift))

      for i in np.arange(0, (nx - 1)+(1)):
         xtmp = xcen[i]
         ytmp = ycen[i]
         outs=djs_photcen(xtmp, ytmp, image, calg=calg,
                          cbox=cbox, cmaxiter=cmaxiter,
                          cmaxshift=cmaxshift, fwhm=fwhm,
                          fixfw=fixfw, ceps=ceps)
         
         for key in outputs.keys():
            if outputs[key] is None:
                outputs[key] = outs[key]
            else:
                outputs[key] = np.append(
                                   np.array(outputs[key]),
                                    outs[key])
      return outputs 
   
   if ((calg is not None) == 0):   
      calg = 'iweight'
   if ((cbox is not None) == 0):   
      cbox = 7
   if ((cmaxiter is not None) == 0):   
      cmaxiter = 10
   if ((cmaxshift is not None) == 0):   
      cmaxshift = 3.0
   if ((ceps is not None) == 0):   
      ceps = 0.0
   
   if ((fwhm is not None) == 0):   
      fwhvec = np.array([1.0, 1.0])
   else:   
      if ((fwhm is not None) == 1):   
         fwhvec = np.array([fwhm, fwhm])
      else:   
         fwhvec = np.array([fwhm(0), fwhm(1)])
   if ((fixfw is not None) == 0):   
      fixfw = 0
   # Return if no centering is to be performed
   if (calg == 'none'):   
      return outputs 
   
   # Use the data array dimensions
   dims = np.shape(image)
   naxis1 = dims[0]
   naxis2 = dims[1]
   
   radius = 0.5 * cbox
   
   # Iterate until one of the following conditions is met:
   #   (1) the maximum number of iterations is reached
   #   (2) no pixels are within the centering box
   #   (3) the maximum shift in both X and Y has been exceeded
   #   (4) the same center position is returned after two iterations,
   #       to within ceps of each other in both X and Y
   iiter = 0
   dcen = 2 * ceps + 1.0 # Set > ceps to prevent stopping at first iteration
   qmaxshift = 0
   while ( (iiter <= cmaxiter) &
           (qmaxshift == 0) &
           (np.max(abs(dcen)) > ceps)):
      if (iiter > 0):   
         dcen = np.array([xcen, ycen])
      
      # Limit computations to pixels in a square centered on (xCen, yCen)
      # which bounds the box with radius Radius.
      
      # Trim the X radius to be the largest symmetric box that falls within
      # the image boundaries.
      xrad = min(
                        radius, 
                        np.maximum(xcen, 0),
                        np.maximum((naxis1 - xcen),0))
         
      istart = np.floor(xcen + 0.5 - xrad).astype('int')
      iend = np.ceil(xcen - 0.5 + xrad).astype('int')
      if ( (istart > naxis1) | (iend < 0)):
         #maybe throw a sys.exit in python?
         print 'Error - No pixels in X range'
         return outputs
      
      ilen = iend - istart + 1
      
      # Trim the Y radius to be the largest symmetric box that falls within
      # the image boundaries.
      yrad = min(
                        radius, 
                        np.maximum(ycen, 0),
                        np.maximum((naxis2 - ycen),0))
      jstart = np.floor(ycen + 0.5 - yrad).astype('int')
      jend = np.ceil(ycen - 0.5 + yrad).astype('int')
      if ( (jstart > naxis2) | (jend < 0)):
         print 'Error - No pixels in Y range'
         return outputs
      jlen = jend - jstart + 1
      
      # Compute pixel corner locations relative to (xCen,yCen)
      xa = istart + np.arange(ilen) - 0.5 - xcen
      xb = xa + 1
      ya = jstart + np.arange(jlen) - 0.5 - ycen
      yb = ya + 1
      
      # Trim pixel sizes to only those coordinates falling within the
      # centering box
      xa = xa.clip(-xrad)
      xb = xb.clip(max=xrad)
      ya = ya.clip(-yrad)
      yb = yb.clip(max=yrad)

      # Determine the indices for the pixel
      npixelx = length(xa)
      npixely = length(ya)
      npixels = npixelx * npixely
      iindx   = (np.arange(npixelx*npixely, dtype=np.uint32)
                        .reshape(npixelx, npixely))%npixelx
       
      jindx   = (np.arange(npixelx*npixely, dtype=np.uint32)
                        .reshape(npixelx, npixely))/npixelx
      #iindx = lindgen(npixelx, npixely) % npixelx
      #jindx = array(lindgen(npixelx, npixely) / npixelx, copy=0).astype(Int32) + 0
      xpixnum = (iindx + istart)
      ypixnum = (jindx + jstart)
      #pixnum = 0 + xpixnum + naxis1 * ypixnum
      #(this line is commented out and simplified in 
      #the python version)
      # Compute contribution of each pixel
      fracs = (xb[iindx[:]] - xa[iindx[:]]) * (yb[jindx[:]] - ya[jindx[:]])
      
      # Create a temporary data array with the box of pixels
      subimg = image[xpixnum, ypixnum]
      
      _expr = calg
      
      #This section needs to fix the gaussian fits
      
      if _expr == 'iweight':   
         # Intensity-weighted centering
         norm = np.sum(subimg[:] * fracs)
         
         # Insist that the total flux is positive
         if (norm > 0):   
            # Work with xPixNum-xcen and yPixNum-ycen for numerical stability
            newi = np.sum(subimg[:] * (xpixnum[:] - xcen) * fracs) / norm + xcen
            newj = np.sum(subimg[:] * (ypixnum[:] - ycen) * fracs) / norm + ycen
            xcen = newi
            ycen = newj


      elif _expr == 'gauss1':   
         # THIS NEEDS TO BE IMPLEMENTED
         # One-dimensional gaussian fit, indepently fit in X and Y
         # Collapse the sub-image into two one-dimensional arrays
         
         # Fit X center
         ptemp = total(subimg, 2)
         fracx = xb - xa
         newi = djs_photcen_gfit1(xpixnum[0,:], ptemp, fracx, fwhvec[0], xcen, fixfw)
         xcen = newi
         
         # Fit Y center
         ptemp = total(subimg, 1)
         fracy = yb - ya
         newj = djs_photcen_gfit1(transpose(ypixnum[:,0]), ptemp, fracy, fwhvec[1], ycen, fixfw)
         ycen = newj
         
      elif _expr == 'gauss2':   
         # THIS NEEDS TO BE IMPLEMENTED
         # Two-dimensional gaussian fit, simultaneously fitting in X and Y
         
         gres = gauss2dfit(subimg, gcoeff, tilt=True)
         xcen = gcoeff[4] + istart
         ycen = gcoeff[5] + jstart
         
      elif _expr == none:   
         asdf = -1
         
      else:   
         print 'Error - Invalid centering algorithm'
         return _ret()
      
      
      # Test to see if maximum shift was reached in either X or Y.
      # If so, then reset coordinates to original center and set the
      # qmaxshift flag.
      if ((cmaxshift is not None)):   
         xold = xcen
         yold = ycen
         xcen = max(min(xcen, xcen+cmaxshift), xcen-cmaxshift)
         ycen = max(min(ycen, ycen+cmaxshift), ycen-cmaxshift)
         #xcen = maximum((minimum(xcen, xcen) + cmaxshift), xcen) - cmaxshift
         #ycen = maximum((minimum(ycen, ycen) + cmaxshift), ycen) - cmaxshift
         if ( (xold != xcen) | (yold != ycen)):
            outputs['qmaxshift'] = 1
      
      if (iiter > 0):   
         dcen = dcen - np.array([xcen, ycen])
      iiter = iiter + 1
      outputs['xcen']=xcen
      outputs['ycen']=ycen 
   
   return outputs

#------------------------------------------------------------------------------
