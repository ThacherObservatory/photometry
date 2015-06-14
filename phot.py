from __future__ import division
from __future__ import print_function
from length import length
import numpy as np
#import ipdb

#------------------------------------------------------------------------------
#+
# NAME:
#   phot
#
# PURPOSE:
#   Driver for aperture photometry with the option of re-centering and
#   sky-subtraction.
#
# CALLING SEQUENCE:
#   flux = phot( xcen, ycen, objrad, skyrad, image, [invvar, $
#    calg=, cbox=, cmaxiter=, cmaxshift=, $
#    fwhm=, fixfw=, ceps=, $
#    salg=, srejalg=, smaxiter=, $
#    lorej=, hirej=, $
#    flerr=, skyval=, skyrms=, skyerr=, peakval=, /quick, /exact ] )
#
# INPUTS:
#   xcen:       X center(s)
#   ycen:       Y center(s)
#   objrad:     Radius for aperture on object, or a vector of such radii.
#   skyrad:     A 2-element array with two radii to define an annulus,
#               or a scalar to define a circular aperture, or undefined
#               for no sky calculation
#   image:      FITS data array, as read from readfits().
#
# OPTIONAL INPUTS:
#   invvar:     Inverse variance image, for computing the errors FLERR
#   ----------- FOR CENTERING ALGORITHM
#   calg:       Centering algorithm.  Choose from iweight, gauss1, gauss2, none.
#                 iweight = intensiAstronomy/Caltech/MINERVA/Observing/Photometry/ty-weighted center, computed independently
#                           in both X and Y
#                 gauss1  = gaussian fit, including a constant term, computed
#                           independently in both X and Y
#                 gauss2  = not implemented
#                 none    = no centering
#               Default to iweight.
#   cbox:       Centering box width.  Default to 7.
#   cmaxiter:   Maximum number of iterations for centering algorithm.
#               Default to 10.
#   cmaxshift:  Maximum center shift.  If this shift is exceeded in either
#               X or Y, then return the center position passed in XCEN,YCEN.
#               A value of 0 imposes no maximum shift.  Default to 3.
#   fwhm:       FWHM for gauss1 or gauss2 centering algorithms.  If this is
#               a scalar, then the same value is used for both X and Y.
#               If this is a vector, then the first two elements are used
#               for X and Y respectively.
#   fixfw:      If set and nonzero, then fix the FWHM for gauss1 or gauss2 fits.
#   ceps:       Stop iterating when relative shifts in both X and Y are less
#               than this value.  Default to 0.
#
#   ----------- FOR SKY FITTING ALGORITHM
#   salg:       Sky fitting algorithm.  Choose from mean, median, mode, none.
#               Default to "mean" if SKYRAD is set, or "none" otherwise.
#   srejalg:    Rejection algorithm.  Choose from none, sigclip, pclip.
#                 sigclip = sigma clipping; reject outliers that are
#                           more than lorej*sigma below skyval or hirej*sigma
#                           above skyval
#                 pclip   = percentile clipping; reject the lowest lorej
#                           fraction and the highest hirej fraction of points
#                 none    = no rejection
#               Default to sigclip
#   smaxiter:   Maximum number of srejalg iterations.  Default to 10.
#   lorej:      If srejalg="sigclip", then the number of standard deviations
#               below skyval to clip (default to 3.0).
#               If srejalg="pclip", then fraction of low pixels to clip
#               (default to 0.05).
#   hirej:      If srejalg="sigclip", then the number of standard deviations
#               above skyval to clip (default to 3.0).
#               If srejalg="pclip", then fraction of high pixels to clip
#               (default to 0.05).
#
# KEYWORDS:
#   exact       Use slow photo-based photfrac algorithm (Not thoroughly tested)
#               Requires image to be centered such that xcen and ycen
#               are integer values. If set, does not recalculate
#               center.
#   quick       Use faster photfrac algorithm (Not thoroughly tested)
#
# OUTPUTS:
#   flux:       Total flux within each circular aperture defined by objrad,
#               minus the sky contribution within each aperture [NOBJ,NRAD].
#   xcen:       Re-centered X position (modified).
#   ycen:       Re-centered X position (modified).
#
# OPTIONAL OUTPUTS:
#   flerr:      Flux error from the sky background uncertainty within
#               the object aperture(s) [NOBJ,NRAD]
#   skyval:     Sky value in counts per unit area [NOBJ].
#   skyrms:     RMS of sky pixels about skyval, again in counts per unit area.
#               This assumes that each unit area is independent.  The RMS
#               is computed for only the values that remain after all the
#               rejection iterations [NOBJ].
#   skyerr:     The error in skyval assuming that each pixel is independent.
#               This is skyrms divided by the square root of the number of
#               pixels [NOBJ].
#   peakval:    Peak pixel value (before sky-subtraction)
#
# COMMENTS:
#   Sub-pixel sampling of the circular apertures is handled exactly.
#   If /exact keyword is set, input xcen, ycen must be integers or
#     the code bombs. See exact_photfrac.pro for more details, but
#     basically exact_photfrac is much simpler to implement if the
#     object is already centered, which doesn't cost you precision.
#   For similar reasons, if /exact is set, phot will not try to
#     recentroid your object.
#
# REVISION HISTORY:
#   28-Nov-1996  Written by D. Schlegel, Durham.
#   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
#                of FITS headers; make IDL 5 compliant (DJS).
#   02-Nov-2000  objrad, skyrad recast as floats (D. Finkbeiner)
#                  If they are ints, 1% errors may arise.
#-
#------------------------------------------------------------------------------
def phot(xcen, ycen, objrad, skyrad, image, invvar=None, calg=None, 
             cbox=None, cmaxiter=None, cmaxshift=None, fwhm=None, fixfw=None, 
             ceps=None, salg=None, srejalg=None, smaxiter=None, lorej=None, 
             hirej=None, flerr=None, skyval=None, skyrms=None, skyerr=None, 
             peakval=None, quick=None, exact=None):
   outputs = {'flux':None, 'fluxerr':None,
              'xcen':None, 'ycen':None,
              'skyval':None,'skyrms':None}
   
   #convert everything to numpy arrays
   xcen     =np.array([xcen]).flatten()
   ycen     =np.array([ycen]).flatten()
   objrad   =np.array([objrad]).flatten()
   skyrad   =np.array([skyrad]).flatten()

   
   nobj = length(xcen)
  
   if (length(ycen) != nobj):   
      print('XCEN and YCEN must have same number of elements')
   
   if ((exact is not None)): 
      # This still needs to be fixed  
      ini = where(ravel(bitwise_or(array(xcen, copy=0).astype(Int32) != xcen, array(ycen, copy=0).astype(Int32) != ycen)))[0]

      if (nni > 0):   
         print('xcen and ycen MUST be integers if /exact keyword is set')
   
   dims = np.shape(image)
   xdimen = dims[0]
   ydimen = dims[1]
   nrad = length(objrad)
   
   # Allocate memory for return values
   flux = np.zeros([nrad, nobj])
   
   if (flerr is not None):   
      flerr = np.zeros([nrad, nobj])
   
   skyval = np.zeros([nobj])
   skyrms = np.zeros([nobj])
   skyerr = np.zeros([nobj])
   
   
   if (peakval is not None):   
      peakval = np.zeros([nrad, nobj])
   
   #----- LOOP THROUGH EACH OBJECT -----
   for iobj in range(0, (nobj - 1)+(1)):
   
   # Center the object
      xcen1 = xcen[iobj]
      ycen1 = ycen[iobj]
   
      if ((exact is None)):
         
         recenter=photcen(xcen1, ycen1, image, calg=calg, cbox=cbox, 
                          cmaxiter=cmaxiter, cmaxshift=cmaxshift)
         xcen1 = recenter['xcen']
         ycen1 = recenter['ycen']
      else:   
         print('Not recentering -- /exact requires correct input center')
         xcen1 = xcen[iobj]
         ycen1 = ycen[iobj]
      
      # Find the sky value
      # Add 0.0 to the skyrad to convert to floating-point
      
      skyvals= photsky(xcen1, ycen1, skyrad + 0.0, image, salg=salg, 
                       srejalg=srejalg, smaxiter=smaxiter, lorej=lorej, 
                       hirej=hirej, quick=quick, exact=exact)
      skyval[iobj] = skyvals['skyval']
      skyrms[iobj] = skyvals['skyrms']
      skyerr[iobj] = skyvals['skyerr']

      # Find total counts in object aperture(s)
      for irad in np.arange(0, (nrad - 1)+(1)):
         onerad = objrad[irad] + 0.0 # Convert to floating-point
         if (onerad > 0):   
            if (quick is not None):   
               quick_photfrac(xcen1, ycen1, onerad, xdimen=xdimen,
                              ydimen=ydimen, pixnum=pixnum, fracs=fracs)
            else:   
               if (exact is not None):   
                  exact_photfrac(xcen1, ycen1, onerad, xdimen=xdimen, 
                                 ydimen=ydimen, pixnum=pixnum, fracs=fracs)
               else:  
                  output_object=photfrac(xcen1, ycen1, onerad, 
                                             xdimen=xdimen, ydimen=ydimen)
                  pixnum=output_object['pixnum']
                  fillfrac=output_object['fillfrac']
                  xpixnum=output_object['xpixnum']
                  ypixnum=output_object['ypixnum']
                  fracs=output_object['fracs']
         else:   
            pixnum = -1
         if (pixnum[0] == -1):   
            flux[irad,iobj] = 0.0  # No object pixels
            if (flerr is not None):
               flerr[irad,iobj] = 0.0
         else:   
            area = np.sum(fracs)
            flux[irad,iobj] = np.sum(image[xpixnum,ypixnum] * fracs) - \
                                    (area * skyval[iobj])
            if (flerr is not None):
               flerr[irad,iobj] = area * skyerr[iobj]
            if (peakval is not None):
               peakval[irad,iobj] = np.max(image[pixnum])
         
         # If INVVAR was set, then measure those pixel errors and add
         # them in quadrature to the sky-subtraction errors.
         if ((invvar is not None) & (flerr is not None)):
            if (pixnum[0] == -1):   
               sigma2 = 0
            else:   
               if (min(invvar[xpixnum, ypixnum])<=0 ):
                  sigma2 = 0
               else:   
                  sigma2 = np.sum(fracs / invvar[xpixnum, ypixnum])
            if (sigma2 <= 0):   
               flerr[irad,iobj] = -1
            else:   
               flerr[irad,iobj] = np.sqrt(flerr[irad,iobj] ** 2 + sigma2)
      
      # Overwrite the center positions with the newly computed ones
      xcen[iobj] = xcen1
      ycen[iobj] = ycen1
      
      outputs['flux']=flux
      outputs['fluxerr']=flerr
      outputs['xcen'] = xcen
      outputs['ycen'] = ycen
      outputs['skyval']=skyval
      outputs['skyrms']=skyrms
   
   return outputs

#------------------------------------------------------------------------------
#+
# NAME:
#   photsky
#
# PURPOSE:
#   Compute the sky value within an annulus.
#
#   At present, fractional pixels are not treated properly; this is because
#   we need a sort routine that carries index numbers, such as NR indexx().
#
# CALLING SEQUENCE:
#   skyval = photsky( xcen, ycen, skyrad, image, $
#    [ salg=, srejalg=, smaxiter=, $
#    lorej=, hirej=, skyrms=, skyerr= ] )
#
# INPUTS:
#   xcen:       X center(s)
#   ycen:       Y center(s)
#   skyrad:     A 2-element array with two radii to define an annulus,
#               or a scalar to define a circular aperature, or undefined or 0
#               for no sky calculation
#   image:      FITS data array, as read from readfits().
#
# OPTIONAL INPUTS:
#   salg:       Sky fitting algorithm.  Choose from mean, median, mode, none.
#               Default to "mean" if SKYRAD is set, or "none" otherwise.
#   srejalg:    Rejection algorithm.  Choose from none, sigclip, pclip.
#                 sigclip = sigma clipping; reject outliers that are
#                           more than lorej*sigma below skyval or hirej*sigma
#                           above skyval
#                 pclip   = percentile clipping; reject the lowest lorej
#                           fraction and the highest hirej fraction of points
#                 none    = no rejection
#               Default to sigclip
#   smaxiter:   Maximum number of srejalg iterations.  Default to 10.
#   lorej:      If srejalg="sigclip", then the number of standard deviations
#               below skyval to clip (default to 3.0).
#               If srejalg="pclip", then fraction of low pixels to clip
#               (default to 0.05).
#   hirej:      If srejalg="sigclip", then the number of standard deviations
#               above skyval to clip (default to 3.0).
#               If srejalg="pclip", then fraction of high pixels to clip
#               (default to 0.05).
#   quick:      Set to use quick_photfrac (much faster)
#   exact:      Set to use exact_photfrac (slower)
#
# OUTPUTS:
#   skyval:     Sky value in counts per pixel.
#   skyrms:     RMS of sky pixels about skyval, again in counts per pixel.
#               This assumes that each pixel is independent.  The RMS
#               is computed for only the values that remain after all the
#               rejection iterations.
#   skyerr:     The error in skyval assuming that each pixel is independent.
#               This is skyrms divided by the square root of the number of
#               pixels.
#
# PROCEDURES CALLED:
#   photfrac
#
# REVISION HISTORY:
#   28-Nov-1996  Written by D. Schlegel, Durham.
#   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
#                of FITS headers; make IDL 5 complient (DJS).
#-
#------------------------------------------------------------------------------
# INTERNAL SUPPORT PROCEDURES:
#
# photsky_compute
# photsky_reject
#------------------------------------------------------------------------------
def photsky_compute(image, salg):

   _expr = salg
   if _expr == 'mean':   
      skyval = np.sum(image)/np.size(image)
      
   elif _expr == 'median':   
      skyval = np.median(image)
      
   elif _expr == 'mode':   
      skyval = 0.0
      # NOT IMPLEMENTED!!!???
      # Can use a KDE to get the mode
      
   elif _expr == None:   
      skyval = 0.0
      
   else:   
      print('Error - Invalid sky algorithm')
      return 0.0
   
   
   return skyval

#------------------------------------------------------------------------------
def photsky_reject(skyval, image, srejalg, lorej, hirej):
   
   _expr = srejalg
   if _expr == 'sigclip':   
      diff = image - skyval
      sigval = np.sqrt(np.sum(diff * diff) / np.size(diff))
      pix = np.nonzero(
               (image > (skyval - lorej*sigval)) &
               (image < (skyval + hirej*sigval)))
      if (np.size(pix)==0):
         newimg = image
      else:   
         newimg = image[pix]
      
   elif _expr == 'pclip':   
      ndata = np.size(image)
      newimg = np.sort(image)
      indx1 = long(lorej*ndata)
      indx2 = ndata - 1 - long(lorej*ndata)
      if (indx2 >= indx1):   
         newimg = newimg[indx1:(indx2)+1]
      
   else:   
      newimg = image
   
   
   return newimg

#------------------------------------------------------------------------------
def photsky(xcen, ycen, skyrad, image, salg=None, srejalg=None, smaxiter=None, lorej=None, hirej=None, quick=None, exact=None):
   """Syntax - result = photsky( xcen, ycen, skyrad, image,  [ salg=, srejalg=, smaxiter=, lorej=, hirej=, skyrms=, skyerr=, /quick, /exact ])"""
#skyval, skyrms, and skyerr are the outputs
   outputs = {'skyval':None, 'skyrms':None, 'skyerr':None}
   
   if (salg is None):   
         salg = 'mean'
   
   if (((salg != 'mean') & (salg !='median') & (salg != 'mode'))|(skyrad is None)):
      outputs['skyval']=None
      outputs['skyrms']=0.0
      outputs['skyerr']=0.0
      return outputs
   
   if ((srejalg is not None) == 0):   
      srejalg = 'sigclip'
   if ((smaxiter is not None) == 0):   
      smaxiter = 10
   
   _expr = srejalg
   if _expr == 'sigclip':   
      if ((lorej is not None) == 0):   
         lorej = 3.0
      if ((hirej is not None) == 0):   
         hirej = 3.0
   elif _expr == pclip:   
      if ((lorej is not None) == 0):   
         lorej = 0.05
      if ((hirej is not None) == 0):   
         hirej = 0.05
   else:   
      if ((lorej is not None) == 0):   
         lorej = 0.0
      if ((hirej is not None) == 0):   
         hirej = 0.0
   
   
   # Find sky pixels and contribution from each
   xdimen = length(image[0,:])
   ydimen = length(image[:,0])
   if (quick is not None):   
   #THIS NEEDS TO BE FIXED
      quick_photfrac(xcen, ycen, skyrad, xdimen=xdimen, ydimen=ydimen, pixnum=pixnum, fracs=fracs, ragged=True)
   else:   
      if (exact is not None): 
         #THIS NEEDS TO BE FIXED  
         exact_photfrac(xcen, ycen, skyrad, xdimen=xdimen, ydimen=ydimen, pixnum=pixnum, fracs=fracs)
      else:   
         output_fracs=photfrac(xcen, ycen, skyrad, xdimen=xdimen, ydimen=ydimen)
         pixnum     = output_fracs['pixnum']
         xpixnum    = output_fracs['xpixnum']
         ypixnum    = output_fracs['ypixnum'] 
         fracs      = output_fracs['fracs']
   
   if (length(pixnum) !=0):   
      
      # Trim to only those pixels more than half filled, unless that rejects
      # all pixels
      filled = np.nonzero(fracs > 0.5)
      count = np.size(filled) 
      if (count != 0):   
         pixnum  = pixnum[filled]
         xpixnum = xpixnum[filled]
         ypixnum = ypixnum[filled]
         fracs   = fracs[filled]
      
      # Iterate until one of the following conditions is met:
      #   (1) the maximum number of iterations is reached
      #   (2) no points are retained after the rejection algorithm
      #   (3) the same sky value is returned after two iterations
      subimg = image[xpixnum,ypixnum]
      iiter = 0
      dsky = 1.0 # Set NE 0 to prevent stopping at first iteration
      while ( (iiter <=smaxiter) & 
              (np.size(subimg) != 0) &
              (dsky != 0.0)):
         if (iiter > 0):   
            dsky = skyval
         skyval = photsky_compute(subimg, salg)
         skydiff = subimg - skyval
         skyrms = np.sqrt(np.sum((skydiff * skydiff)/np.size(skydiff)))
         skyerr = skyrms/np.sqrt( np.size(subimg))
         #ipdb.set_trace()   
         subimg = photsky_reject(skyval, image[xpixnum, ypixnum], srejalg, lorej, hirej)
         
         if (iiter > 0):   
            dsky = dsky - skyval
         iiter = iiter + 1
   else:   
      
      skyval = 0.0  # No sky pixels exist in the given annulus
      skyrms = 0.0
      skyerr = 0.0
      
   outputs['skyval']=skyval
   outputs['skyrms']=skyrms
   outputs['skyerr']=skyerr
   return outputs

#------------------------------------------------------------------------------
# NAME:
#   photcen
#
# PURPOSE:
#   Recenter an object position.
#
# CALLING SEQUENCE:
#   photcen, xcen, ycen, image, $
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
# photcen_gfunc1
# photcen_gfit1()
#------------------------------------------------------------------------------
# Function of the form
#   F(x) = a * exp(-e^2) + c
# where e=(x-x0)/sigma, and the fitting parameters avec=[a,x0,c,sigma].
# This function can be passed to curvefit().

def photcen_gfunc1(xvec, avec, fval, pderiv):
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
def photcen_gfit1(xvec, yvec, wvec, fwhm, xcen, fixfw):

   n_params = 6
   
   global fwhcom
   
   # If too few data points, then return with no change to xcen
   if (array(xvec, copy=0).nelements() < 4):   
      print('Error - Too few data points for gaussian fit')
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
   
   yfit = curvefit(xvec, yvec, wvec, avec, function_name='photcen_gfunc1')
   
   newx = avec[1]
   return newx


#------------------------------------------------------------------------------
def photcen(xcen, ycen, image, calg=None, cbox=None, cmaxiter=None, cmaxshift=None, fwhm=None, fixfw=None, ceps=None, qmaxshift=None):
   """Syntax - result = photcen( xcen, ycen, image, $'
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
      print('Dimensions of NX and NY do not agree')
      pass
   
   if (nx >1):   
      #outputs['xcen']=np.zeros(length(xcen))
      #outputs['ycen']=np.zeros(length(ycen))
      #outputs['qmaxshift']=np.zeros(length(qmaxshift))

      for i in np.arange(0, (nx - 1)+(1)):
         xtmp = xcen[i]
         ytmp = ycen[i]
         outs=photcen(xtmp, ytmp, image, calg=calg,
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
         print('Error - No pixels in X range')
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
         print('Error - No pixels in Y range')
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
         newi = photcen_gfit1(xpixnum[0,:], ptemp, fracx, fwhvec[0], xcen, fixfw)
         xcen = newi
         
         # Fit Y center
         ptemp = total(subimg, 1)
         fracy = yb - ya
         newj = photcen_gfit1(transpose(ypixnum[:,0]), ptemp, fracy, fwhvec[1], ycen, fixfw)
         ycen = newj
         
      elif _expr == 'gauss2':   
         # THIS NEEDS TO BE IMPLEMENTED
         # Two-dimensional gaussian fit, simultaneously fitting in X and Y
         
         gres = gauss2dfit(subimg, gcoeff, tilt=True)
         xcen = gcoeff[4] + istart
         ycen = gcoeff[5] + jstart
         
      elif _expr == None:   
         asdf = -1
         
      else:   
         print('Error - Invalid centering algorithm')
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
# NAME:
#   photfrac
#
# PURPOSE:
#   Create a list of pixel numbers and their fractional contribution to
#   an annular region.
#
# CALLING SEQUENCE:
#   photfrac, xcen, ycen, Rvec, xdimen=, ydimen=, $
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
#   ceil()
#   floor()
#
# REVISION HISTORY:
#   Written D. Schlegel, 27 November 1996, Durham
#   Bug identified - 2 Nov 2000, D. Finkbeiner
#-
#------------------------------------------------------------------------------
# INTERNAL SUPPORT PROCEDURES:
#
# photfrac_intcirc
#------------------------------------------------------------------------------
# Return the integral of the circle with radius Radius within the pixel
# defined by x=[xA,xB] and y=[yA,yB].  Add to this the area beneath the
# pixel, i.e. the value (xB-xA)*yA.
#
# This function assumes that xB > xA >= 0 and yB > yA >= 0, and that
# the circle passes through the pixel:
#   xA^2 + yA^2 <= Radius^2 <= xB^2 + yB^2

def photfrac_intcirc(xa, xb, ya, yb, radius):
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
      gg[pix2] = np.maximum(xanorm[pix2], np.sqrt(1.0 - ybnorm[pix2] * \
                                                  ybnorm[pix2]))
   
   hh = np.minimum(xbnorm, np.sqrt(1.0 - yanorm * yanorm))
   
   result = radius * radius * ((gg - xanorm) * ybnorm + (xbnorm - hh) * yanorm + 0.5 * (hh * np.sqrt(1.0 - hh * hh) + np.arcsin(hh) - gg * np.sqrt(1.0 - gg * gg) - np.arcsin(gg)))
   return result

#------------------------------------------------------------------------------
def photfrac(xcen, ycen, rvec, xdimen=None, ydimen=None, xpixnum=None, ypixnum=None, pixnum=None, fracs=None, fillfrac=None):
   """Syntax - photfrac, xcen, ycen, Rvec, $
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
      print('Error - No pixels in X range')
      pass
   
   jstart = np.maximum(0, 
                np.floor(ycen + 0.5 - radius2).astype('int'))
   jend = np.minimum((ydimen - 1), 
                np.ceil(ycen - 0.5 + radius2).astype('int'))
   jlen = jend - jstart + 1
   if ((jstart > ydimen) |
       (jend < 0) |
       (jlen < 1)):
      print('Error - No pixels in Y range')
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
      integral[pix2] = photfrac_intcirc(xap[iindx[pix2]], xbp[iindx[pix2]], yap[jindx[pix2]], ybp[jindx[pix2]], radius2)
   
   pix1 = np.nonzero((sqradiusa >= sqradius1) & qpix0==True)
   count= np.size(pix1)
   if (count != 0):   
      integral[pix1] = integral[pix1] - (xbp[iindx[pix1]] - xap[iindx[pix1]]) * yap[jindx[pix1]]
   
   pix2 = np.nonzero((sqradiusa < sqradius1) & qpix0==True)
   count = np.size(pix2)
   if (count != 0):   
      integral[pix2] = integral[pix2] - photfrac_intcirc(xap[iindx[pix2]], xbp[iindx[pix2]], yap[jindx[pix2]], ybp[jindx[pix2]], radius1)
   
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
