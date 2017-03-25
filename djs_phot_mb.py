from __future__ import division
from __future__ import print_function
import djs_photsky_mb as djsphotsky
import djs_photcen_mb as djsphotcen
import djs_photfrac_mb as djsphotfrac
from length import length
import numpy as np
import pdb

#------------------------------------------------------------------------------
#+
# NAME:
#   djs_phot
#
# PURPOSE:
#   Driver for aperture photometry with the option of re-centering and
#   sky-subtraction.
#
# CALLING SEQUENCE:
#   flux = djs_phot( xcen, ycen, objrad, skyrad, image, [invvar, $
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
#                 iweight = intensity-weighted center, computed independently
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
#   For similar reasons, if /exact is set, djs_phot will not try to
#     recentroid your object.
# PROCEDURES CALLED:
#   djs_photcen
#   djs_photfrac
#   djs_photsky()
#
# REVISION HISTORY:
#   28-Nov-1996  Written by D. Schlegel, Durham.
#   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
#                of FITS headers; make IDL 5 compliant (DJS).
#   02-Nov-2000  objrad, skyrad recast as floats (D. Finkbeiner)
#                  If they are ints, 1% errors may arise.
#-
#------------------------------------------------------------------------------
def djs_phot(xcen, ycen, objrad, skyrad, image, invvar=None, calg=None, cbox=None, cmaxiter=None, cmaxshift=None, fwhm=None, fixfw=None, ceps=None, salg=None, srejalg=None, smaxiter=None, lorej=None, hirej=None, flerr=None, skyval=None, skyrms=None, skyerr=None, peakval=None, quick=None, exact=None):
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
      message('XCEN and YCEN must have same number of elements')
   
   if ((exact is not None)): 
      #this still needs to be fixed  
      ini = where(ravel(bitwise_or(array(xcen, copy=0).astype(Int32) != xcen, array(ycen, copy=0).astype(Int32) != ycen)))[0]

      if (nni > 0):   
         message('xcen and ycen MUST be integers if /exact keyword is set')
   
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
         
         recenter=djsphotcen.djs_photcen(xcen1, ycen1, image, calg=calg, cbox=cbox, cmaxiter=cmaxiter, cmaxshift=cmaxshift)
         xcen1 = recenter['xcen']
         ycen1 = recenter['ycen']
      else:   
         print('Not recentering -- /exact requires correct input center')
         xcen1 = xcen[iobj]
         ycen1 = ycen[iobj]
      
      # Find the sky value
      # Add 0.0 to the skyrad to convert to floating-point
      
      skyvals= djsphotsky.djs_photsky(xcen1, ycen1, skyrad + 0.0, image, salg=salg, srejalg=srejalg, smaxiter=smaxiter, lorej=lorej, hirej=hirej, quick=quick, exact=exact)
      skyval[iobj] = skyvals['skyval']
      skyrms[iobj] = skyvals['skyrms']
      skyerr[iobj] = skyvals['skyerr']

      # Find total counts in object aperture(s)
      for irad in np.arange(0, (nrad - 1)+(1)):
         onerad = objrad[irad] + 0.0 # Convert to floating-point
         if (onerad > 0):   
            if (quick is not None):   
               quick_photfrac(xcen1, ycen1, onerad, xdimen=xdimen, ydimen=ydimen, pixnum=pixnum, fracs=fracs)
            else:   
               if (exact is not None):   
                  exact_photfrac(xcen1, ycen1, onerad, xdimen=xdimen, ydimen=ydimen, pixnum=pixnum, fracs=fracs)
               else:  
                  output_object=djsphotfrac.djs_photfrac(xcen1, ycen1, onerad, xdimen=xdimen, ydimen=ydimen)
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
            flux[irad,iobj] = np.sum(image[xpixnum,ypixnum] * fracs) - (area * skyval[iobj])
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
