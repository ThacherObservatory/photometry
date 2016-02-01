from length import length
import pdb
import numpy as np
import djs_photfrac_mb as djsphotfrac
#------------------------------------------------------------------------------
#+
# NAME:
#   djs_photsky
#
# PURPOSE:
#   Compute the sky value within an annulus.
#
#   At present, fractional pixels are not treated properly; this is because
#   we need a sort routine that carries index numbers, such as NR indexx().
#
# CALLING SEQUENCE:
#   skyval = djs_photsky( xcen, ycen, skyrad, image, $
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
#   djs_photfrac
#
# REVISION HISTORY:
#   28-Nov-1996  Written by D. Schlegel, Durham.
#   01-Jun-2000  Major revisions: change XYCEN to XCEN,YCEN; remove use
#                of FITS headers; make IDL 5 complient (DJS).
#-
#------------------------------------------------------------------------------
# INTERNAL SUPPORT PROCEDURES:
#
# djs_photsky_compute
# djs_photsky_reject
#------------------------------------------------------------------------------
def djs_photsky_compute(image, salg):

   n_params = 2
   
   _expr = salg
   if _expr == 'mean':   
      skyval = np.sum(image)/np.size(image)
      
   elif _expr == 'median':   
      skyval = np.median(image)
      
   elif _expr == 'mode':   
      skyval = 0.0
      # NOT IMPLEMENTED!!!???
      
   elif _expr == None:   
      skyval = 0.0
      
   else:   
      print 'Error - Invalid sky algorithm'
      return 0.0
   
   
   return skyval

#------------------------------------------------------------------------------
def djs_photsky_reject(skyval, image, srejalg, lorej, hirej):
   
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
def djs_photsky(xcen, ycen, skyrad, image, salg=None, srejalg=None, smaxiter=None, lorej=None, hirej=None, quick=None, exact=None):
   """Syntax - result = djs_photsky( xcen, ycen, skyrad, image,  [ salg=, srejalg=, smaxiter=, lorej=, hirej=, skyrms=, skyerr=, /quick, /exact ])"""
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
         output_fracs=djsphotfrac.djs_photfrac(xcen, ycen, skyrad, xdimen=xdimen, ydimen=ydimen)
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
         skyval = djs_photsky_compute(subimg, salg)
         skydiff = subimg - skyval
         skyrms = np.sqrt(np.sum((skydiff * skydiff)/np.size(skydiff)))
         skyerr = skyrms/np.sqrt( np.size(subimg))
         #ipdb.set_trace()   
         subimg = djs_photsky_reject(skyval, image[xpixnum, ypixnum], srejalg, lorej, hirej)
         
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
