import pandas as pd
import matplotlib.pyplot as pl
import datetime
import numpy as np
import ephem
from scipy import interpolate
from time import strptime, gmtime, mktime
from astropy.coordinates import Angle
import astropy.units as u

class observatory(ephem.observe):

    def __init__(self,long,lat,elev,dir):
        super(Observatory,self).__init__()
        self.long = -119.1773417
        self.lat = 34.467028
        self.elev = 504.4
        self.data2013 = pd.read_table(dir+'Landolt_Standards_2013.txt',sep='|', header=62)
        self.data2009 = pd.read_table(dir+'Landolt_Standards_2009.txt',sep='|', header=66)

    def findStandardRADec(RA,Dec,dist=None,**kwargs):
        """
        *Args:
        1. ra_dec: specify ra and dec:
            if dist=None, specify the range of ra and dec, list, in dms,[min_ra,max_ra,min_dec,max_dec]
            else, one values of ra and dec, list, in dms, [ra,dec]
        2. dist specified in dms string format, default is None

        **Kwargs:
        1. mag: list of 2 ints, in mags
        2. V-I: list of 2 ints
        3. data: data2009 as default, data 2013 as an option

        Returns:
        Pandas data frame opf standard stars meeting the specified criteria
        """

    def findStandardAltAz(Alt,Az,):
        """
        *Args:
        Local datetime entered in the form of 'YY MMM DD hh:mm:ss'
        Altitude, Azimuth restritions, list, [min_value, max_value]

        **Kwargs:
        Magnitude and V-I restritions, list, [min_value, max_value]
        data: data2009 as default, data 2013 as an option

        Returns:
        Pandas data frame opf standard stars meeting the specified criteria
        """
        
