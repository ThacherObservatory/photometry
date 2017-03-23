import pandas as pd
import matplotlib.pyplot as pl
import datetime
import numpy as np
import ephem
from scipy import interpolate
from time import strptime, gmtime, mktime
from astropy.coordinates import Angle
import astropy.units as u

class Observatory(ephem.Observer):

    def __init__(self, long=-119.1773417, lat=34.467028, elev=504.4,
                 dir="./"):
        """
        Observatory object that utilizes pyephem's Observer class to create a
        custom observatory given a longitude, latitude, and elevation

        long: float (degrees)
            the longitude of your observatory
        lat: float (degrees)
            the latitude of your observatory
        elev: float (meters)
            the elevation from sea level of your observatory
        dir: str
            the path which contains any files you may read in
            * defaulted to your present working directory
        """
        super(Observatory, self).__init__()
        self.long = long
        self.lat = lat
        self.elev = elev
        self.data2013 = pd.read_table(dir+'Landolt_Standards_2013.txt',sep='|', header=62)
        self.data2009 = pd.read_table(dir+'Landolt_Standards_2009.txt',sep='|', header=66)

    def findStandardRADec(RA,Dec,dist=None,**kwargs):
        """
        Parameters:
        1. RA/Dec: specify right ascension and declination:
            if dist=None, specify the range of ra and dec, list, in dms,[min_ra,max_ra],[min_dec,max_dec]
            else, one values of ra and dec, list, in dms, [ra,dec]
        2. dist specified in degrees-minutes-seconds string format
            * default is None

        **Kwargs:
        1. mag: list of 2 ints, in mags
        2. V-I: list of 2 ints
        3. data: data2009 as default, data 2013 as an option

        Returns:
        Pandas data frame opf standard stars meeting the specified criteria
        """
        data = kwargs.get('data', d=self.data2009)
        namel = []; magl = []; VIl = []; ral = []; decl = []; selected=[]
        if dist:
            ra,ra_dist = Angle([RA, dist[0]], unit = u.hour); dec, dec_dist = Angle([Dec, dist[1]], unit = u.degree)
            for i in range(len(data)-1):
                if ((ra - ra_dist) < Angle(data['RAJ2000'][i], unit = u.hour) < (ra + ra_dist)) and \
                        ((dec-dec_dist) < Angle(data['DEJ2000'][i], unit=u.degree) < (dec + dec_dist)): selected.append(i)
        else:
            min_ra,max_ra = Angle([RA[0],RA[1]], unit = u.hour)
            min_dec,max_dec = Angle([Dec[0], Dec[1]], unit = u.degree)
            for i in range(len(data-1)):
                if data[min_ra < data['RAJ2000'][i] < max_ra and min_dec < data['DEJ2000'][i] < max_dec]: selected.append(i)
        return

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
        return
