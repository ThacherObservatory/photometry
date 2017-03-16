# -*- coding: utf-8 -*-
"""
Created on Sun Feb 05 10:58:47 2017

@author: yao
"""

import pandas as pd
import matplotlib.pyplot as pl
import datetime
import numpy as np
import ephem
from scipy import interpolate
from time import strptime, gmtime, mktime
from astropy.coordinates import Angle
import astropy.units as u

data2013 = pd.read_table('/Users/yaosarayin/python/photometry/Landolt_Standards_2013.txt',sep='|', header=62)
data2009 = pd.read_table('/Users/yaosarayin/python/photometry/Landolt_Standards_2009.txt',sep='|', header=66)

thob = ephem.Observer()
thob.long = ephem.degrees("-119.1773417")
thob.lat = ephem.degrees("34.467028")
thob.elevation = 504.4

def radian_to_dms(radian):
    totalSeconds = radian * 360. * 60. * 60. / (2. * np.pi)
    seconds = totalSeconds % 60
    minutes = (totalSeconds / 60.) % 60
    degrees = totalSeconds / (60. * 60.)
    return '%.2f:%.2f:%.2f' % (degrees,minutes,seconds)

def dms_to_radian(dms,sep=' '):
    """
    "dd sep mm sep ss"
    """
    degrees = np.sum(np.fromstring(dms, sep=sep)*[1.0, 1/60.0, 1/3600.0])
    return np.radians(degrees)

def filter_data(data,datetime, mag=False, alt=False, az=False, VI=False, ra_dec=False, dist=False):
    """
    Parameters:
    1. data base of standard stars
    2. local datetime entered in the form of 'YY MMM DD hh:mm:ss'
    3. magnitude, altitude, azimuth and V-I restritions, list, [min_value, max_value]
    4. ra_dec: specify ra and dec:
        if dist=False, specify the range of ra and dec, list, in dms,[min_ra,max_ra,min_dec,max_dec]
        if dist==True, one values of ra and dec, list, in dms, [ra,dec]
    5. dist: if True, specify distance from one value od ra and dec, list , in dms, [ra_range, dec_range]
    Returns:
    Pandas DataFrame for the information of the stars given the ristrictions
    """
    thob.date = gmtime(mktime(strptime(datetime,'%y %b %d %H:%M:%S')))[:6]
    namel = []
    azl = []
    altl = []
    magl = []
    VIl = []
    ral = []
    decl = []
    selected=[]
        #apply function to series pd.Series.apply or pd.DataFrame.apply data['RAJ2000'] = data['RAJ2000'].apply(dms_to_radian)
    if ra_dec and dist:
        try:
            ra,ra_dist=Angle([ra_dec[0],dist[0]],unit=u.hour); dec,dec_dist=Angle([ra_dec[1],dist[1]],unit=u.degree)
            for i in range(len(data)-1):
                if ((ra-ra_dist)<Angle(data['RAJ2000'][i],unit=u.hour)<(ra+ra_dist)) and \
                        ((dec-dec_dist)<Angle(data['DEJ2000'][i],unit=u.degree)<(dec+dec_dist)):selected.append(i)
        except:
            print 'ValueError: ra_dec must be list of 2 values [ra,dec] if dist==True'
    elif ra_dec and not dist:
        try:
            min_ra,max_ra = Angle([min_ra,max_ra],unit=u.hour)
            min_dec,max_dec = Angle([min_dec,max_dec],unit=u.degree)
            for i in range(len(data-1)):
                if data[min_ra<data['RAJ2000'][i]<max_ra and min_dec<data['DEJ2000'][i]<max_dec]: selected.append(i)
        except:
            print 'ValueError: ra_dec must be list of 4 values [min_ra,max_ra,min_dec,max_dec] if dist==False'
    else: selected = range(len(data)-1)
    data = data.loc[selected,:]
    for i in selected:
        dataline = data['SimbadName'][i]+','+'f'+','+data['RAJ2000'][i]+','+data['DEJ2000'][i]+','+str(data['Vmag'][i])
        body = ephem.readdb(dataline)
        body.compute(thob)
        if mag and alt and az and VI:
            if np.radians(alt[0]) < body.alt < np.radians(alt[1]) and \
                np.radians(az[0]) < body.az < np.radians(az[1]) and \
                mag[0]< body.mag < mag[1] and \
                VI[0]< data['V-I'][i] <VI[1]:
                namel.append(body.name);azl.append(body.az);altl.append(body.alt)
                magl.append(body.mag);VIl.append(data['V-I'])
                ral.append(radian_to_dms(body.ra));decl.append(radian_to_dms(body.dec))
        else:
            namel.append(body.name);azl.append(body.az);altl.append(body.alt)
            magl.append(body.mag);VIl.append(data['V-I'])
            ral.append(radian_to_dms(body.ra));decl.append(radian_to_dms(body.dec))
    filtered_data = pd.DataFrame.from_dict(dict([('Name',namel), ('Az',azl), ('Alt',altl),('Ra',ral),('Dec',decl),('Mag',magl),('V-I',VIl)]))
    return filtered_data

#filtered_data = filter_data(data2009,"17 Mar 16 20:00:00",ra_dec=['13:26:14','-07:28:09'],dist=['10:05:00','22:00:00'])
#filtered_data = filter_data(data2009,"17 Mar 16 20:00:00",mag=[8,12],alt=[20,80],az=[0,360],VI=[0.5,1])
