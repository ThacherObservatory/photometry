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

# Yao wins!!!!
data2013 = pd.read_table('/Users/sara/python/photometry/Landolt_Standards_2013.txt',sep='|', header=62)
data2009 = pd.read_table('/Users/sara/python/photometry/Landolt_Standards_2009.txt',sep='|', header=66)

thob = ephem.Observer()
thob.long = ephem.degrees("-119.1773417")
thob.lat = ephem.degrees("34.467028")
thob.elevation = 504.4

def filter_data(data,datetime, mag, alt, az, VI):
    """
    Parameters:
    1. data base of standard stars
    2. local datetime entered in the form of 'YY MMM DD hh:mm:ss'
    3. magnitude, altitude, azimuth and V-I restritions, list, [min_value, max_value]

    Returns:
    Pandas DataFrame for the information of the stars given the ristrictions
    """
    thob.date = gmtime(mktime(strptime(datetime,'%y %b %d %H:%M:%S')))[:6]
    namel = []
    azl = []
    altl = []
    magl = []
    VIl = []
    for i in range(len(data)-1):
        dataline = data['Name'][i]+','+'f'+','+data['RAJ2000'][i]+','+data['DEJ2000'][i]+','+str(data['Vmag'][i])
        body = ephem.readdb(dataline)
        body.compute(thob)
        if alt[0] < body.alt < alt[1] and \
            az[0] < body.az < az[1] and \
            mag[0]< body.mag < mag[1] and \
            VI[0]< data['V-I'][i] <VI[1]:
            namel.append(body.name)
            azl.append(body.az)
            altl.append(body.alt)
            magl.append(body.mag)
            VIl.append(data['V-I'])
    filtered_data = pd.DataFrame.from_dict(dict([('Name',namel), ('Az',azl), ('Alt',altl),('Mag',magl),('V-I',VIl)]))
    return filtered_data

filter_data(data2009,"17 Feb 05 16:22:56",[10,15],[0.5,1],[3,6],[0.2,1])
