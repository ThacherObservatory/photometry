import pandas as pd
import matplotlib.pyplot as pl
import datetime
import numpy as np
import ephem
from scipy import interpolate

# Yao wins!!!!
data2013 = pd.read_table('Landolt_Standards_2013.txt',sep='|', header=62)
data2009 = pd.read_table('Landolt_Standards_2009.txt',sep='|', header=66)

thob = ephem.Observer()
thob.long = ephem.degrees("-119.1773417")
thob.lat = ephem.degrees("34.467028")
thob.elevation = 504.4

def read_data(data):
    for i in len(data)
        dataline = data['Name'][i]+'f'+data['RAJ2000'][i]+data['DEJ2000'][i]+data['Vmag'][i]
    body = ephem.readdb(dataline)

def body_info(date, starttime, endtime, position, ra, dec, plot=False):
    """
    This function returns the altitude and azimuth of stars of certain magnitude within a time period
    """
    thob.date = date+time #both should be strings in format '1984/5/30 16:22:56'
    body =[]
    for i in len(data):
        dataline = data['Name'][i]+'f'+data['RAJ2000'][i]+data['DEJ2000'][i]+data['Vmag'][i]
        body = ephem.readdb(dataline)
        body = ephem.FixedBody(ra=ra,dec=dec)
        body.compute(thob)
        risetime = body.next_rising(start='date')
        settime = body.next_setting(start='date')
        if risetime > starttime and settime < endtime and body.mag > mag:
            body.append()
    az = body.az
    alt = body.alt
    return risetime, settime, az, alt
