# -*- coding: utf-8 -*-
"""
Created on Sun Feb 05 10:58:47 2017

@author: yao
"""
"""
keys:
Name
Vmag
e_vmag
b-v
u-b
v-r
r-i
v-i
ra de J2000
pmra and pmde
SimbadName
"""
import pandas as pd
import matplotlib.pyplot as pl
import numpy as np
import ephem
from scipy import interpolate
from time import strptime, gmtime, mktime
from astropy.coordinates import Angle
import astropy.units as u
import pdb

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

def get_Landolt():
    data2013 = pd.read_table('Landolt_Standards_2013.txt',sep='|', header=62)
    data2009 = pd.read_table('Landolt_Standards_2009.txt',sep='|', header=66)
    data2009 = data2009[['Name','Vmag','e_Vmag','B-V','U-B','V-R','R-I','V-I','RAJ2000','DEJ2000','pmRA','pmDE','SimbadName']]
    data2013 = data2013[['Name','Vmag','e_Vmag','B-V','U-B','V-R','R-I','V-I','RAJ2000','DEJ2000','pmRA','pmDE','SimbadName']]
    return pd.concat([data2009,data2013],ignore_index=True)


def filter_data(data,datetime, mag=False, alt=False, az=False, VI=False, ra_dec=False, dist=False):
    """
    Parameters:
    1. data base of standard stars
    2. local datetime entered in the form of 'YY MMM DD hh:mm:ss'
    3. magnitude, altitude, azimuth and V-I restritions, list, [min_value, max_value]
    4. ra_dec: specify ra and dec AS TUPLE:
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
            ra,ra_dist=Angle([ra_dec[0],dist[0]],unit=u.hour)
            dec,dec_dist=Angle([ra_dec[1],dist[1]],unit=u.degree)
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
    data = data.loc[selected,:]
    for i in selected:
        dataline = data['SimbadName'][i]+','+'f'+','+data['RAJ2000'][i]+','+data['DEJ2000'][i]+','+str(data['Vmag'][i])
        body = ephem.readdb(dataline)
        body.compute(thob)

        if mag or alt or az or VI:
            namel.append(body.name)
            ral.append(radian_to_dms(body.ra))
            decl.append(radian_to_dms(body.dec))
            if alt and az:
                if np.radians(alt[0]) < body.alt < np.radians(alt[1]) and \
                    np.radians(az[0]) < body.az < np.radians(az[1]):
                    azl.append(body.az)
                    altl.append(body.alt)
                else:
                    azl.append('NaN')
                    altl.append('NaN')
            elif mag:
                if mag[0]< body.mag < mag[1]:
                    magl.append(body.mag)
                else:
                    magl.append('NaN')
            elif VI:
                if VI[0]< data['V-I'][i] <VI[1]:
                    VIl.append(data['V-I'][i])
                else:
                    VIl.append('NaN')

        else:
            namel.append(body.name)
            azl.append(body.az)
            altl.append(body.alt)
            magl.append(body.mag)
            VIl.append(data['V-I'])
            ral.append(radian_to_dms(body.ra))
            decl.append(radian_to_dms(body.dec))

    filtered_data = pd.DataFrame.from_dict(dict([('Name',namel), ('Az',azl), ('Alt',altl),('Ra',ral),('Dec',decl),('Mag',magl),('V-I',VIl)]))

    return filtered_data

#filtered_data = filter_data(data2009,"17 Feb 15 20:00:00",ra_dec=['13:26:14','-07:28:09'],dist=['10:05:00','22:00:00'])

# reorganizing the previous function for better ease of use
def filters(data,ra,dec,datetime,dist=10,mag=False,alt_az=False,VI=False):
    """
    Parameters:
    data: pandas dataframe of standard stars. eg  get_Landolt()
    ra: string of the ra. eg '13:26:14'
    dec: string of the dec. eg '-07:28:09'
    datetime: local datetime in form of "YY Mmm DD hh:mm:ss". eg "17 Feb 15 20:00:00"
    dist: angular distance from ra and dec. is preset to 10 but can change based on needs
    mag: list of lower lim of mag, then upper lim. eg [9,13]
    alt_az: list of alt az lims in the form of ra and dec, [alt_lower,alt_upper,az_lower,az_upper]
    VI: same form as mag
    """

    ra_lower = Angle(ra,unit=u.hour) - Angle(dist,unit=u.hour)
    ra_upper = Angle(ra,unit=u.hour) + Angle(dist,unit=u.hour)

    dec_lower = Angle(dec,unit=u.degree) - Angle(dist,unit=u.degree)
    dec_upper = Angle(dec,unit=u.degree) + Angle(dist,unit=u.degree)

    thob.date = gmtime(mktime(strptime(datetime,'%y %b %d %H:%M:%S')))[:6]

    for i in range(len(data)-1):
        """
        dataline = data['SimbadName'][i]+','+'f'+','+data['RAJ2000'][i]+','+data['DEJ2000'][i]+','+str(data['Vmag'][i])
        body = ephem.readdb(dataline)
        body.compute(thob)
        """

        if not (Angle(data['RAJ2000'][i],unit=u.hour) > ra_lower and Angle(data['RAJ2000'][i],unit=u.hour) < ra_upper):
            data = data.drop(data.index[[i]])
        if not (Angle(data['DEJ2000'][i],unit=u.degree) > dec_lower and Angle(data['DEJ2000'][i],unit=u.degree) < dec_upper):
            data = data.drop(data.index[[i]])
        if mag:
            if not (data['Vmag'][i] > mag[0] and data['Vmag'][i] < mag[1]):
                data = data.drop(data.index[[i]])
        if VI:
            if not (data['V-I'][i] > VI[0] and data['V-I'][i] < VI[1]):
                data = data.drop(data.index[[i]])
        """
        elif not (Angle(body.alt,unit=u.hour) > Angle(alt_az[0],unit=u.hour) and Angle(body.alt,unit=u.hour) <Angle(alt_az[1],unit=u.hour)):
            data = data.drop(data.index[[i]])
        elif not (Angle(body.az,unit=u.degree) > Angle(alt_az[2],unit=u.degree) and Angle(body.az,unit=u.degree) <Angle(alt_az[3],unit=u.degree)):
            data = data.drop(data.index[[i]])
        """
    return data



