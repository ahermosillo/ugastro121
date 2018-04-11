import numpy as np
import matplotlib.pyplot as plt
import astropy as ap
from astropy.time import Time 
from astropy import coordinates
from astropy import units as u
import re #regex, if needed
import time
from threading import Timer
from datetime import datetime
from pytz import timezone
import pytz
#import pandas as pd
# %matplotlib inline
import ugradio
import leuschner

def PST_to_uni(yr, mon, day, hr = 0, minute = 0,sec = 0):
    return int(time.mktime(time.strptime(
        str(yr) +'-'+str(mon)+'-'+str(day)+' '+str(hr)+':'+str(minute)+':'+str(sec), 
        '%Y-%m-%d %H:%M:%S')))
        
uni_to_jul = lambda unix_t: Time(unix_t, format='unix').jd

def get_altaz(ra,dec,jd=None,lat = 37.9183, lon = -122.1067, alt = 304, equinox='J2000'):
    """
    Return the altitude and azimuth of an object whose right ascension 
    and declination are known.
    Parameters
    ----------
    ra : float, right ascension in degrees
    dec: float, declination in degrees
    jd : float, Julian Date, default=now
    lat: float, latitude in degrees
    lon: float, longitude in degrees
    alt: float, altitude in m, 
    equinox : string, equinox of ra/dec coordinates.  default='J2000'
    Returns
    -------
    alt : float, altitude in degrees
    az : float, azimuth in degrees

    """
    if jd: t = ap.time.Time(jd,format='jd')
    else: t = ap.time.Time(time.time(),format='unix')
    l = ap.coordinates.EarthLocation(lat=lat*u.deg,
                        lon=lon*u.deg,height=alt*u.m)
    f = ap.coordinates.AltAz(obstime=t,location=l)
    c = ap.coordinates.SkyCoord(ra, dec, frame='fk5',unit='deg',equinox=equinox)
    altaz = c.transform_to(f)
    return altaz.alt.deg, altaz.az.deg
    
def convert_galactic2000_to_ra_dec(l,b):
    x0 = np.cos(np.radians(b))*np.cos(np.radians(l))
    x1 = np.cos(np.radians(b))*np.sin(np.radians(l))
    x2 = np.sin(np.radians(b))

    R = np.matrix([[-0.054876,-0.873437,-0.483835],[0.494109,-0.444830,0.746982],[-0.867666,-0.198076,0.455984]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(np.transpose(R),x)

    ra =  float(np.degrees(np.arctan2(xp[1],xp[0])))
    dec = float(np.degrees(np.arcsin(xp[2])))

    return ra, dec
    
    
def collect_data(ra,dec,unix,Nspectra,dt,fileName,fitName,noise=False):
    """
    This code assumes that you have already found the conversion between galactic coordinates
    
    ra: array of RA's (ra corresponds to a different galactic coordinate)
    dec: array of DEC's (dec corresponds to a different galactic coordinate)
    unix: start time for what we are observing
    Nsepctra: number of spectra to collect for each point
    dt: time to stay at each point
    noise: bool, default False, if True we collect noise data
    """
    with open('{}'.format(fileName), 'w') as pointFile:
        pointFile.write('{}'.format('agilent'))
            
        alt, az = get_altaz(ra[0],dec[0],jd =uni_to_jul(unix), lat=37.9183, lon=-122.1067, alt =304)
        LeuschTelescope.point(alt,az)
        print(LeuschTelescope.get_pointing())

        if noise:
            ugradio.leusch.LeuschNoise()
            LeuschNoise.on()
            
        ugradio.agilent.SynthClient(host='127.0.0.1')
        pointFile.write('{}'.format(SynthClient.get_frequency()))
        
        #initialize spectrometer thing
        leuschner.Spectrometer('10.0.1.2')
        
        for r,d in zip(ra,dec):
            obsv_time = uni_to_jul(time.time())
            alt,az = get_altaz(ra[0],dec[0], jd=obsv_time, lat=37.9183, lon=-122.1067, alt = 304)
            LeuschTelescope.point(alt,az)
            currentAlt, currentAz = leusch.get_pointing()
            print('alt: {} , az: {}'.format(currentAlt, currentAz))
            Spectrometer.read_spec('{}_{}_r_d.fits'.format(unix,fitName), Nspec, (r,d), 'eq')

ra_120,dec_120 = convert_galactic2000_to_ra_dec(120,0)

def test_collect(ra,dec,unix,Nspec):
    alt,az = get_altaz(ra,dec)
    leuTel = ugradio.leusch.LeuschTelescope()
    spec = leuschner.Spectrometer('10.0.1.2')
    leuTel.point(alt,az)

    currentAlt, currentAz = leuTel.get_pointing()
    print('alt: {} , az: {}'.format(currentAlt, currentAz))

    spec.read_spec('{}_{}_{}.fits'.format(unix,ra,dec), Nspec, (ra,dec), 'eq')


N = 50

test_collect(ra_120,dec_120,time.time(),N)
