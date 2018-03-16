import numpy as np
import astropy as ap
from astropy import units as u
from astropy.time import Time 
from astropy import coordinates

import matplotli.pyplot as plt
import ugradio
import time
import time
from threading import _Timer # Need this for python 2. Just Timer for python 3

#Kevin put these but I'm not using them.. 
import re #regex
from datetime import datetime
from pytz import timezone
import pytz
#Smaller functions for time conversions. All functions take one numerical values

###Seconds and hours
hrs_to_sec = lambda x: x*3600
sec_to_hrs = lambda x: x*1.0/3600

## Class to make a timer of when a function will run. Don't Actually use this 
## fearing that the computer will timeout/crash, but it works lol
class CustomTimer(Timer):
    '''
    Creates a custom Timer from the threading package. Works the same way as Timer
    but has an additional join function
    
    '''
    def __init__(self, interval, function, args=[], kwargs={}):
        self._original_function = function
        super(CustomTimer, self).__init__(
            interval, self._do_execute, args, kwargs)

    def _do_execute(self, *a, **kw):
        self.result = self._original_function(*a, **kw)

    def join(self):
        super(CustomTimer, self).join()
        return self.result

###PST to Unix
def PST_to_uni(yr, mon, day, hr = 0, minute = 0,sec = 0):
    """
    Returns unix time for given PST
    yr, mon, day, hr, minute, sec: (ints or floats) Int that order
    hr, minute, sec = 0 
    """
    return int(time.mktime(time.strptime(
        str(yr) +'-'+str(mon)+'-'+str(day)+' '+str(hr)+':'+str(minute)+':'+str(sec), 
        '%Y-%m-%d %H:%M:%S')))


#Skeleton of original conversions from the web
#time.mktime(time.strptime('2018-03-14 3:0:0', '%Y-%m-%d %H:%M:%S')) #Convert PST to Unix
#datetime.fromtimestamp(time.time()
#).astimezone(timezone('US/Pacific')).strftime('%Y-%m-%d %H:%M:%S %Z%z') #Convert Unix to PST

def findAltAz(ra, dec, jdT = None, precession = False):
    """
    Gives the altitude and azimuth given a certain ra, dec and Julian day.
    If no Julian day is given, will find it at the current JD
    
    Arguments: 
    ra: (float) ra of source 
    dec: (float) dec of source
    jd_array: (float) Julian day
    precession: (bool) If True will get precessed ra and dec 

    Returns:
    two numbers: altitude, azimuth
    """
    if precession: 
        ra, dec = ugradio.coord.precess(ra, dec, jd = jdT)
        al , az = ugradio.coord.get_altaz(ra, dec, jd = jdT)
        return al, az
    else:
        return ugradio.coord.get_altaz(ra, dec, jd = jdT)

def sunAltAz(jdT = None):
    """
    Gives the altitude and azimuth of the sun given Julian day. If no 
    Julian day is given find alt, az at current time

    Arguments: 
    jdT: (float) Julian day

    Returns: 
    Two numbers: altitude, azimuth
    """
        #Sun position in ra, dec
    if jdT: raS, decS = ugradio.coord.sunpos(jd = jdT)
    else: raS, decS = ugradio.coord.sunpos()
    return ugradio.coord.get_altaz(raS, decS, jd = jdT)
    
def moonAltAz(jdT = None):
    """
    Give the altitude and azimuth of the moon at current or given Julian day
    
    Arguments:
    jdT: (float) Julian day
    
    Return:
    Two numbers: altitude, azimuth
    """
    #Moon position in ra, dec
    if jdT: 
        raMoon, decMoon = ugradio.coord.moonpos(jd = jdT)
    else: 
        raMoon, decMoon = ugradio.coord.moonpos()
    return ugradio.coord.get_altaz(raMoon, decMoon, jd = jdT) #Translates ra, dec to al, az

def corrections(al, az):
    '''
    Corrects altitude and azimuth if not in the range of the interferometers.
    Altitude Range [5,175]
    Azimuth Range [90,300]

    al: float, the altitude of a given input
    az: float, the azimuth of a given input

    returns: two floats, al, az, which are the corrected values for our interferometers 
    to point at. Altitude is done manually here by us (the range in which we
    collect the data).
    '''
    assert al > 5 and al < 175, "Altitude not in range!!!"
    if az < 90:
        az = az + 180
        al = abs(180 - al)
    elif az > 300:
        az = az - 180
        al = abs(180 - al)
    return al, az
    
def collect_data(ra,dec,unix,end_t,dt,n,dataFile,pointfileName ):
        '''
        Collects data for a point source whose ra and dec is already known. 

        ra: (float) the right ascension of our galactic coordinate
        dec: (float) the declination of our galactic coordinate
        unix: (float) start time of observing, in seconds (unix)
        end_t: (float) end time of observing, in seconds (unix)
        dt: (float) updating the interferometer new coordinates by every dt seconds
        n: (float) how fast we're collecting data from the multimeter
        dataFile: (str) what we want to name our voltage times series. This will be a npz file
        pointfileName: (str) Name of info txt file. 
        Has Julian Day, measured Alt, measured Az, corrected Alt, corrected Az, 
        w Alt, w Az, e Alt, e Az

        Returns:
        A new file every 15 minutes and a file at the very end
        Files have name 'dataFile_unix' and 'dataFile_lat' and 'pointfileName'

        '''
    with open('{}'.format(pointfileName), 'w') as pointFile:
        pointFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Julian Date','Measured Alt','Measured Az',
                                    'Correct Alt','Correct Az','westAltitude', 'westAzimuth', 'eastAltitude', 'eastAzimuth', ))

        alts, azs = findAltAz(ra,dec,uni_to_jul(unix))
        al, az = corrections(alts,azs)
        ifm.point(al,az)
        print(ifm.get_pointing()) # where we start pointing. before we collect data
        hpm.start_recording(n)
        point_time = time.time()
        
        while unix < end_t:
            unix = time.time()
            jul = ugradio.timing.julian_date(unix)
            alts, azs = findAltAz(ra,dec,jul)
            al, az = corrections(alts,azs)
            ifm.point(al,az)
            get_point = ifm.get_pointing()
            valAltAz = get_point.values()
            pointFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(jul, alts, azs, 
                                    al, az, valAltAz[0][0], valAltAz[0][1], valAltAz[1][0], valAltAz[1][1]))
            print('{}\t{}\t{}\t{}\n'.format(valAltAz[0][0], valAltAz[0][1], valAltAz[1][0], valAltAz[1][1]))

            time.sleep(dt)
            if unix - point_time >= 15*60: 
                point_time = time.time()
                volts, times = hpm.get_recording_data()
                np.savez('{}_{}'.format(dataFile, jul), v = volts, tt = times)

    vend, tend = hpm.end_recording()
    np.savez('{}_last'.format(dataFile), vv = vend, ttt = tend) 
    return volts, times
    
def collect_data_sun(unix,end_t,dt,n,dataFile,pointfileName ):
    '''
    Collects voltage data from the sun using interferometer and multimeter 
    
    unix: (float) start time of observing, unix seconds
    end_t: (float) end time of observing, unix seconds
    dt: (float) updating the interferometer new coordinates by every dt seconds
    n: (float) how fast we're collecting data from the multimeter
    dataFile: (str) what we want to name our voltage times series
    pointfileName: str, Name of file with time and alts, azs measurements from measured,
    correct, west and east interferometer
    '''
    with open('{}'.format(pointfileName), 'w') as pointFile:
        pointFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Unix Time','Measured Alt','Measured Az',
                                    'Correct Alt','Correct Az','westAltitude', 'westAzimuth', 'eastAltitude', 'eastAzimuth', ))

        alts, azs = sunAltAz()
        al, az = corrections(alts,azs)
        ifm.point(al,az)
        print(ifm.get_pointing())
        hpm.start_recording(n)

        point_time = time.time()

        while unix < end_t:
            unix = time.time()
            alts, azs = sunAltAz()
            al, az = corrections(alts,azs)
            ifm.point(al,az)
            get_point = ifm.get_pointing()
            valAltAz = get_point.values()
            pointFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(unix, alts, azs, 
                                    al, az, valAltAz[0][0], valAltAz[0][1], valAltAz[1][0], valAltAz[1][1]))
            print('{}\t{}\t{}\t{}\n'.format(valAltAz[0][0], valAltAz[0][1], valAltAz[1][0], valAltAz[1][1]))

            time.sleep(dt)
            if unix - point_time >= 15*60: 
                point_time = time.time()
                volts, times = hpm.get_recording_data()
                np.savez('{}'+str(unix).format(dataFile), v = volts, tt = times)

    vend, tend = hpm.end_recording()
    np.savez('{}_last'.format(dataFile), vv = vend, ttt = tend) 
    return volts, times
    
def collect_data_moon(unix,end_t,dt,n,dataFile,pointfileName ):
        '''
        Collects voltage data for the moon using interferometers on campbell hall and multimeter functions

        unix: float, start time of observing
        end_t: float, end time of observing
        dt: float, updating the interferometer new coordinates by every dt seconds
        n: float, how fast we're collecting data from the multimeter
        datafile: str, what we want to name our voltage times series
        pointfile: str, what we want to name our coordinates 

        '''
    with open('{}'.format(pointfileName), 'w') as pointFile:
        pointFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('Julian Date','Measured Alt','Measured Az',
                                    'Correct Alt','Correct Az','westAltitude', 'westAzimuth', 'eastAltitude', 'eastAzimuth', ))

        alts, azs = moonAltAz()
        al, az = corrections(alts,azs)
        ifm.point(al,az)
        print(ifm.get_pointing())
        hpm.start_recording(n)


        point_time = time.time()

        while unix < end_t:
            unix = time.time()
            alts, azs = moonAltAz()
            al, az = corrections(alts,azs)
            ifm.point(al,az)
            get_point = ifm.get_pointing()
            valAltAz = get_point.values()
            pointFile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(unix, alts, azs, 
                                    al, az, valAltAz[0][0], valAltAz[0][1], valAltAz[1][0], valAltAz[1][1]))
            print('{}\t{}\t{}\t{}\n'.format(valAltAz[0][0], valAltAz[0][1], valAltAz[1][0], valAltAz[1][1]))

            time.sleep(dt)
            if unix - point_time >= 15*60: 
                point_time = time.time()
                volts, times = hpm.get_recording_data()
                np.savez('{}_{}'.format(dataFile, unix), v = volts, tt = times)

    vend, tend = hpm.end_recording()
    np.savez('{}_last'.format(dataFile), vv = vend, ttt = tend) 
    return volts, times


def set_timer(unix_t, func, *args):
    '''
    Create a function
    
    unix_t: float, unix time when you want function to run
    func: func, function you want to run at unix_t
    args: tuple, arguments for your function. Must have a comma
    after your last argument
    
    Example:
    set_timer(time.time() + 1, sum, (np.arange(4),) )
    returns 6
    
    set_timer(time.time() + 2, lambda x,y: 7 + x, (6,'maybe?',))
    returns 13
    
     
    Prints returned value in your function arguments. Also runs anything
    inside your function like np.save or additional print statements inside
    your function.
    
    if you're a dumbshit and your unix_t < time.time(), use variable.is_alive() 
    to see if Timer is still running and variable.cancel to cancel the Timer
    '''

    current_time = time.time()
    sec_till = unix_t - current_time
    variable = CustomTimer(sec_till, func, *args)
    variable.start()
    print(variable.join())
    

def julianDayArray(hrs, dt, initUnixT):
    """
    Makes an array of Julian dates for a specified interval

    Arguments:
    hrs: (float) Number of hours to record Julian dates. Won't get exact 
    number of hours because length of array will depend on interval time
    dt: (float) Number of minutes before recording next Julian day
    initUnixT: (float) Unix time for the first initial Julian day

    Returns:
    Two arrays: Julian day array with interval dt, Hours array initialized to 0

    """
    interval = hrs*60./dt
    secInt = dt*60.
    jds = []
    for i in range(int(interval)):
        jds.append(ugradio.timing.julian_date(initUnixT + secInt*i))
    # actual julian data to put into alt az functions
    jds = np.array(jds)
    # start from 0 and convert to hours for plot
    jdHrs = (jds - jds[0])*24
    return jds, jdHrs

def getAltAzArray(ra, dec, timeArr):
    """
    Makes an array of altitudes and azimuths for several times given an array. This uses the function findAltAz
    ra: (float) right ascention
    dec: (float) declination
    timeArr: (ndarray) array of times (julian days)
    """
    altitude = []
    azimuth = []
    for tt in timeArr:
        alt, az = findAltAz(ra, dec, jdT = tt) 
        altitude.append(alt)
        azimuth.append(az)
    return np.array(altitude), np.array(azimuth)