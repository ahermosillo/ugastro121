import numpy as np
import astropy as ap
from astropy import units as u
import matplotli.pyplot as plt
import ugradio
import time

ifm = ugradio.interf.Interferometer()
hpm = ugradio.hp_multi.HP_Multimeter()

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


def findAltAz(ra, dec, jd_array):
    """
    Give array of altitude and azimuth given a certain ra, dec and Julian day array
    
    Arguments: 
    ra: (float) ra of source 
    dec: (float) dec of source
    jd_array: (ndarray) array of Julian days to get alt and az for 

    Returns:
    Two arrays. Alt, az
    """
    altArr = []
    azArr = []
    for j in jd_array:
        raP, decP = ugradio.coord.precess(ra, dec, jd = j)
        alt, az = ugradio.coord.get_altaz(raP, decP, jd =j)
        altArr.append(alt)
        azArr.append(az)
    return np.array(altArr), np.array(azArr)


def moonAltAz(jd_array):
    """
    Gives array of altitude and azimuth of the moon given Julian days (from lat, long of Campbell hall)

    Arguments:
    jd_array: (ndarray) Julian day array to get alt and az at those times

    Returns:
    Two Arrays. Alt, Az
    """
    altArr = []
    azArr = []
    for j in jd_array:
        raM,decM = ugradio.coord.moonpos(jd = j)
        alt, az = ugradio.coord.get_altaz(raM, decM, jd = j)
        altArr.append(alt)
        azArr.append(az)
    return np.array(altArr), np.array(azArr)

def sunAltAz(jd_array):
    """
    Gives array of altitude and azimuth of the moon given Julian days 

    Arguments: 
    jd_array: (ndarray) Julian day array to get alt and az at those times

    Returns: 
    Two arrays. Alt, Az
    """
    altArr = []
    azArr = []
    for j in jd_array:
        raS, decS = ugradio.coord.sunpos(jd = j)
        alt, az = ugradio.coord.get_altaz(raS, decS, jd = j)
        altArr.append(alt)
        azArr.append(az)
    return np.array(altArr), np.array(azArr)

def corrections(altitude, azimuth):
    """
    Corrects the altitude and azimuth of source since telescope range is
    Alt: 5 - 175; Az: 90 - 300

    Arguments:
    altitude: (ndarray) 
    azimuth: (ndarray)

    Returns:
    Two corrected arrays. Alt, Az
    """
    for al, az, i in zip(altitude, azimuth, range(len(altitude))):
        if az < 90:
            azimuth[i] = az + 180
            altitude[i] = abs(180 - altitude[i])
        elif az > 300:
            azimuth[i] = az - 180
            altitude[i] = abs(180 - altitude[i])
    return altitude, azimuth

def collect_data(ra, dec, hrs, dt, n, dataFile, pointfileName):
    """
    Collects voltage data from multimeter using ugradio interferometer and multimeter modules
    
    Arguments:
    ra: (float) ra of source 
    dec: (float) dec of source
    hrs: (float) hours to collect data for
    dt: (float) minutes before interferometer moves to new alt, az
    n: (float) input in start_recording. Seconds before new data value is collected
    dataFile: (string) Name of data files with voltages
    pointfileName: (string) Name of file with coordinates where interf. pointed

    Returns:
    File with coordinates, files with voltages at each Julian day
    Print Statement each time telescope moves to new location.
    """


    with open('{}'.format(pointfileName), 'w') as pointFile:
        pointFile.write('{}\t{}\t{}\t{}\n'.format('westAltitude', 'westAzimuth', 
'eastAltitude', 'eastAzimuth'))
        
        unix = ugradio.timing.unix_time()
        jds, jdsHr = julianDayArray(hrs, dt, unix)
        alts, azs = findAltAz(ra, dec, jds)
        correctAlt, correctAz = corrections(alts, azs)

        ct = 0
        for al, az in zip(correctAlt, correctAz):
            ifm.point(al, az)
            get_point = ifm.get_pointing()
            valAltAz = get_point.values()
            pointFile.write('{}\t{}\t{}\t{}\n'.format(valAltAz[0][0], valAltAz[0][1], 
valAltAz[1][0], valAltAz[1][1]))
            print('Now we are looking here:\n{}\t{}\t{}\t{}\n'.format(valAltAz[0][0], 
valAltAz[0][1], valAltAz[1][0], valAltAz[1][1]))

            hpm.start_recording(n)
            start_time = time.time()
            time.sleep(dt*60)
            current_time = time.time()
            if current_time - start_time >= dt*60:
                volts, times = hpm.end_recording()
                np.savez('{}_{}'.format(dataFile, ct), v = volts, tt = times)
                print('After {} seconds, data saved to {}_{}'.format(current_time - start_time, dataFile, ct))   
                ct += 1
