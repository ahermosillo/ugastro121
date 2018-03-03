import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import astropy as ap
import ugradio

def collect_data(v, div, dualmode = True, nbloc = 1, nsamp = 16000,  saveFile = False, fileName = 'data'):
    '''
    Collects and saves data from A and B port of pico sampler. Uses ugradio.pico.capture_data().

    Arguments:

    v: (string) volt_range same as pico.capture_data
    div: (int) Divide the 62.5 MHz sample clock by this number for sampling
    dual_mode: (bool) samples from both A and B. If False, samples only from A.
    nbloc: (int) number of blocks (each with nsample data points)
    nsamp: (int) number of samples. default is 16000
    fileName: (string) the name of file
    saveFile: (bool) if True, will save data to two files fileName_A and fileName_B
    
    Returns:
    Statement that data is saved to name of files
    '''

    A = ugradio.pico.capture_data(v, divisor = div, dual_mode = dualmode, nsamples = nsamp, nblocks = nbloc)[:nsamp]
    B = ugradio.pico.capture_data(v, divisor = div, dual_mode = dualmode, nsamples = nsamp, nblocks = nbloc)[nsamp:]
    if saveFile:
        np.savetxt('{}_A.txt'.format(fileName), A)
        np.savetxt('{}_B.txt'.format(fileName), B)
        print('Saved A port data to {}_A.txt'.format(fileName))
        print('Saved B port data to {}_B.txt'.format(fileName))
    return [A, B]

def make_complex(data, nblocs):
    '''
    Takes in data from pico sampler and returns complex array 

    Arguments: 
    data: (array) 1D array of format [A1, A2, ... , AN, B1, ... , BN]
    nblocs: (int) number of blocks set when collecting data with pico sampler

    Returns:
    2D array with shape (nblocs, nsamples)
    '''
    split = np.array(np.split(data, 2*nblocs))
    real =  np.array(np.arange(0, 2*nblocs) >= nblocs, dtype = bool)
    imag = np.invert(real)
    return split[real] + 1j*split[imag]

def transform_2D_3D(data, nblocs, nchunks, nsamples = 16000):
    '''
    Transforms 2d array to 3d array

    Arguments: 
    data: (array) 2D array 
    nblocs: (int) number blocks set when collecting data with pico sampler
    nchunks: (int) number of chunks to break up each row into
    N: (int) number of samples

    Returns:
    3D array with shape (nblocs, nchunks, N)
    '''
    if nsamples % nchunks != 0:
        raise Exception("Cannot split into these number of chunks. Try a divisor of {}".format(nsamples))
    else:
        new_comp = data.reshape((nblocs, nchunks, nsamples/nchunks))
    return new_comp

def frequencies(div, N):
    '''
    Array of frequencies for given sampling frequency and N number of points
    Arguments: 
    div: (int) divisor chosen when collecting data with pico sampler to determing sampling frequency
    N: (int) number of points 

    Returns:
    Array of frequencies of length N
    '''

    vsamp = 62.5/div
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./N), N)
    return freq

def power(data,N):
    '''
    Calculate the power spectra
    Arguments:

    data: (array) 
    N: (int) number of samples to use

    Returns: 
    Power shifted using np.fft.fftshift
    '''
    fft_ = np.fft.fft(data[:,:N])
    power = abs(fft_)**2
    return np.fft.fftshift(power)

def avg_power(compData, div, N):
    vsamp = 62.5/div
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./N),N)
    sum_power = 0
    for i in range(len(compData)):
        fft_ = np.fft.fft(compData[i][:N])
        sum_power += abs(fft_)**2
    sum_power = sum_power/len(compData)
    return freq, np.fft.fftshift(sum_power)

def avg_3d_power(data, nblock, nchunk, nsample, N, div = 1):
    samp_len = nsample/nchunk

    vsamp = 62.5/div
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./samp_len), samp_len)

    comp_3D = transform_2D_3D(data, nblock, nchunk, nsamples=nsample)
    newpower = power(comp_3D, N)
    power_avg = (newpower.sum(axis = 1)/nchunk).sum(axis =0)/nblock
    return freq, power_avg

def avg_3d_power_plot(data, nblock, nchunk, nsample, N, color, lab, div = 1):
    '''
    Calculates averaged power spectra for data split into chunks in a 3d array

    Arguments: 
    data: (array) 3D array 
    nblock: (int) number of blocks of data
    nchunk: (int) number of chunks to split up data into 
    nsample: (int) number of samples originally
    N: (int) number of samples to use
    div: (int) divisor set when collecting the data
    
    Returns:
    plot of the averaged Power using the 3D array
    '''

    samp_len = nsample/nchunk
    
    vsamp = 62.5/div
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./samp_len), samp_len)
    
    comp_3D = transform_2D_3D(data, nblock, nchunk, nsamples=nsample)
    newpower = power(comp_3D, N)
    power_avg = (newpower.sum(axis = 1)/nchunk).sum(axis =0)/nblock
    plt.plot(freq, power_avg, '{}.-'.format(color), label = '{}'.format(lab))
    plt.xlabel('Frequency MHz', fontsize = 15)
    plt.ylabel('Power', fontsize = 15)
    plt.tick_params(labelsize = 15)
    plt.grid()

def avg_power_plot(compData, div, N, color, lab):
    '''
    Calculate the average power for data
    
    Arguments:
    compData: (array)
    div: (int) divisor chosen when collecting data set
    N: (int) number of samples to use

    Returns:
    Average Power spectra 
    '''
    
    vsamp = 62.5/div
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./N),N)
    sum_power = 0
    for i in range(len(compData)):
        fft_ = np.fft.fft(compData[i][:N])
        sum_power += abs(fft_)**2
    sum_power = sum_power/len(compData)
    plt.plot(freq, np.fft.fftshift(sum_power), '{}.-'.format(color), label = '{}'.format(lab))
    plt.xlabel('Frequency MHz', fontsize = 15)
    plt.ylabel('Power', fontsize = 15)
    plt.tick_params(labelsize = 15)
    plt.grid()
    
    print("Sampling at a frequency {}".format(vsamp))
    print("Sampling with {} number of samples".format(N))

def plot_power(compData, div, N):
    '''
    Plots the Power spectra for each block and average
    
    Arguments:
    compData: (array)
    div: (int) divisor chosen when collecting data set
    N: (int) number of samples to use
    
    Returns:                                                                                                                         
    Average Power spectra
    '''

    vsamp = 62.5/div
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./N),N)
    sum_power = 0
    for i in range(len(compData)):
#     freq, notest_fft = ugradio.dft.dft(notest_comp[i][:N], time, freq, vsamp = 62.5)
        fft_ = np.fft.fft(compData[i][:N])
        sum_power += abs(fft_)**2
        power = abs(fft_)**2
        plt.plot(freq, np.fft.fftshift(power), label = 'block {}'.format(i))
    plt.plot(freq, np.fft.fftshift(sum_power)/len(compData), label = 'average')
    plt.xlabel('Frequencies $MHz$')
    plt.ylabel('Power')
    
    print("Sampling at a frequency {}".format(vsamp))
    print("Sampling with {} number of samples".format(N))

## ======================= Coordinate Conversions ==========================================

'''
az = azimuth 
alt = altitude
LST = Local Sidereal TIme 
Phi = Terrestrial altitude 
ha = LST - ra
dec = declination 
ra = right ascension 
l = longitude galactic coordinate 
b = langitude galactic coordinate
'''


def convert_az_alt_to_ha_dec(az,alt,phi):
    x0 = np.cos(np.radians(alt))*np.cos(np.radians(az))
    x1 = np.cos(np.radians(alt))*np.sin(np.radians(az))
    x2 = np.sin(np.radians(alt))

    R = np.matrix([[-np.sin(np.radians(phi)),0,np.cos(np.radians(phi))],[0,-1,0],[np.cos(np.radians(phi)),0,np.sin(np.radians(phi))]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(np.transpose(R),x)

    ha = float(np.degrees(np.arctan2(xp[1],xp[0]))) + 360
    dec = float(np.degrees(np.arcsin(xp[2]))) 

    return ha, dec

def convert_ha_dec_to_az_alt(ha,dec,phi):
    x0 = np.cos(np.radians(dec))*np.cos(np.radians(ha))
    x1 = np.cos(np.radians(dec))*np.sin(np.radians(ha))
    x2 = np.sin(np.radians(dec))

    R = np.matrix([[-np.sin(np.radians(phi)),0,np.cos(np.radians(phi))],[0,-1,0],[np.cos(np.radians(phi)),0,np.sin(np.radians(phi))]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(R,x)

    az = float(np.degrees(np.arctan2(xp[1],xp[0]))) 
    alt = float(np.degrees(np.arcsin(xp[2]))) 

    return az, alt

def convert_ra_dec_to_ha_dec(ra,dec,LST):
    x0 = np.cos(np.radians(dec))*np.cos(np.radians(ra))
    x1 = np.cos(np.radians(dec))*np.sin(np.radians(ra))
    x2 = np.sin(np.radians(dec))

    R = np.matrix([[np.cos(np.radians(LST)),np.sin(np.radians(LST)),0],[np.sin(np.radians(LST)),-np.cos(np.radians(LST)),0],[0,0,1]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(R,x)

    ha = float(np.degrees(np.arctan2(xp[1],xp[0])))
    dec = float(np.degrees(np.arcsin(xp[2])))
    return ha, dec

def convert_ha_dec_to_ra_dec(ha,dec,LST):
    x0 = np.cos(np.radians(dec))*np.cos(np.radians(ha))
    x1 = np.cos(np.radians(dec))*np.sin(np.radians(ha))
    x2 = np.sin(np.radians(dec))
        
    R = np.matrix([[np.cos(np.radians(LST)),np.sin(np.radians(LST)),0],[np.sin(np.radians(LST)),-np.cos(np.radians(LST)),0],[0,0,1]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(np.transpose(R),x)

    ra = float(np.degrees(np.arctan2(xp[1],xp[0])))
    dec = float(np.degrees(np.arcsin(xp[2])))
    return ra, dec

def convert_ra_dec_to_galactic1950(ra,dec): 
    x0 = np.cos(np.radians(dec))*np.cos(np.radians(ra))
    x1 = np.cos(np.radians(dec))*np.sin(np.radians(ra))
    x2 = np.sin(np.radians(dec))

    R = np.matrix([[-0.066989,-0.872756,-0.483539],[0.492728,-0.450347,0.744585],[-0.867601,-0.188375,0.460200]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(R,x)

    l = float(np.degrees(np.arctan2(xp[1],xp[0]))) + 360
    b = float(np.degrees(np.arcsin(xp[2])))

    return l ,b

def convert_galactic1950_to_ra_dec(l,b): 
    x0 = np.cos(np.radians(b))*np.cos(np.radians(l))
    x1 = np.cos(np.radians(b))*np.sin(np.radians(l))
    x2 = np.sin(np.radians(b))

    R = np.matrix([[-0.066989,-0.872756,-0.483539],[0.492728,-0.450347,0.744585],[-0.867601,-0.188375,0.460200]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(np.transpose(R),x)

    ra =  float(np.degrees(np.arctan2(xp[1],xp[0])))
    dec = float(np.degrees(np.arcsin(xp[2])))

    return ra, dec

def convert_ra_dec_to_galactic2000(ra,dec): 
    x0 = np.cos(np.radians(dec))*np.cos(np.radians(ra))
    x1 = np.cos(np.radians(dec))*np.sin(np.radians(ra))
    x2 = np.sin(np.radians(dec))

    R = np.matrix([[-0.054876,-0.873437,-0.483835],[0.494109,-0.444830,0.746982],[-0.867666,-0.198076,0.455984]])

    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(R,x)

    l = float(np.degrees(np.arctan2(xp[1],xp[0])) ) + 360
    b = float(np.degrees(np.arcsin(xp[2])))

    return l,b

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


def convert_ra_dec_to_az_alt(ra,dec,LST,phi):
    x0 = np.cos(np.radians(dec))*np.cos(np.radians(ra))
    x1 = np.cos(np.radians(dec))*np.sin(np.radians(ra))
    x2 = np.sin(np.radians(dec))

    R0 = np.matrix([[np.cos(np.radians(LST)),np.sin(np.radians(LST)),0],[np.sin(np.radians(LST)),-np.cos(np.radians(LST)),0],[0,0,1]])
    R1 = np.matrix([[-np.sin(np.radians(phi)),0,np.cos(np.radians(phi))],[0,-1,0],[np.cos(np.radians(phi)),0,np.sin(np.radians(phi))]])



    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(R1,np.dot(R0,x))

    az =  float(np.degrees(np.arctan2(xp[1],xp[0])))
    alt = float(np.degrees(np.arcsin(xp[2])))

    return az, alt


def convert_az_alt_to_ra_dec(az,alt,LST,phi):
    x0 = np.cos(np.radians(alt))*np.cos(np.radians(az))
    x1 = np.cos(np.radians(alt))*np.sin(np.radians(az))
    x2 = np.sin(np.radians(alt))

    R0 = np.matrix([[np.cos(np.radians(LST)),np.sin(np.radians(LST)),0],[np.sin(np.radians(LST)),-np.cos(np.radians(LST)),0],[0,0,1]])
    R1 = np.matrix([[-np.sin(np.radians(phi)),0,np.cos(np.radians(phi))],[0,-1,0],[np.cos(np.radians(phi)),0,np.sin(np.radians(phi))]])



    x = np.matrix([[x0],[x1],[x2]])

    xp = np.dot(np.transpose(R0),np.dot(np.transpose(R1),x))

    ra =  float(np.degrees(np.arctan2(xp[1],xp[0])))
    dec = float(np.degrees(np.arcsin(xp[2])))

    return ra, dec
