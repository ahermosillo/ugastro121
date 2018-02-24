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

def transform_2D_3D(data, nblocs, nchunks, N):
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
    if len(data[0,:]) % nchunks != 0:
        raise Exception("Cannot split into these number of chunks. Try a divisor of {}".format(len(data[0,:])))
    else:
        new_comp = data.reshape((nblocs, nchunks, N/nchunks))
    return new_comp

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


def avg_3d_power(data, nblock, nchunk, nsample, N, div = 1):
    '''
    Calculates averaged power spectra for data split into chunks in a 3d array

    Arguments: 
    data: (array) 3D array 
    nblock: (int) number of blocks of data
    nchunk: (int) number of chunks to split up data into 
    nsample: (int) number of samples originally
    N: (int) number of samples to use
    div: (int) divisor set when collecting the data
    '''

    samp_len = nsample/nchunk
    
    vsamp = 62.5/div
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./samp_len), samp_len)
    
    comp_3D = transform_2D_3D(d1420_comp, nblock, nchunk, nsample)
    newpower = power(comp_3D,1, N)
    power_avg = (newpower.sum(axis = 1)/nchunk).sum(axis =0)/nblock
    plt.plot(freq, np.fft.fftshift(power_avg), label = 'average')
    plt.xlim(-3,3)


def avg_power(compData, div, N):
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
    plt.plot(freq, np.fft.fftshift(sum_power), 'k.-')
    plt.xlabel('Frequencies MHz')
    plt.ylabel('Power')
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
