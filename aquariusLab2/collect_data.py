import numpy as np
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

def plot_fft_power(v, div, N, arrA = True, arrB = True, dualmode = True, nbloc = 1, nsamp = 16000, saveFile = False, fileName = 'data'):

    '''
    Arguments:
    v: (string) volt_range same as pico.capture_data
    div: (int) Divide the 62.5 MHz sample clock by this number for sampling
    arrA: (bool, ndarray) Optional. If no ndarray input, it will take data from pico.capture_data in real time. 
    if there is an array input, it will use that data instead
    arrB: (bool, ndarray) Same explanation as arrA                          
    dual_mode: (bool) samples from both A and B. If False, samples only from A.                                                                        
    nbloc: (int) number of blocks (each with nsample data points)                                                                                      
    nsamp: (int) number of samples. default is 16000                                                                                                   
    fileName: (string) the name of file                                                                                                                
    saveFile: (bool) if True, will save data to two files fileName_A and fileName_B
    

    Returns:
    Plots for A and B of fft and power
    '''


    ll = collect_data(v, div, dualmode, nbloc, nsamp)
    if not isinstance(arrA, ndarray):
        arrA = ll[0]
    if not isinstance(arrB, ndarray):
        arrB= ll[1]
        
    fftA_real = np.fft.fft(arrA[:N]).real
    fftB_real = np.fft.fft(arrB[:N]).real
    fftA_imag = np.fft.fft(arrA[:N]).imag
    fftB_imag = np.fft.fft(arrB[:N]).imag
    

    vsamp = 62.5/div
    time = np.linspace((-N/2.)/vsamp, (N/2. - 1)/vsamp,N)
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./N),N)


    fig, axes = plt.subplots(1,2, figsize = (12,10))
    
    axes[1,1].plot(time, arrA[:N], label = 'port A', alpha = 0.5)
    axes[1,1].plot(time, arraB[:N], label = 'port B', alpha = 0.5)

    return [arrA,arrB]
