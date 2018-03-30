# Calculate the range of local fringe frequencies that I should see in the data. 
# Look at fourier transform of data to check then do some filtering. Large freq/smaller freq?? 

def calculateFringes():


    return 0


def frequencies(vsamp, N):
    '''
    Array of frequencies for given sampling frequency and N number of points
    Arguments: 
    div: (int) divisor chosen when collecting data with pico sampler to determing sampling frequency
    N: (int) number of points 

    Returns:
    Array of frequencies of length N
    '''
    freq = np.linspace(-vsamp/2, (vsamp/2)*(1-2./N), N)
    return freq

def power(data):
    '''
    Calculate the power spectra
    Arguments:

    data: (array) 
    N: (int) number of samples to use

    Returns: 
    Power shifted using np.fft.fftshift
    '''
    fft_ = np.fft.fft(data)
    power = abs(fft_)**2
    return np.fft.fftshift(power)

def avg_power(data, nChunk):
    #split data up
    

def plot_FFT():

    plt.plot()

def plot_POWER():

    plt.plot()


