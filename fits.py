# -*- coding: utf-8 -*-
"""
"""
__author__          =   "Alireza Panna"
__email__           =   "apanna1236@gmail.com"
__status__          =   "Development"
__date__            =   "10/22/14"
__version__         =   "0.2"

import numpy as np

def linear(x, intercept, slope):
    """
    A linear fit function:
    y-intercept                 : intercept
    slope                       : slope
    """
    return intercept+slope*x

def gaussian(x, offs, slope, peak, mu, sigma):
    """
    A gaussian fit function:
    Offset                      : offset
    baseline                    : slope
    Amplitude                   : peak
    Mean                        : mu
    Standard deviation          : sigma
    """
    return -offs-slope*x+\
            peak*np.exp(-(x-mu)**2/(2*sigma**2))
              
def three_gaussian(x, offs, slope, peak_1, mu_1, sigma_1, peak_2, mu_2, sigma_2, peak_3, mu_3, sigma_3):
    return -offs-slope*x+\
            peak_1*np.exp(-(x-mu_1)**2/(2*sigma_1**2)) + \
            peak_2*np.exp(-(x-mu_2)**2/(2*sigma_2**2)) + \
            peak_3*np.exp(-(x-mu_3)**2/(2*sigma_3**2))

def log_gaussian(x, peak, mu, sigma):        
    """
    compute the a, b, c that parameterise the parabola
    log(gauss(x, A, mu, sigma)) ~ a + bx**2 + cx.
    """
    return np.log(peak) - (x-mu)**2/(2*sigma**2)

def lorentzian(x, offs, slope, peak, mu, fwhm):
    """
    A lorentzian peak with:
    Offset                      : p[0]
    slope                       : p[1]
    Amplitude                   : p[2]
    Mean                        : p[3]
    FWHM                        : p[4]
    """
    return -offs-slope*x+\
           (peak*fwhm**2/4)*(1/((x-mu)**2+(fwhm/2)**2))

def power(x, norm, offs, const):
    """
    A power law fit with:
    Normalization               : p[0]
    Offset                      : p[2]
    Constant                    : p[3]
    """
    return norm*(x-offs)**offs + const      

def residuals(meas, fit):
    res=[]
    for i,j in zip(meas, fit):
        res.append(i-j)
    return res

def chi_squared(meas, fit, meas_sigma):
    c2=0
    for i, j, k in zip(meas, fit, meas_sigma):
        c2+=((i-j)/k)**2
    return c2

def r_squared(meas, fit, meas_mean):
    ss_res=0
    ss_tot=0
    for i, j in zip(meas, fit):
        ss_res+=(i-j)**2
        ss_tot+=(i-meas_mean)**2
    return 1-(ss_res/ss_tot)
    
def linear_baseline_correction(x, y, pts):
    """
    Perform linear baseline correction on data to flatten the tails
    """
    corrected = []
    pts = int(pts)
    slope = (np.mean(y[-pts:]) - np.mean(y[0:pts]))/(np.mean(x[-pts:]) - np.mean(x[0:pts]))
    offset = (np.mean(y[0:pts]) * np.mean(x[-pts:])  - \
              np.mean(y[-pts:]) * np.mean(x[0:pts]))/(np.mean(x[-pts:]) - np.mean(x[0:pts]));
    for i, j in enumerate(y):
        corrected.append(j-offset-slope*x[i])
        
    return np.array(corrected)

def first_diff(x, y):
    """
    Performs first order difference
    
    Parameters
    ----------
    Returns
    -------
    """
    fd = []
    fd.append(0)
    for y1, y2, x1, x2 in zip(y[0:-1], y[1:], x[0:-1], x[1:]):
        if (x2 - x1) != 0:
            fd.append((y2 - y1)/(x2 - x1))
    return np.array(fd)
    
def find_peaks(x, y, **kwargs):
    """
    Finds peaks in dataset.
    """
    try:
        peak_seperation = int(kwargs.get('dist'))
        peak_threshold  = kwargs.get('threshold')    
    except:
        peak_seperation = 15
        peak_threshold  = 150
    ct                  = 0
    previous            = 0
    current             = 0
    local_maxima_x      = []
    maxima_y            = []
    max_y               = []
    maxima_x            = []
    local_maxima_y      = []
#    x_interp = np.linspace(min(x), max(x), num = 1000)
#    y_interp = np.interp(x_interp, x, y)
    # first order difference
    diff = first_diff(x, y)
#    diff2 = first_diff(x, diff)
#    print diff2
    for i, j in zip(diff[0:-1], diff[1:]):
        # look for maxima.
        previous = current
        if ((i > 0 and j < 0) and y[ct] > peak_threshold):
            current = ct;
            if (current - previous > peak_seperation):
                local_maxima_x.append(x[ct])
                local_maxima_y.append(y[ct])
        ct = ct + 1
       
    local_maxima_y.sort()
    # Get the highest peaks of the gaussians
    if (len(local_maxima_y) >= 3):
        maxima_y = [local_maxima_y[-1], local_maxima_y[-2], local_maxima_y[-3]]
    elif (len(local_maxima_y) == 2):
        maxima_y = [local_maxima_y[-1], local_maxima_y[-2]]
    elif (len(local_maxima_y) == 1):
        maxima_y = [local_maxima_y[-1]]
    else:
        maxima_y = [0]
    # Get the peak position corresponding to the peaks
    for i in maxima_y:
        for k, j in zip(x, y):
            if (i == j):
                maxima_x.append(k)  
    # Resort according to peak position
    maxima_x.sort(reverse=False)
    for i in maxima_x:
        max_y.append(y[i])
    return  maxima_x, max_y

if __name__=="__main__":
    print "There is no place like main!"
6