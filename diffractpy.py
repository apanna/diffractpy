# -*- coding: utf-8 -*-
#!/usr/bin/env python
"""
Plots grating diffraction profiles and fits to gaussian
__author__   =   Alireza Panna
__status__   =   Development
__version__  =   1-0
__to-do__    =   better estimation of sigma.
__date__     =   5/25/2016

CHANGELOG:
05/26/2016       
"""
from libtiff import TIFF
from PIL import Image
import numpy as np
import os, time, sys
import easygui
import scipy.optimize, scipy.special, scipy.stats, scipy.signal
import fits
import plot
import math

PEAK_SEPERATION = 5
PEAK_HEIGHT = 75

def filter_savgol(y):
    """
    Savitzky-Golay filter
    """
    return scipy.signal.savgol_filter(np.array(y), window_length=11, polyorder=2)

def gaussian_fit(x, y, guess):
    """
    Do the gaussian fit on data
    """
    g_fit = []
    try:
        popt, pcov = scipy.optimize.curve_fit(fits.three_gaussian,\
        np.array(x), np.array(y), guess, \
        maxfev=1000*(len(x)+1))
    except:
        popt, pcov = guess, None

    for i in fits.three_gaussian(np.array(x), *popt):
        g_fit.append(i)
    return g_fit, popt, pcov
    
def write_file(o, c, p, r2, path, fname):
    """
    write optimal fit params
    """
    with open(path + fname + '.txt', 'w+') as f:
        f.write("All peak Values = " + str(p) + "\n")
        f.write("Amplitudes from fit = " + str(o[2]) + " " + str(o[5]) + " " + str(o[8]) + "\n")
        f.write("Mean values from fit = " + str(o[3]) + " " + str(o[6]) + " " + str(o[9]) + "\n")
        f.write("Sigma values from fit = " + str(o[4]) + " " + str(o[7]) + " " + str(o[10]) + "\n")
        f.write("FWHM values from fit = " + str(2.35482*o[4]) + " " + str(2.35482*o[7]) + " " + str(2.35482*o[10]) + "\n") 
        f.write("Area under gaussian curve = " + str(np.sqrt(2*math.pi)*o[2]*o[4]) + " " + \
                str(np.sqrt(2*math.pi)*o[5]*o[7]) + " " + str(np.sqrt(2*math.pi)*o[8]*o[10]) + "\n")
        f.write("R^2 of fit = " + str(r2))
    f.close()
    return
            
def finish_dialog(st):
    end_time = time.time() - st   
    easygui.msgbox("          Analysis Complete     \n" + "Time Taken: " + str(end_time) + ' s')
   
if __name__ == '__main__':
    fieldNames = ["Upper Left x:", "Upper Left y:", "Lower Right x:", "Lower Right y:"]
    crop_coords = []
    crop_str = ''
    pix = []
    avrg = []
    norm_avrg = []
    ct = 0
    directory= ''
    dir_and_file = easygui.fileopenbox(msg="Select _*.tif image", title="File Explorer", default='*', filetypes='.tif', multiple=False)
    directory_list = dir_and_file.split("\\")[:-1]
    for dir in directory_list:
        directory += dir + os.sep
    fname = dir_and_file.split("\\")[-1]
    fname_extension = fname.split(".")[-1]
    fname_no_ext = fname.split(".")[0]
    fname_no_underscore = fname_no_ext.split('_')[0]
    if directory == None:
    	sys.exit()
     
    crop_coords = easygui.multenterbox(msg='Enter Crop Co-Ordinates.', \
                                       title='Crop-Co-ordinates (Leave blank for last saved)', \
                                       fields=fieldNames, values=())
    if not os.path.exists(directory + os.sep + 'Figures'):
        os.makedirs(directory + os.sep + 'Figures')
    
    if crop_coords == ['','','','']:
    	try:
    		with open(os.getcwd() + os.sep + "crop_corners.txt", 'r') as f:
    			crop_str = f.readline()
    		crop_coords = [crop_str.split(" ")[0], crop_str.split(" ")[1], crop_str.split(" ")[2], crop_str.split(" ")[3].split("\n")[0]]
    	except IOError:
    		pass
    else:
    	with open(os.getcwd() + os.sep + "crop_corners.txt", 'w+') as f:
    		f.write(crop_coords[0] + " " + crop_coords[1] + " " + crop_coords[2] + " " + crop_coords[3])
    start_time = time.time()
    for files in os.listdir(directory):
         extension = os.path.splitext(directory + os.sep + files)[1]
         if extension == '.'+ fname_extension and files.startswith(fname_no_underscore) == True:
             print files
             pix = []
             avrg = []
             norm_avrg = []
             corr_avrg = []
             tif = TIFF.open(directory + os.sep + files, mode='r')
             im_array = tif.read_image()
             im = Image.fromarray(im_array)
             im = im.crop((int(crop_coords[0]), int(crop_coords[1]), int(crop_coords[2]), int(crop_coords[3])))
             for i, rows in enumerate(np.array(im)):
                 avrg.append(np.mean(rows))
                 pix.append(i)
#             corr_avrg = fits.linear_baseline_correction(pix, avrg, 5)
#             print corr_avrg
             for i, avr in enumerate(filter_savgol(avrg)):
                 norm_avrg.append(avr)
#             for i, corr in enumerate(corr_avrg):
#                 norm_avrg.append(corr)
             # Estimate peak and peak positions for initial guess.   
             indexes, peaks = fits.find_peaks(pix, norm_avrg, dist=PEAK_SEPERATION, threshold=PEAK_HEIGHT)
             print indexes
             if (len(indexes) == 3 and len(peaks) == 3):
                 mean_1 = indexes[0]
                 mean_2 = indexes[1]
                 mean_3 = indexes[2]
                 peak_1 = peaks[0]
                 peak_2 = peaks[1]
                 peak_3 = peaks[2]
                 sigma_1 = np.sqrt(mean_1)
                 sigma_2 = np.sqrt(mean_2) - sigma_1
                 sigma_3 = np.sqrt(mean_3) - sigma_2 - sigma_1
                 
             elif (len(indexes) == 2 and len(peaks) == 2):
                 mean_1 = indexes[0]
                 mean_2 = indexes[1]
                 mean_3 = 0
                 peak_1 = peaks[0]
                 peak_2 = peaks[1]
                 peak_3 = 0
                 sigma_1 = np.sqrt(mean_1)
                 sigma_2 = np.sqrt(mean_2) - sigma_1
                 sigma_3 = 0
                 
             elif (len(indexes) == 1 and len(peaks) == 1):
                 mean_1 = indexes[0]
                 mean_2 = 0
                 mean_3 = 0
                 peak_1 = peaks[0]
                 peak_2 = 0
                 peak_3 = 0
                 sigma_1 = np.sqrt(mean_1)
                 sigma_2 = 0
                 sigma_3 = 0
             else:
                 mean_1 = 0
                 mean_2 = 0
                 mean_3 = 0
                 peak_1 = 0
                 peak_2 = 0
                 peak_3 = 0
                 sigma_1 = 0
                 sigma_2 = 0
                 sigma_3 = 0
                
             p_guess = [0, 0, peak_1, mean_1, sigma_1, \
                          peak_2, mean_2, sigma_2, \
                          peak_3, mean_3, sigma_3]
             gauss_fit, opt, cov = gaussian_fit(pix, norm_avrg, p_guess)
             r2 = fits.r_squared(norm_avrg, gauss_fit, np.mean(norm_avrg))
             file_path = directory + os.sep + 'Figures' + os.sep 
             write_file(opt, cov, peaks, r2, file_path, files.split('.')[0])
             # Create the figure...
             plot.plot(pix, avrg, pix, gauss_fit, indexes, peaks, \
                       save='Yes', grid='Yes', x_label='Pixels', \
                       y_label='Intensity (A.U.)', \
                       save_path=directory + os.sep + 'Figures' + os.sep + \
                       files.split('.')[0] + ".jpg")
             ct = ct+1
             print '---------------------------'
    finish_dialog(start_time)

