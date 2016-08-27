# -*- coding: utf-8 -*-
"""
Created on Fri May 27 22:06:57 2016
@author: Ali
"""
import matplotlib.pyplot as plt

def plot (x, y, *args,  **kwargs):
    """
    Plots upto three (x, y) datasets on a single figure. 
    """
    try:
        x1          = args[0]
        y1          = args[1]
        x2          = args[2]
        y2          = args[3]
        x_label     = kwargs.get('x_label')
        y_label     = kwargs.get('y_label')
        grid        = kwargs.get('grid')
        save        = kwargs.get('save')
        save_path   = kwargs.get('save_path')
    except:
        x1      = []
        y1      = []
        x2      = []
        y2      = []
        x_label = ''
        y_label = ''
        grid    = 'No'
        save    = 'No'
        
    plt.plot(x, y, 'r+', x1, y1, 'b-', linewidth=1.0, ms=5)
    plt.plot(x2, y2, 'ko', ms=5)
    plt.xlabel(x_label, fontsize=22)
    plt.ylabel(y_label, fontsize=22)
    plt.legend(loc='upper left')
    if grid == 'Yes' or grid == 'yes':
        plt.grid(linestyle='--', color='gray')
    if (save == 'Yes' or save == 'yes'):
        plt.savefig(save_path, bbox_inches='tight', dpi=80)
        plt.close()
    else:
        plt.show()