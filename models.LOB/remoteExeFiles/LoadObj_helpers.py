'''
Created on Feb 7, 2015

@author: weiluo
'''

import pickle, sys
import numpy as np
import matplotlib.pyplot as plt 
def two_array_hist_plot(first, second, k =50, label_first="first", label_second="second", show_xlabel=""):
    upperLimit = int(max(np.max(np.array(first)),  np.max(np.array(second))))+1
    bins = np.linspace(0, upperLimit, k*upperLimit)
    plt.hist(first,bins, alpha = 0.5, color='b', label = label_first)
    plt.hist(second,bins,  alpha = 0.5, color='r', label = label_second)   
    plt.legend(loc="best")
    plt.xlabel(show_xlabel)
    plt.show()

class LoadObj(object):
    '''
    classdocs
    '''
    

    def __init__(self,  dir = "./pickleObjs/"):
        self.name = None
        self.data = None
        self.dir = dir
    def readPiclObj(self, fileName, relative = True):
        if relative:
            fileName = self.dir + fileName
        try:
            fileHandler = open(fileName, 'r')
            self.data = pickle.load(fileHandler)
        except:
            print sys.exc_info()[0]
        