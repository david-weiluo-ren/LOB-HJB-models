'''
Created on Feb 7, 2015

@author: weiluo
'''

import pickle, sys
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.pyplot import xlabel
def two_array_hist_plot(first, second, k =50, label_first="first", label_second="second", show_xlabel=""):
    upperLimit = int(max(np.max(np.array(first)),  np.max(np.array(second))))+1
    bins = np.linspace(0, upperLimit, k*upperLimit)
    plt.hist(first,bins, alpha = 0.5, color='b', label = label_first)
    plt.hist(second,bins,  alpha = 0.5, color='r', label = label_second)   
    plt.legend(loc="best")
    plt.xlabel(show_xlabel)
    plt.show()

def draw_q(data, jump = 200):
    index = np.arange(0, len(data.non_zero_beta_obj.simulated_q_mean), jump)
    ax = plt.gca()
    try:
        time_step = data.non_zero_beta_obj.options["delta_t"]
    except:
        time_step = 0.01
        
    try:
        num_steps = data.non_zero_beta_obj.options["num_time_step"]
    except:
        num_steps = 100
    
    ax.errorbar(np.linspace(0, num_steps*time_step, len(data.non_zero_beta_obj.simulated_q_mean))[index],\
            data.non_zero_beta_obj.simulated_q_mean[index],\
            yerr=np.vstack([np.sqrt(data.non_zero_beta_obj.simulated_q_var)[index],
                              np.sqrt(data.non_zero_beta_obj.simulated_q_var)[index]]),  label="with price impact")

    index = np.arange(0, len(data.zero_beta_obj.simulated_q_mean), 200)

    ax.errorbar(np.linspace(0,num_steps*time_step, len(data.zero_beta_obj.simulated_q_mean))[index],\
            data.zero_beta_obj.simulated_q_mean[index],\
            yerr=np.vstack([np.sqrt(data.zero_beta_obj.simulated_q_var)[index],
                              np.sqrt(data.zero_beta_obj.simulated_q_var)[index]]), color='g', label="no price impact")
    plt.legend(loc="best",prop={'size':9})
    plt.xlabel("t")
    plt.ylabel("inventory") 


    plt.draw()
def draw_s(data, jump = 200):
    index = np.arange(0, len(data.non_zero_beta_obj.simulated_s_mean), jump)
    ax = plt.gca()
    try:
        time_step = data.non_zero_beta_obj.options["delta_t"]
    except:
        time_step = 0.01
        
    try:
        num_steps = data.non_zero_beta_obj.options["num_time_step"]
    except:
        num_steps = 100
    
    ax.errorbar(np.linspace(0, num_steps*time_step, len(data.non_zero_beta_obj.simulated_s_mean))[index],\
            data.non_zero_beta_obj.simulated_s_mean[index],\
            yerr=np.vstack([np.sqrt(data.non_zero_beta_obj.simulated_s_var)[index],
                              np.sqrt(data.non_zero_beta_obj.simulated_s_var)[index]]), label="with price impact")

    index = np.arange(0, len(data.zero_beta_obj.simulated_s_mean), 200)

    ax.errorbar(np.linspace(0,num_steps*time_step, len(data.zero_beta_obj.simulated_s_mean))[index],\
            data.zero_beta_obj.simulated_s_mean[index],\
            yerr=np.vstack([np.sqrt(data.zero_beta_obj.simulated_s_var)[index],
                              np.sqrt(data.zero_beta_obj.simulated_s_var)[index]]), color='g', label="no price impact")
    plt.legend(loc="best",prop={'size':9})
    plt.xlabel("t")
    plt.ylabel("price") 
    plt.draw()
   
class LoadObj(object):
    '''
    classdocs
    '''
    

    def __init__(self,  directory = "./pickleObjs/"):
        self.name = None
        self.data = None
        self.directory = directory
    def readPiclObj(self, fileName, relative = True):
        if relative:
            fileName = self.dir + fileName
        try:
            fileHandler = open(fileName, 'r')
            self.data = pickle.load(fileHandler)
        except:
            print sys.exc_info()[0]
        