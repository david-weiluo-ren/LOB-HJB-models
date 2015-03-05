'''
Created on Mar 1, 2015

@author: weiluo
'''
import unittest

from poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC, Poisson_expUtil_implicit_sameSlopeBC,Poisson_expUtil_implicit_priceSineDrift_NeumannBC
from pylab import plot, show
import numpy as np
import matplotlib.pyplot as plt
import time
class Test(unittest.TestCase):
    
    
    def test_control_adaption(self):
        myObj = Poisson_expUtil_implicit_priceSineDrift_NeumannBC\
        (theta = 10, sigma_s = 0.00001, num_time_step=500, beta=0, A=0.1,\
        N=16.0, gamma=0.5, kappa=0.4, q_0=3)
        myObj.run()
        myObj.simulate_forward()
       

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()