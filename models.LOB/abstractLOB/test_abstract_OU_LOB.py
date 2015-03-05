'''
Created on Mar 3, 2015

@author: weiluo
'''
import unittest
from poisson_OU_explicit import *
from pylab import plot, show, ylim
import numpy as np
from matplotlib.pyplot import xlim
class Test(unittest.TestCase):

    def testPoisson_OU_LOB2(self):
        myObj = Poisson_explicit_OU_LOB(num_time_step=1000, N=15, q_0=8, \
                                 alpha = 10, delta_t=0.0001, half_S=1.5, half_I_S=45, kappa=1.5, beta=0.02, A=20, gamma=0.01, sigma_s=2)
        myObj.run()
        
        myObj.simulate_forward()
        plot(myObj.s, 'g')
        

        show()
        

        plot(np.array(myObj.simulate_control_a), 'r')
        plot(np.array(myObj.simulate_control_b), 'b')
        ylim([-2,12])
        show()
        plot(myObj.q)
        show()



    @unittest.SkipTest
    def testPoisson_OU_LOB(self):
        myObj = Poisson_explicit_OU_LOB(num_time_step=100000, N=5, q_0=0, \
                                 alpha = 10, delta_t=0.00001, kappa=1.5, beta=0.05, A=20, gamma=1, sigma_s=3)
        myObj.run()
        
        myObj.simulate_forward(100000)
        plot(myObj.s, 'g')
        #plot(np.array(myObj.simulate_control_a)+np.array(myObj.s)[1:], 'r')
        #plot(np.array(myObj.s)[1:] - np.array(myObj.simulate_control_b), 'b')

        show()
        #plot(np.array(myObj.s)-4, 'g')

        plot(np.array(myObj.simulate_control_a), 'r')
        plot(np.array(myObj.simulate_control_b), 'b')

        show()
        plot(myObj.q)
        show()
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()