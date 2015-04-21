'''
Created on Apr 17, 2015

@author: weiluo
'''
import unittest

from poisson_OU_explicit_marketOrder import Poisson_OU_explicit_marketOrder
from pylab import plot, show
import numpy as np
class Test(unittest.TestCase):

    @unittest.SkipTest
    def testOneStepBack(self):
        myObj = Poisson_OU_explicit_marketOrder()
        myObj.one_step_back(self.v_init, 0)
    @unittest.SkipTest
    def testRun(self):
        myObj = Poisson_OU_explicit_marketOrder(l=0.001)
        myObj.run()
    @unittest.SkipTest    
    def testRun2(self):
        myObj = Poisson_OU_explicit_marketOrder(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=10000,\
                                  half_S=3, delta_t=0.0001, N=8.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0, l=0.01,\
                                  half_I=80)
        myObj.run()
        myObj.simulate_forward()
        plot( myObj.simulate_accumulate_xi)
        show()
        plot(myObj.simulate_control_a)
        plot(myObj.simulate_control_b)
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.s, 'g')
        show()
    @unittest.SkipTest
    def testRun3(self):
        myObj = Poisson_OU_explicit_marketOrder(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=100000,\
                                  half_S=3, delta_t=0.00001, N=8.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0, l=0.1,\
                                  half_I=80)
        myObj.run()
        myObj.simulate_forward()
        plot( myObj.simulate_accumulate_xi)
        show()
        plot(myObj.simulate_control_a)
        plot(myObj.simulate_control_b)
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.s, 'g')
        show()    
    @unittest.SkipTest
    def testRun4(self):   
        obj8_tmp = Poisson_OU_explicit_marketOrder(A=20, s_0=0, kappa=1.5, sigma_s=1.5, num_time_step=20000,\
                                  half_S=3.0, delta_t=0.00005, N=6.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=2.5,gamma=1.0, s_long_term_mean=0, l=0.05,\
                                  half_I=80)
        obj8_tmp.run()    

    def testRun_checkInventoryUnderflow(self):
        objNormal = Poisson_OU_explicit_marketOrder(A=20, s_0=2.0, kappa=1.5, sigma_s=3/np.sqrt(2), num_time_step=10000,\
                                  half_S=4.0, delta_t=0.0001, N=6.0, half_I_S=60,\
                                  beta=0, q_0=0.0, alpha=5,gamma=1.0, s_long_term_mean=2.0, l=0.05,\
                                  half_I=60)
        objNormal.run()
        tmp=plot(np.array(objNormal._a_control[7000]).T)
        show()
        tmp=plot(np.array(objNormal._a_control[7000]))
        show()
        objNormal.simulate_forward(q_0=1.0)
if __name__ == "__main__":

    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()