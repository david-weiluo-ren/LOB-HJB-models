'''
Created on Mar 29, 2015

@author: weiluo
'''
import unittest
from poisson_OU_explicit import Poisson_explicit_OU_LOB
from poisson_OU_implicit import Poisson_OU_implicit_truncateControlAtZero 
from pylab import plot, show
import numpy as np
import pickle
class Test(unittest.TestCase):
    def test_OU_sameRandomness_moreAccurate(self):
        explicitObj =Poisson_explicit_OU_LOB(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=20000,\
                                  half_S=2.5, delta_t=0.00001, N=10.0, half_I_S=125,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0)
        randomSource = explicitObj.generate_random_source()
        implicitObj = Poisson_OU_implicit_truncateControlAtZero(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=20000,\
                                  half_S=2.5, delta_t=0.00001, N=10.0, half_I_S=125,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0, new_weight=0.1, \
                                  abs_threshold_power = -8, rlt_threshold_power = -7)
        explicitObj.run()
        implicitObj.run()
        explicitObj.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        implicitObj.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        explicitHandler = open('explicitObj1', 'w')
        pickle.dump(explicitObj, explicitHandler)
        implicitHandler = open('implicitObj1', 'w')
        pickle.dump(implicitObj, implicitHandler)

        
        plot(explicitObj.s, 'r')
        plot(implicitObj.s, 'b')
        show()
        plot(explicitObj.simulate_control_a, 'r')
        plot(implicitObj.simulate_control_a, 'b')
        show()
        plot(implicitObj.simulate_control_a, 'b')
        show()
        plot(explicitObj.simulate_control_b, 'r')
        plot(implicitObj.simulate_control_b, 'b')
        show()
        plot(implicitObj.simulate_control_b, 'b')
        show()
        plot(np.array(explicitObj.simulate_control_a) + np.array(explicitObj.s[:-1]), 'r')
        plot(np.array(implicitObj.simulate_control_a) + np.array(implicitObj.s[:-1]), 'b')
        show()
        plot(np.array(implicitObj.simulate_control_a) + np.array(implicitObj.s[:-1]), 'b')
        show()
        plot(-np.array(explicitObj.simulate_control_b) + np.array(explicitObj.s[:-1]), 'r')
        plot(-np.array(implicitObj.simulate_control_b) + np.array(implicitObj.s[:-1]), 'b')
        show()
        plot(-np.array(implicitObj.simulate_control_b) + np.array(implicitObj.s[:-1]), 'b')
        show()

    @unittest.SkipTest
    def test_OU_sameRandomness(self):
        explicitObj =Poisson_explicit_OU_LOB(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=20000,\
                                  half_S=2.5, delta_t=0.00001, N=10.0, half_I_S=250,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0)
        randomSource = explicitObj.generate_random_source()
        implicitObj = Poisson_OU_implicit_truncateControlAtZero(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=20000,\
                                  half_S=2.5, delta_t=0.00001, N=10.0, half_I_S=250,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0, new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        explicitObj.run()
        implicitObj.run()
        explicitObj.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        implicitObj.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        #explicitHandler = open('explicitObj', 'w')
        #pickle.dump(explicitObj, explicitHandler)
        #implicitHandler = open('implicitObj', 'w')
        #pickle.dump(implicitObj, implicitHandler)

        
        plot(explicitObj.s, 'r')
        plot(implicitObj.s, 'b')
        show()
        plot(explicitObj.simulate_control_a, 'r')
        plot(implicitObj.simulate_control_a, 'b')
        show()
        plot(implicitObj.simulate_control_a, 'b')
        show()
        plot(explicitObj.simulate_control_b, 'r')
        plot(implicitObj.simulate_control_b, 'b')
        show()
        plot(implicitObj.simulate_control_b, 'b')
        show()
        plot(np.array(explicitObj.simulate_control_a) + np.array(explicitObj.s[:-1]), 'r')
        plot(np.array(implicitObj.simulate_control_a) + np.array(implicitObj.s[:-1]), 'b')
        show()
        plot(np.array(implicitObj.simulate_control_a) + np.array(implicitObj.s[:-1]), 'b')
        show()
        plot(-np.array(explicitObj.simulate_control_b) + np.array(explicitObj.s[:-1]), 'r')
        plot(-np.array(implicitObj.simulate_control_b) + np.array(implicitObj.s[:-1]), 'b')
        show()
        plot(-np.array(implicitObj.simulate_control_b) + np.array(implicitObj.s[:-1]), 'b')
        show()
        
        
    @unittest.SkipTest
    def test_OU_sameRandomness2(self):
        explicitObj =Poisson_explicit_OU_LOB(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=10000,\
                                  half_S=3, delta_t=0.0001, N=10.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0)
        randomSource = explicitObj.generate_random_source()
        implicitObj = Poisson_OU_implicit_truncateControlAtZero(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=10000,\
                                  half_S=3, delta_t=0.0001, N=10.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0, new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        explicitObj.run()
        implicitObj.run()
        explicitObj.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        implicitObj.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        plot(explicitObj.s, 'r')
        plot(implicitObj.s, 'b')
        show()
        plot(explicitObj.simulate_control_a, 'r')
        plot(implicitObj.simulate_control_a, 'b')
        show()
        plot(implicitObj.simulate_control_a, 'b')
        show()
        plot(explicitObj.simulate_control_b, 'r')
        plot(implicitObj.simulate_control_b, 'b')
        show()
        plot(implicitObj.simulate_control_b, 'b')
        show()
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_OU_sameRandomness']
    unittest.main()