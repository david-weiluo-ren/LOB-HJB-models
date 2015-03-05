'''
Created on Feb 17, 2015

@author: weiluo
'''
import unittest
from abstractLOB.poisson_expUtil_implicit import Poisson_expUtil_implicit_priceSineDrift_NeumannBC

from SimulateSingle import *
from pylab import plot, show
class Test(unittest.TestCase):
    
    def test_average_q(self):
        obj8 = Poisson_expUtil_implicit_priceSineDrift_NeumannBC\
(theta = 1,  sigma_s = 0.01, num_time_step=50, beta=0, A=0.1,\
 N=16.0, gamma=0.5, kappa=0.4, q_0=3)
        obj8.run()
        result8 = LoadSingleData(summary_mean_var_helper(obj8, 5,None,"",False))
        plot(result8.simulated_q_mean)
        show()

    @unittest.SkipTest
    def test_s_drift(self):
        options={"BC": "sameslope", "type":"poisson", "num_time_step":100}
        
        summary_mean_var(options, 20, "")
    @unittest.SkipTest
    def test_simulate_q(self):
        options={"BC": "sameslope", "type":"poisson", "num_time_step":100, "q_0":2.0}
        result = summary_mean_var(options, 200, "")
        plot(result[3][3])
        show()
        
    @unittest.SkipTest    
    def test_simulate_q2(self):
        options={'A': 0.1, 'BC': 'sameslope','N': 16.0,\
                 'beta': 0.08,'delta_t': 0.01,'gamma': 0.5,\
                 'kappa': 0.4,'num_time_step': 100,'q_0': 2.0, 'sigma_s': 0.04,\
                 'type': 'poisson'}
        result = summary_mean_var(options, 200, "")
        show()
        plot(result[3][3])
        show()
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_s_drift']
    unittest.main()