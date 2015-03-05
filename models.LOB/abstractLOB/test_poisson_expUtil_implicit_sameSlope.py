'''
Created on Feb 16, 2015

@author: weiluo
'''
import unittest
from poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC, Poisson_expUtil_implicit_sameSlopeBC


class Test(unittest.TestCase):
    
    
    def test_s_drift(self):
        params = {"num_time_step": 500}
        myObj = Poisson_expUtil_implicit_sameSlopeBC(**params)
        myObj.combine_solve_forward()
        
    @unittest.SkipTest
    def test_negative_value_function(self):
        params = {"A":0.5,"delta_t":0.02, "sigma_s": 0.04,\
                  "kappa": 0.4, "N":30, "beta":0.1, "gamma":0.5,\
                  "q_0":5.0, "num_time_step":40000}
        myObj = Poisson_expUtil_implicit_sameSlopeBC(**params)
        myObj.run()

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()