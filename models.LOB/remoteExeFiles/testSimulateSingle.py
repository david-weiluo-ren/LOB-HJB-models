'''
Created on Feb 17, 2015

@author: weiluo
'''
import unittest

from SimulateSingle import summary_mean_var
from pylab import plot, show
class Test(unittest.TestCase):

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
        
        
    def test_simulate_q2(self):
        options={'A': 0.1, 'BC': 'sameslope','N': 16.0,\
                 'beta': 0.08,'delta_t': 0.01,'gamma': 0.5,\
                 'kappa': 0.4,'num_time_step': 100,'q_0': 2.0, 'sigma_s': 0.04,\
                 'type': 'poisson'}
        result = summary_mean_var(options, 200, "")
        plot(result[3][3])
        show()
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_s_drift']
    unittest.main()