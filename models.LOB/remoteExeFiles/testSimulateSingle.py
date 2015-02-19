'''
Created on Feb 17, 2015

@author: weiluo
'''
import unittest

from SimulateSingle import summary_mean_var
class Test(unittest.TestCase):


    def test_s_drift(self):
        options={"BC": "sameslope", "type":"poisson", "num_time_step":100}
        
        summary_mean_var(options, 20, "")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_s_drift']
    unittest.main()