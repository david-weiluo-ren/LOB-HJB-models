'''
Created on Feb 18, 2015

@author: weiluo
'''
import unittest
from SimulateComparison import LoadMultiData
import pickle
class Test(unittest.TestCase):


    def test_index_ourRange(self):
        
        fileHandler = open("pickleObjs/A_0.1_simulate_num_20000_kappa_0.4_BC_sameslope_sigma_s_0.04_num_time_step_40000_N_10.0_delta_t_0.005_beta_0.04_q_0_2.0_type_poisson_gamma_0.5_comparison",\
'r')
        myData = LoadMultiData(pickle.load(fileHandler))

    def test_index_ourRange2(self):
        
        fileHandler = open("pickleObjs/A_0.1_simulate_num_2000_kappa_0.4_random_q_0_true_BC_sameslope_sigma_s_0.04_num_time_step_10000_N_16.0_delta_t_0.01_beta_0.04_q_0_2.0_type_poisson_gamma_0.5_comparison",\
                   'r')
        myData = LoadMultiData(pickle.load(fileHandler))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()