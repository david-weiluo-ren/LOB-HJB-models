'''
Created on Feb 7, 2015

@author: weiluo
'''
import unittest
from LoadObj_helpers import LoadObj, two_array_hist_plot


class Test(unittest.TestCase):


    def test_load1(self):
        myObj = LoadObj()
        myObj.readPiclObj("gmm_1_sgm01_kpp3_A05_N_7_beta1")
        two_array_hist_plot(myObj.data[0], myObj.data[2])
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_load1']
    unittest.main()