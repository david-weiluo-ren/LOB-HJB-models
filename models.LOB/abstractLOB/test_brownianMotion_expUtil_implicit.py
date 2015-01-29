'''
Created on Jan 28, 2015

@author: weiluo
'''
import unittest
from brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC
from pylab import plot, show
class Test(unittest.TestCase):


    def test_result1(self):
        myObj = BrownianMotion_ExpUtil_Implicit_NeumannBC()
        myObj.run()
        self.assertEqual(len(myObj.q_space), len(myObj.result[-1]))
        
        self.assertEqual(len(myObj.q_space), len(myObj.a_control[-1]))
        tmp = myObj.result
        plot(tmp[-1])
        show()

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()