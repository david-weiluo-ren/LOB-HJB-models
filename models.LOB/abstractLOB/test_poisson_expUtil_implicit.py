'''
Created on Jan 28, 2015

@author: weiluo
'''
import unittest
from poisson_expUtil_implicit import Poisson_expUtil_implicit
from brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC
from pylab import plot, show
class Test(unittest.TestCase):

    @unittest.skip
    def test_result1(self):
        myObj = Poisson_expUtil_implicit()
        myObj.run()
        
       
        self.assertEqual(len(myObj.q_space), len(myObj.result[-1]))
        
        self.assertEqual(len(myObj.q_space), len(myObj.a_control[-1]))
        
        plot(myObj.q_space, myObj.result[-1])
        show()
        plot(myObj.q_space, myObj.a_control[-1])
        show()
        
        plot(myObj.q_space, myObj.b_control[-1])
        show()
    def test_interface1(self):
        myPoisson = Poisson_expUtil_implicit()
        myBM = BrownianMotion_ExpUtil_Implicit_NeumannBC()
        '''
        for myObj in [myPoisson, myBM]:
            myObj.run()
        
       
            self.assertEqual(len(myObj.q_space), len(myObj.result[-1]))
        
            self.assertEqual(len(myObj.q_space), len(myObj.a_control[-1]))
        
            plot(myObj.q_space, myObj.result[-1])
            show()
            plot(myObj.q_space, myObj.a_control[-1])
            show()
        
            plot(myObj.q_space, myObj.b_control[-1])
            show()
        '''
        for myObj in [myPoisson, myBM]:
            myObj.run()
        
       
            self.assertEqual(len(myObj.q_space), len(myObj.result[-1]))
            self.assertEqual(len(myObj.q_space), len(myObj.a_control[-1]))
            plot(myObj.q_space, myObj.a_control[-1])
        
            plot(myObj.q_space, myObj.b_control[-1])
        show()

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_result1']
    unittest.main()