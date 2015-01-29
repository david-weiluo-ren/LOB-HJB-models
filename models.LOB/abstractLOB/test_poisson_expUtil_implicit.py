'''
Created on Jan 28, 2015

@author: weiluo
'''
import unittest
from poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC
from brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC
from pylab import plot, show
class Test(unittest.TestCase):

    def test_constructor(self):
        params = {"new_weight":0.2, "beta":0.0001}
        myObj = Poisson_expUtil_implicit_NeumannBC(extend_space=5, **params)
        self.assertEqual(0.2, myObj.new_weight)
        self.assertEqual(0.0001, myObj.beta)
        self.assertEqual(5, myObj.extend_space)
        myObj.run(use_cache=True)

        plot(myObj.q_space, myObj.a_control[-1])
        show()
        plot(myObj._data[myObj._index_a_control_2darray][-1])
        show()
        plot(myObj._data[myObj._index_result_2darray][-1])
        show()

        print myObj.implement_q_space[:10]
        print myObj.implement_q_space[-10:]
        
        
    @unittest.skip("testing skipping")
    def test_result1(self):
        myObj = Poisson_expUtil_implicit_NeumannBC()
        myObj.run()
        
       
        self.assertEqual(len(myObj.q_space), len(myObj.result[-1]))
        
        self.assertEqual(len(myObj.q_space), len(myObj.a_control[-1]))
        
        plot(myObj.q_space, myObj.result[-1])
        show()
        plot(myObj.q_space, myObj.a_control[-1])
        show()
        
        plot(myObj.q_space, myObj.b_control[-1])
        show()
    @unittest.skip("testing skipping")
    def test_interface1(self):
        myPoisson = Poisson_expUtil_implicit_NeumannBC()
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