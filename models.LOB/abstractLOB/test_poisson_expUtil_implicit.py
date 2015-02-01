'''
Created on Jan 28, 2015

@author: weiluo
'''
import unittest
from poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC
from brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC
from pylab import plot, show
import numpy as np
class Test(unittest.TestCase):
    @unittest.SkipTest
    def test_seterr(self):
        params = {"delta_t": 0.0001, "num_time_step": 1000, "gamma": 1.0,\
                 "A": 20, "kappa":1.5, "sigma_s": 1.0, "N": 10}
        myObj = Poisson_expUtil_implicit_NeumannBC(**params)
        myObj.run(1)
        x = np.true_divide(1.0, 0.0)
        print x
        np.seterr(all='raise')
        self.assertRaisesRegexp(FloatingPointError, lambda x:  np.true_divide(1.0, 0.0))
    @unittest.SkipTest    
    def test_vary_beta1(self):
        params = {"delta_t": 0.0001, "num_time_step": 1000, "gamma": 1.0,\
                 "A": 20, "kappa":1.5, "sigma_s": 1.0, "N": 10}
        
        for curr_beta_power in np.arange(-6, -2):
            myObj = Poisson_expUtil_implicit_NeumannBC(beta = 10**curr_beta_power, **params)
            myObj.run()
            
            myObj_zero_beta = Poisson_expUtil_implicit_NeumannBC(beta=0, **params)
            myObj_zero_beta.run()
            plot(myObj.q_space, myObj.a_control[0], 'r')
            plot(myObj_zero_beta.q_space, myObj_zero_beta.a_control[0], 'b')
            show()
            
            plot(myObj.q_space, myObj.result[-1], 'r')
            plot(myObj_zero_beta.q_space, myObj_zero_beta.result[-1], 'b')
            show()  
    
    
    def test_cannot_converge_params2(self):
        
        params = {"delta_t": 0.0001, "num_time_step": 1000, "gamma": 1.0,\
                 "A": 20, "kappa":1.5, "sigma_s": 1.0, "N": 10, "beta": 0.001}
        myObj = Poisson_expUtil_implicit_NeumannBC(**params)
        myObj.run()
        
        params_zero_beta = params.copy()
        params_zero_beta['beta'] = 0
        myObj_zero_beta = Poisson_expUtil_implicit_NeumannBC(**params_zero_beta)
        myObj_zero_beta.run()
        plot(myObj.q_space, myObj.a_control[0], 'r')
        plot(myObj_zero_beta.q_space, myObj_zero_beta.a_control[0], 'b')
        show()
        
        plot(myObj.q_space, myObj.result[-1], 'r')
        plot(myObj_zero_beta.q_space, myObj_zero_beta.result[-1], 'b')
        show()
    @unittest.SkipTest
    def test_cannot_converge_params(self):
        
        params = {"delta_t": 0.0001, "num_time_step": 1000, "gamma": 1.0,\
                 "A": 20, "kappa":1.5, "sigma_s": 1.0, "N": 10, "beta": 0.000001}
        myObj = Poisson_expUtil_implicit_NeumannBC(**params)
        myObj.run()
        
        params_zero_beta = params.copy()
        params_zero_beta['beta'] = 0
        myObj_zero_beta = Poisson_expUtil_implicit_NeumannBC(**params_zero_beta)
        myObj_zero_beta.run()
        plot(myObj.q_space, myObj.a_control[0], 'r')
        plot(myObj_zero_beta.q_space, myObj_zero_beta.a_control[0], 'b')
        show()
        
        plot(myObj.q_space, myObj.result[-1], 'r')
        plot(myObj_zero_beta.q_space, myObj_zero_beta.result[-1], 'b')
        show()
    
    
    @unittest.skip("testing skipping")
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