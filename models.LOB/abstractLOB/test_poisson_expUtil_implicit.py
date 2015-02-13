'''
Created on Jan 28, 2015

@author: weiluo
'''
import unittest
from poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC, Poisson_expUtil_implicit_sameSlopeBC
from brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC
from pylab import plot, show
import numpy as np
import matplotlib.pyplot as plt
import time

class Test(unittest.TestCase):
    '''
    @unittest
    def test_large_space_sameSlope(self):
        params = {"num_time_step": 1000, "beta": 0.4, "gamma": 0.1,\
                 "A": 0.1, "kappa":0.2, "sigma_s": 0.1, "N": 20}
        myObj = Poisson_expUtil_implicit_sameSlopeBC(**params)
        myObj.run()
        for arr in myObj.b_control:
            plot(arr)
        show()
    '''
    
    def test_large_space(self):
        params = {"num_time_step": 1000, "beta": 0.4, "gamma": 0.1,\
                 "A": 0.1, "kappa":0.2, "sigma_s": 0.5, "N": 20}
        myObj = Poisson_expUtil_implicit_sameSlopeBC(**params)
        myObj.run()
        
        for arr in myObj.b_control:
            plot(arr)
        for arr in myObj.a_control:
            plot(arr)
        show()
    
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
    
    @unittest.SkipTest
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
    def test_simulate1(self):
        myObj = Poisson_expUtil_implicit_NeumannBC(gamma = 1, beta = 0.05, A=5, num_time_step=2000)
        myObj.run()
        myObj.simulate_forward()
        
        plt.subplot(221)
        plt.plot(myObj.simulate_control_a, 'r')
        plt.plot(myObj.simulate_control_b, 'b')
        
        plt.subplot(222)
        plt.plot(myObj.q)
        
        plt.subplot(223)
        plt.plot(myObj.x)
        
        plt.subplot(224)
        plt.plot(myObj._a_control[-1])
        plt.plot(myObj._b_control[-1])
        
        plt.show()
    @unittest.SkipTest
    def test_simulate_sameSlopeBC_bug1(self):
        myObj = Poisson_expUtil_implicit_sameSlopeBC(gamma = 1, sigma_s=0.5, beta = 0.0005, A=5, num_time_step=1000, iter_max=500)
        myObj.run()
        
        for arr in myObj.a_control:
            plot(arr)
        show()
    @unittest.SkipTest
    def test_simulate_sameSlopeBC_bug2(self):
        myObj = Poisson_expUtil_implicit_sameSlopeBC(gamma = 1, sigma_s=0.5, beta = 0, A=5, num_time_step=1000, iter_max=500)
        myObj.run()
        
        for arr in myObj.a_control:
            plot(arr)
        show()
    @unittest.SkipTest   
    def test_simulate_sameSlopeBC_bug_tryNeumanBC(self):
        myObj = Poisson_expUtil_implicit_NeumannBC(gamma = 1, sigma_s=0.5, beta = 0, A=5, num_time_step=1000, iter_max=500)
        myObj.run()
        
        for arr in myObj.a_control:
            plot(arr)
        show()
    @unittest.SkipTest    
    def test_simulate_sameSlopeBC_bug_tryNeumanBC2(self):
        myObj = Poisson_expUtil_implicit_NeumannBC(gamma = 1, sigma_s=0.5, beta = 0.0005, A=5, num_time_step=500, iter_max=500)
        myObj.run()
        
        for arr in myObj.a_control:
            plot(arr)
        show()
        
    @unittest.SkipTest    
    def test_simulate_Gueant_params_smallBeta(self):
        myObj = Poisson_expUtil_implicit_NeumannBC(gamma = 0.1,\
                     sigma_s=0.3, kappa = 0.3, beta = 0.005, A=0.9, num_time_step=20000, iter_max=500)
        myObj.run()
        
        for arr in myObj.a_control:
            plot(arr)
        show()
        
        for i in xrange(len(myObj.a_control.T)):
            if not i%5:
                plot(myObj.a_control.T[i])
        show()
    @unittest.SkipTest
    def test_simulate_Gueant_params(self):
        myObj = Poisson_expUtil_implicit_NeumannBC(gamma = 0.1,\
                     sigma_s=0.3, kappa = 0.3, beta = 0, A=0.9, num_time_step=20000, iter_max=500)
        myObj.run()
        
        for arr in myObj.a_control:
            plot(arr)
        show()
        
        for i in xrange(len(myObj.a_control.T)):
            if not i%5:
                plot(myObj.a_control.T[i])
        show()
    @unittest.SkipTest    
    def test_simulate_Gueant_params_not_truncated_at_zero(self):
        myObj = Poisson_expUtil_implicit_NeumannBC(gamma = 0.1,\
                     sigma_s=0.3, kappa = 0.3, beta = 0, A=0.9, num_time_step=5000, iter_max=500,truncated_at_zero=False)
        myObj.run()
        
        for arr in myObj.a_control:
            plot(arr)
        show()
        
        for i in xrange(len(myObj.a_control.T)):
            if not i%5:
                plot(myObj.a_control.T[i])
        show()
    @unittest.SkipTest    
    def test_simulate_sameSlopeBC1(self):
        myObj = Poisson_expUtil_implicit_sameSlopeBC(gamma = 1, beta = 0.05, A=5, num_time_step=2000)
        myObj.run()
        myObj.simulate_forward()
        
        plt.subplot(221)
        plt.plot(myObj.simulate_control_a, 'r')
        plt.plot(myObj.simulate_control_b, 'b')
        
        plt.subplot(222)
        plt.plot(myObj.q)
        
        plt.subplot(223)
        plt.plot(myObj.x)
        
        plt.subplot(224)
        plt.plot(myObj._a_control[-1])
        plt.plot(myObj._b_control[-1])
        
        plt.show()   
    @unittest.SkipTest     
    def test_simulate_sameSlopeBC_noSparse1(self):
        myObj = Poisson_expUtil_implicit_sameSlopeBC(use_sparse = False, gamma = 1, beta = 0.05, A=5, num_time_step=2000)
        myObj.run()
        myObj.simulate_forward()
        
        plt.subplot(221)
        plt.plot(myObj.simulate_control_a, 'r')
        plt.plot(myObj.simulate_control_b, 'b')
        
        plt.subplot(222)
        plt.plot(myObj.q)
        
        plt.subplot(223)
        plt.plot(myObj.x)
        
        plt.subplot(224)
        plt.plot(myObj._a_control[-1])
        plt.plot(myObj._b_control[-1])
        
        plt.show()     
    @unittest.SkipTest
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
    @unittest.SkipTest
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
        
    def setUp(self):
        self.startTime = time.time()

    def tearDown(self):
        t = time.time() - self.startTime
        print "%s: %.3f" % (self.id(), t)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_result1']
    #unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(Test)
    unittest.TextTestRunner(verbosity=1).run(suite)