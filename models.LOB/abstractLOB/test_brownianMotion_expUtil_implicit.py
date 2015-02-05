'''
Created on Jan 28, 2015

@author: weiluo
'''
import unittest
from brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC
from pylab import plot, show
import matplotlib.pyplot as plt
class Test(unittest.TestCase):
    @unittest.SkipTest
    def test_notConverge2(self):
        params = {"delta_t": 0.0001, "num_time_step": 1000, "gamma": 1.0,\
                 "A": 20, "kappa":1.5, "sigma_s": 1.0, "N": 3, "beta": 0.1,\
                 "extend_space":7}
        myObj = BrownianMotion_ExpUtil_Implicit_NeumannBC(**params)
        myObj.run(8000)
    
        plot(myObj.result[-1])
        show()
        plot(myObj._a_control[-1])
        show()
    
    @unittest.SkipTest
    def test_notConverge1(self):
        params = {"delta_t": 0.0001, "num_time_step": 1000, "gamma": 5.0,\
                 "A": 20, "kappa":1.5, "sigma_s": 1.0, "N": 10, "beta": 0.01}
        myObj = BrownianMotion_ExpUtil_Implicit_NeumannBC(**params)
        myObj.run()
    
        plot(myObj.result[-1])
        show()
        plot(myObj._a_control[-1])
        show()
    @unittest.SkipTest
    def test_extend_space1(self):
        myObj = BrownianMotion_ExpUtil_Implicit_NeumannBC(extend_space = 5, beta = 0.0001, sigma_s = 0.1)
        self.assertEqual(5 * 100, myObj.extend_space)
        myObj.run()
        plot(myObj.a_control[-1])
        show()
        plot(myObj._a_control[-1])
        show()
    
    
    @unittest.SkipTest
    def test_extend_space2(self):
        myObj = BrownianMotion_ExpUtil_Implicit_NeumannBC()
        self.assertEqual(2 * 100, myObj.extend_space)
        myObj.run()
        plot(myObj.a_control[-1])
        show()
        plot(myObj._a_control[-1])
        show()  
        
    def test_forward1(self):
        myObj = BrownianMotion_ExpUtil_Implicit_NeumannBC()
        self.assertEqual(2 * 100, myObj.extend_space)
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
        plt.plot(myObj.a_control[-1])
        plt.plot(myObj.b_control[-1])
        
        plt.show()

    @unittest.SkipTest
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