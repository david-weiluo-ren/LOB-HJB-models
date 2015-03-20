'''
Created on Mar 17, 2015

@author: weiluo
'''
import unittest
from poisson_OU_implicit import Poisson_OU_implicit, Poisson_OU_implicit_truncateControlAtZero
import numpy as np
from pylab import plot, show
class Test(unittest.TestCase):

    @unittest.SkipTest
    def test_implicit2(self):
        myObj= Poisson_OU_implicit(A=20, s_0=5, kappa=1.5, sigma_s=3.0, num_time_step=1000,\
                                  half_S=5, delta_t=0.0001, N=10.0, half_I_S=150,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=5.0, new_weight=0.1)
        myObj.run()
        
        
    @unittest.SkipTest  
    def test_implicit_truncateControlZero_nonZeroBeta(self):
        myObj= Poisson_OU_implicit_truncateControlAtZero(A=20, s_0=0, kappa=1.5, sigma_s=3.0, num_time_step=200,\
                                  half_S=3, delta_t=0.001, N=10.0, half_I_S=50,\
                                  beta=0.4, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=0.0)
        myObj.run()
        myObj.simulate_forward()
        plot(myObj.s)
        show()
        plot(myObj.simulate_control_a,'r')
        print myObj.simulate_control_a

        plot(myObj.simulate_control_b,'b')
        print myObj.simulate_control_b
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.s[:-1], 'g')
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.simulate_price_a_test, 'g')
        plot(myObj.simulate_price_b_test, 'y')
        plot(myObj.s[:-1], 'black')

        show()
        plot(myObj.q[:-1])
        show()    
        
        
    def test_implicit_price(self):
        myObj= Poisson_OU_implicit(A=20, s_0=0, kappa=1.5, sigma_s=3.0, num_time_step=200,\
                                  half_S=3, delta_t=0.001, N=10.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=0.0)
        myObj.run()
        myObj.simulate_forward()
        plot(myObj.s)
        show()
        plot(myObj.simulate_control_a,'r')
        print myObj.simulate_control_a

        plot(myObj.simulate_control_b,'b')
        print myObj.simulate_control_b
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.s[:-1], 'g')
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.simulate_price_a_test, 'g')
        plot(myObj.simulate_price_b_test, 'y')
        plot(myObj.s[:-1], 'black')

        show()
        plot(myObj.q[:-1])
        show()
    @unittest.SkipTest
    def test_implicit_truncateControlZero(self):
        myObj= Poisson_OU_implicit_truncateControlAtZero(A=20, s_0=0, kappa=1.5, sigma_s=3.0, num_time_step=200,\
                                  half_S=3, delta_t=0.001, N=10.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=0.0)
        myObj.run()
        myObj.simulate_forward()
        plot(myObj.s)
        show()
        plot(myObj.simulate_control_a,'r')
        print myObj.simulate_control_a

        plot(myObj.simulate_control_b,'b')
        print myObj.simulate_control_b
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.s[:-1], 'g')
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.simulate_price_a_test, 'g')
        plot(myObj.simulate_price_b_test, 'y')
        plot(myObj.s[:-1], 'black')

        show()
        plot(myObj.q[:-1])
        show()
        
        
    @unittest.SkipTest
    def test_implicit1(self):
        myObj= Poisson_OU_implicit(A=20, s_0=0, kappa=1.5, sigma_s=3.0, num_time_step=200,\
                                  half_S=3, delta_t=0.001, N=10.0, half_I_S=50,\
                                  beta=0, q_0=0.0, alpha=10.0,gamma=1.0, s_long_term_mean=0.0)
        myObj.run()
        myObj.simulate_forward()
        plot(myObj.s)
        show()
        plot(myObj.simulate_control_a,'r')
        plot(myObj.simulate_control_b,'b')
        show()
        plot(myObj.simulate_price_a, 'r')
        plot(myObj.simulate_price_b, 'b')
        plot(myObj.s[:-1], 'g')
        show()
        plot(myObj.q[:-1])
        show()

    @unittest.SkipTest
    def test_expNegOptimalControl(self):
        myObj = Poisson_OU_implicit(N=2, half_I_S=1, extend_space=0)
        
        myObj.run()
        print myObj._result[-1]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_expNegOptimalControl']
    unittest.main()