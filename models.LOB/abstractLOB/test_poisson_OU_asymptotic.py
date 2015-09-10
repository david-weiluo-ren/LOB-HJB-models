'''
Created on May 26, 2015

@author: weiluo
'''
import unittest
from poisson_OU_asymptotic import Poisson_OU_asymptotic, Poisson_OU_asymptotic_backwards, Poisson_OU_asymptotic_fullBackwards_average, Poisson_OU_asymptotic_RK4
from poisson_OU_implicit import Poisson_OU_implicit
from pylab import *
class Test(unittest.TestCase):

    @unittest.SkipTest
    def test_shapeOfAsymptotic(self):
        myObj_normalAlpha_midA_noBeta_longTermMean3_longTime = Poisson_OU_implicit(A=10, s_0=3, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=2, delta_t=0.001, N=10.0, half_I_S= 20000,\
                                  beta=-0, q_0=0.0, alpha=10,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        asymptotic_obj = Poisson_OU_asymptotic(myObj_normalAlpha_midA_noBeta_longTermMean3_longTime,0,0,5.67)
        asymptotic_obj.run()
        plot(asymptotic_obj.value_function[1000:-1000])
        #plot(asymptotic_obj.derivative_atGreaterThanMu[:30])
        show()
        #print asymptotic_obj.value_function
        print asymptotic_obj.derivative_atGreaterThanMu
    
    
   
    @unittest.SkipTest
    #I think it's pretty good.
    def test_shapeOfAsymptotic_backwards(self):
        myObj_normalAlpha_midA_noBeta_longTermMean3_longTime = Poisson_OU_implicit(A=10, s_0=3, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=1, delta_t=0.001, N=10.0, half_I_S=50000,\
                                  beta=-0, q_0=0.0, alpha=10,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        asymptotic_obj = Poisson_OU_asymptotic_backwards(myObj_normalAlpha_midA_noBeta_longTermMean3_longTime,-0.05,0,5.67)
        asymptotic_obj.run()
        plot(asymptotic_obj.value_function[1000:-1000])
        #plot(asymptotic_obj.derivative_atGreaterThanMu[:30])
        show()
        #print asymptotic_obj.value_function
        print asymptotic_obj.derivative_atGreaterThanMu
    @unittest.SkipTest
    def test_shapeOfAsymptotic_backwards2(self):
        myObj_normalAlpha_midA_noBeta_longTermMean3_longTime = Poisson_OU_implicit(A=10, s_0=3, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=2, delta_t=0.001, N=10.0, half_I_S=20000,\
                                  beta=-0, q_0=0.0, alpha=10,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        asymptotic_obj = Poisson_OU_asymptotic_backwards(myObj_normalAlpha_midA_noBeta_longTermMean3_longTime,-0.045,0,5.6)
        asymptotic_obj.run()
        plot(asymptotic_obj.implement_s_space[5000:-5000],asymptotic_obj.value_function[5000:-5000])
        #plot(asymptotic_obj.derivative_atGreaterThanMu[:30])
        show()
        #print asymptotic_obj.value_function
        print asymptotic_obj.derivative_atGreaterThanMu
        
    @unittest.SkipTest
    #Pretty nice as well
    def test_shapeOfAsymptotic_fullBackwards_average(self):
        myObj_normalAlpha_midA_noBeta_longTermMean3_longTime = Poisson_OU_implicit(A=10, s_0=3, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=1, delta_t=0.001, N=10.0, half_I_S=50000,\
                                  beta=-0, q_0=0.0, alpha=10,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        asymptotic_obj = Poisson_OU_asymptotic_fullBackwards_average(myObj_normalAlpha_midA_noBeta_longTermMean3_longTime,-0.05,0,5.67)
        asymptotic_obj.run()
        plot(asymptotic_obj.value_function[1000:-1000])
        #plot(asymptotic_obj.derivative_atGreaterThanMu[:30])
        show()
        #print asymptotic_obj.value_function
        #print asymptotic_obj.derivative_atGreaterThanMu    
    @unittest.SkipTest
    def test_shapeOfAsymptotic_fullBackwards_average2(self):
        myObj_normalAlpha_midA_noBeta_longTermMean3_longTime = Poisson_OU_implicit(A=10, s_0=3, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=1.5, delta_t=0.001, N=10.0, half_I_S=100000,\
                                  beta=-0, q_0=0.0, alpha=10,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        asymptotic_obj = Poisson_OU_asymptotic_fullBackwards_average(myObj_normalAlpha_midA_noBeta_longTermMean3_longTime,-0.01,0,5.60)
        asymptotic_obj.run()
        plot(asymptotic_obj.value_function[1:-1])
        #plot(asymptotic_obj.derivative_atGreaterThanMu[:30])
        show()
        #print asymptotic_obj.value_function
        #print asymptotic_obj.derivative_atGreaterThanMu      
        
    @unittest.SkipTest
    #Pretty Nice
    def test_shapeOfAsymptotic_RK4(self):
        myObj_normalAlpha_midA_noBeta_longTermMean3_longTime = Poisson_OU_implicit(A=10, s_0=3, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=1, delta_t=0.001, N=10.0, half_I_S=50000,\
                                  beta=-0, q_0=0.0, alpha=10,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.1, \
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        asymptotic_obj = Poisson_OU_asymptotic_RK4(myObj_normalAlpha_midA_noBeta_longTermMean3_longTime,-0.05,0,5.67)
        asymptotic_obj.run()
        plot(asymptotic_obj.value_function[1000:-1000])
        #plot(asymptotic_obj.derivative_atGreaterThanMu[:30])
        show()
        
    
    def test_shapeOfAsymptotic_RK4_2(self):
        myObj_largeAlpha_largeA_noBeta_longTermMean3_longTime_finerTimeForAsymptotic = Poisson_OU_implicit(A=100, s_0=3, kappa=1.5, sigma_s=6.0, num_time_step=10000,\
                                  half_S=1.7, delta_t=0.001, N=8.0, half_I_S=100000,\
                                  beta=-0, q_0=0.0, alpha=40,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.02, lambda_tilde=0.2,iter_max=2000,\
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
        l=1000
        asymptotic_obj_test = Poisson_OU_asymptotic_RK4(myObj_largeAlpha_largeA_noBeta_longTermMean3_longTime_finerTimeForAsymptotic,0.0,0,52.968)
        asymptotic_obj_test.run()
        plot(asymptotic_obj_test.implement_s_space[l:-l], asymptotic_obj_test.value_function[l:-l])
        show()
    
    
    
    
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()