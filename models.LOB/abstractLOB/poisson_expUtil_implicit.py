'''
Created on Jan 28, 2015

@author: weiluo
'''
from abstractDerivedLOB import AbstractImplicitLOB_NeumannBC
import numpy as np 
from numpy import exp
from scipy.special import lambertw
import scipy as sp
class Poisson_expUtil_implicit_NeumannBC(AbstractImplicitLOB_NeumannBC):
    '''
    
    Poisson Model using exponential utility and implicit method.
    '''
    def compute_q_space(self):
        super(Poisson_expUtil_implicit_NeumannBC, self).compute_q_space()
        self._q_space  = np.linspace(-self.N, self.N, self.I)
        self._delta_q = 1
    
    @property
    def half_I(self):
        return self.N
    
    @property
    def q_space(self):
        return self._q_space
    
    @property
    def delta_q(self):
        return self._delta_q
    @property
    def result(self):
        if len(self._result) == 0:
            return self._result
        return self.user_friendly_list_of_array(self._result, 0, call_func = lambda x: -x)
    
    
    @property
    def a_control(self):
        return self._data_helper(self._index_a_control_2darray, self.extend_space-1, self.extend_space-1)
    
    @property
    def b_control(self):
        return self._data_helper(self._index_b_control_2darray, self.extend_space-1, self.extend_space-1)
   
    
    def terminal_condition(self):
        super(Poisson_expUtil_implicit_NeumannBC, self).terminal_condition()
        return np.ones(len(self.implement_q_space))

    def F_1(self, x, beta1, beta2, beta3):

        return 0 if x == self.control_upper_bound else self.A*exp(-self.kappa*x)*(-beta1 * x + beta2 * exp(-self.gamma * x) - beta3)

    def feedback_control(self, v):
        a_curr = np.zeros(self.implement_I)
        b_curr = np.zeros(self.implement_I)
        for i in xrange(self.implement_I):
            if(i == 0 or i == self.implement_I - 1):
                continue
            q = i - self.half_implement_I
            a_beta1 = self.gamma * q *self.beta * v[i]
            a_beta2 = v[i - 1]
            a_beta3 = v[i]
            
            b_beta1 = -self.gamma * q *self.beta * v[i]
            b_beta2 = v[i + 1]
            b_beta3 = v[i]
                
           
                

            
            #When exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1)) overflow. Try compute manully or another way to compute
            #without using the lambertW function.
            if abs(a_beta1)<10**(-3):
                a_curr[i] = np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)) + np.log(v[i-1]) - np.log(v[i]))
                b_curr[i] = np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)) + np.log(v[i+1]) - np.log(v[i]))
            elif q > 0:
                a_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(a_beta3, a_beta1) + np.true_divide(1, self.gamma)\
                 * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*a_beta1)\
                                       * self.gamma * a_beta2 *exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1))))
                if np.true_divide(self.gamma + self.kappa, self.kappa*b_beta1)  * self.gamma * b_beta2 *exp(-np.true_divide(self.gamma*(b_beta1 - self.kappa*b_beta3),self.kappa*b_beta1)) <= - exp(-1):
                    b_curr[i] = self.control_upper_bound
                else:
                    b_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(b_beta3, b_beta1) + np.true_divide(1, self.gamma) * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*b_beta1)  * self.gamma * b_beta2 *exp(-np.true_divide(self.gamma*(b_beta1 - self.kappa*b_beta3),self.kappa*b_beta1)), -1))


            else:
                b_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(b_beta3, b_beta1) + np.true_divide(1, self.gamma) * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*b_beta1)  * self.gamma * b_beta2 *exp(-np.true_divide(self.gamma*(b_beta1 - self.kappa*b_beta3),self.kappa*b_beta1))))
            

                if np.true_divide(self.gamma + self.kappa, self.kappa*a_beta1)  * self.gamma * a_beta2 *exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1)) <= - exp(-1):
                    a_curr[i] = self.control_upper_bound
                else:
                    a_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(a_beta3, a_beta1) + np.true_divide(1, self.gamma) * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*a_beta1)  * self.gamma * a_beta2 *exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1)), -1))
                      
            a_curr[i] = 0 if a_curr[i] < 0 else a_curr[i]
            b_curr[i] = 0 if b_curr[i] < 0 else b_curr[i]
            a_curr[i] = self.control_upper_bound if self.F_1(a_curr[i], a_beta1, a_beta2, a_beta3)>0 else a_curr[i]
            b_curr[i] = self.control_upper_bound if self.F_1(b_curr[i], b_beta1, b_beta2, b_beta3)>0 else b_curr[i]
        return [a_curr[1:-1], b_curr[1:-1]]
    
    def linear_system(self, v_curr, curr_control):
        super(Poisson_expUtil_implicit_NeumannBC, self).linear_system( v_curr, curr_control)
        a_curr, b_curr = curr_control
        co_left = self.A * self.delta_t * np.exp(-(self.kappa + self.gamma) * a_curr)
        co_right = self.A * self.delta_t * np.exp(-(self.kappa + self.gamma) * b_curr)
        co_mid = 1 + self.delta_t * ( - 0.5 * self.sigma_s**2 * self.gamma**2 * self.implement_q_space[1:-1]**2 + self.gamma * self.beta * self.A * self.implement_q_space[1:-1] \
                                      * (np.exp(-self.kappa* a_curr) * a_curr - np.exp(-self.kappa * b_curr) * b_curr)\
                                       + self.A * (np.exp(-self.kappa* a_curr) + np.exp(-self.kappa* b_curr)))
        
        
        
        eq_right = v_curr.copy()
        eq_right[0] = 0
        eq_right[-1] = 0
        return [-self.coef_at_minus_one_helper(co_left), -self.coef_at_plus_one_helper(co_right), self.coef_at_curr_helper(co_mid), eq_right]
        
    def run(self, K=None, use_cache=False):
        old_settings = np.seterr(all='raise')
        super(Poisson_expUtil_implicit_NeumannBC, self).run( K=K, use_cache=use_cache)
        np.seterr(**old_settings)
        
        
            
    def __init__(self, *args, **kwargs):
        super(Poisson_expUtil_implicit_NeumannBC, self).__init__(*args, **kwargs)
        self.half_implement_I = self.half_I + self.extend_space
        