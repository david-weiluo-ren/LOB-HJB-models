'''
Created on Jan 28, 2015

@author: weiluo
'''
from abstractDerivedLOB import AbstractImplicitLOB_NeumannBC
import numpy as np 
from numpy import exp
from scipy.special import lambertw
import scipy as sp
class Poisson_expUtil_implicit(AbstractImplicitLOB_NeumannBC):
    '''
    
    Poisson Model using exponential utility and implicit method.
    '''
    def compute_q_space(self):
        super(Poisson_expUtil_implicit, self).compute_q_space()
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
    
    def terminal_condition(self):
        super(Poisson_expUtil_implicit, self).terminal_condition()
        return -1 * np.ones(len(self.implement_q_space))

    def F_1(self, x, beta1, beta2, beta3):

        return 0 if x == self.control_upper_bound else self.A_*exp(-self.kappa_*x)*(-beta1 * x + beta2 * exp(-self.gamma_ * x) - beta3)

    def feedback_control(self, v):
        a_curr = np.zeros(self.implement_I)
        b_curr = np.zeros(self.implement_I)
        for i in xrange(self.spaceGridNum_total_):
            if(i == 0 or i == self.spaceGridNum_total_ - 1):
                continue
            q = i - self.spaceGridNum_oneSide_
            a_beta1 = self.gamma_ * q *self.beta_ * v[i]
            a_beta2 = v[i - 1]
            a_beta3 = v[i]
            
            b_beta1 = -self.gamma_ * q *self.beta_ * v[i]
            b_beta2 = v[i + 1]
            b_beta3 = v[i]
                



            if abs(a_beta1)<10**(-3):
                a_curr[i] = np.true_divide(1, self.gamma_)*(np.log(1+np.true_divide(self.gamma_, self.kappa_)) + np.log(v[i-1]) - np.log(v[i]))
                b_curr[i] = np.true_divide(1, self.gamma_)*(np.log(1+np.true_divide(self.gamma_, self.kappa_)) + np.log(v[i+1]) - np.log(v[i]))
            elif q > 0:
                a_curr[i] = np.true_divide(1, self.kappa_) - np.true_divide(a_beta3, a_beta1) + np.true_divide(1, self.gamma_) * sp.real(lambertw( np.true_divide(self.gamma_ + self.kappa_, self.kappa_*a_beta1)  * self.gamma_ * a_beta2 *exp(-np.true_divide(self.gamma_*(a_beta1 - self.kappa_*a_beta3),self.kappa_*a_beta1))))
                if np.true_divide(self.gamma_ + self.kappa_, self.kappa_*b_beta1)  * self.gamma_ * b_beta2 *exp(-np.true_divide(self.gamma_*(b_beta1 - self.kappa_*b_beta3),self.kappa_*b_beta1)) <= - exp(-1):
                    b_curr[i] = self.control_upper_bound
                else:
                    b_curr[i] = np.true_divide(1, self.kappa_) - np.true_divide(b_beta3, b_beta1) + np.true_divide(1, self.gamma_) * sp.real(lambertw( np.true_divide(self.gamma_ + self.kappa_, self.kappa_*b_beta1)  * self.gamma_ * b_beta2 *exp(-np.true_divide(self.gamma_*(b_beta1 - self.kappa_*b_beta3),self.kappa_*b_beta1)), -1))


            else:
                b_curr[i] = np.true_divide(1, self.kappa_) - np.true_divide(b_beta3, b_beta1) + np.true_divide(1, self.gamma_) * sp.real(lambertw( np.true_divide(self.gamma_ + self.kappa_, self.kappa_*b_beta1)  * self.gamma_ * b_beta2 *exp(-np.true_divide(self.gamma_*(b_beta1 - self.kappa_*b_beta3),self.kappa_*b_beta1))))
            

                if np.true_divide(self.gamma_ + self.kappa_, self.kappa_*a_beta1)  * self.gamma_ * a_beta2 *exp(-np.true_divide(self.gamma_*(a_beta1 - self.kappa_*a_beta3),self.kappa_*a_beta1)) <= - exp(-1):
                    a_curr[i] = self.control_upper_bound
                else:
                    a_curr[i] = np.true_divide(1, self.kappa_) - np.true_divide(a_beta3, a_beta1) + np.true_divide(1, self.gamma_) * sp.real(lambertw( np.true_divide(self.gamma_ + self.kappa_, self.kappa_*a_beta1)  * self.gamma_ * a_beta2 *exp(-np.true_divide(self.gamma_*(a_beta1 - self.kappa_*a_beta3),self.kappa_*a_beta1)), -1))
                      
            a_curr[i] = 0 if a_curr[i] < 0 else a_curr[i]
            b_curr[i] = 0 if b_curr[i] < 0 else b_curr[i]
            a_curr[i] = self.control_upper_bound if self.F_1(a_curr[i], a_beta1, a_beta2, a_beta3)>0 else a_curr[i]
            b_curr[i] = self.control_upper_bound if self.F_1(b_curr[i], b_beta1, b_beta2, b_beta3)>0 else b_curr[i]
        return [a_curr, b_curr]
    def     
    def __init__(self, params):
        '''
        Constructor
        '''
        