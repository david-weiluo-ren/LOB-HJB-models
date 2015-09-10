'''
Created on May 26, 2015

@author: weiluo
'''
import numpy as np
class Poisson_OU_asymptotic(object):
    '''
    classdocs
    '''


    def __init__(self, target_object = None, derivative_at_mu = 0, valueFunction_at_mu = 0, ergodic_const = 0):
        self.initialized = False
        if target_object is not None:
            self.initialized = True
            self.A = target_object.A
            self.gamma = target_object.gamma
            self.kappa = target_object.kappa
            self.lambda_tilde = target_object.lambda_tilde            
            self.alpha = target_object.alpha
            self.s_long_term_mean = target_object.s_long_term_mean
            self.sigma_s = target_object.sigma_s
            self.half_S = target_object.half_S
            self.half_I_S = target_object.half_I_S
            self.I_S = target_object.I_S
            self.implement_s_space = np.copy(target_object.implement_s_space)
            self.implement_S = target_object.implement_S
            self.half_implement_S = (self.implement_S - 1) / 2
            self.delta_s = target_object.delta_s
            self.M = np.true_divide(self.A, self.kappa + self.gamma) * (1 + np.true_divide(self.gamma, self.kappa)) ** (-np.true_divide(self.kappa, self.gamma))
            self.value_function = np.zeros(self.implement_S)
            self.value_function[self.half_implement_S] = valueFunction_at_mu
        self.derivative_atGreaterThanMu = []
        self.derivative_atLessThanMu = []

        self.derivative_at_mu = derivative_at_mu 
        self.valueFunction_at_mu = valueFunction_at_mu
        self.ergodic_const = ergodic_const

    def run(self):
        if not self.initialized:
            raise Exception("The Asymptotic Object hasn't been initialized by a Poisson_OU object.")
        self.run_greaterThanMu()
        self.run_lessThanMu()
    def updatePart_helper(self, old_value, current_s):
            return self.gamma * old_value ** 2 - 2 * np.true_divide(self.alpha, self.sigma_s**2) * (self.s_long_term_mean - current_s) * old_value\
                - np.true_divide(2 * self.M, self.sigma_s ** 2) * (np.exp(-self.kappa * (self.s_long_term_mean - current_s)) + np.exp(self.kappa * (self.s_long_term_mean - current_s))) \
                + np.true_divide(2 * self.ergodic_const, self.sigma_s ** 2)
    def run_greaterThanMu(self):
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S+i]
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atGreaterThanMu[-1]
            update_part = self.updatePart_helper(old_value, current_s)
            new_value = old_value + self.delta_s * update_part
            self.derivative_atGreaterThanMu.append(new_value)
        self.value_function[self.half_implement_S+1:] = self.valueFunction_at_mu + np.cumsum(self.delta_s * np.array(self.derivative_atGreaterThanMu))
    def run_lessThanMu(self):
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S  - i]
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atLessThanMu[-1]
            update_part = self.updatePart_helper(old_value, current_s)
            new_value = old_value - self.delta_s * update_part
            self.derivative_atLessThanMu.append(new_value)
        self.value_function[:self.half_implement_S] = (self.valueFunction_at_mu - np.cumsum(self.delta_s * np.array(self.derivative_atLessThanMu)))[::-1]
            
            
            
class Poisson_OU_asymptotic_backwards(Poisson_OU_asymptotic):
    def run_greaterThanMu(self):
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S+1+i] + 0.5 * self.delta_s
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atGreaterThanMu[-1]
            tmp = old_value + self.gamma * self.delta_s * old_value ** 2 - \
                np.true_divide(2 * self.M * self.delta_s, self.sigma_s ** 2) * (np.exp(-self.kappa * (self.s_long_term_mean - current_s)) + np.exp(self.kappa * (self.s_long_term_mean - current_s)))\
                + np.true_divide(2 * self.ergodic_const * self.delta_s, self.sigma_s ** 2)
            self.derivative_atGreaterThanMu.append(np.true_divide(tmp, 1 + self.delta_s * np.true_divide(2*self.alpha*(self.s_long_term_mean - current_s), self.sigma_s ** 2)))
        self.value_function[self.half_implement_S+1:] = self.valueFunction_at_mu + np.cumsum(self.delta_s * np.array(self.derivative_atGreaterThanMu))
    def run_lessThanMu(self):
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S - 1 - i] - 0.5 * self.delta_s
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atLessThanMu[-1]
            tmp = -old_value + self.gamma * self.delta_s * old_value ** 2 - \
                np.true_divide(2 * self.M * self.delta_s, self.sigma_s ** 2) * (np.exp(-self.kappa * (self.s_long_term_mean - current_s)) + np.exp(self.kappa * (self.s_long_term_mean - current_s)))\
                + np.true_divide(2 * self.ergodic_const * self.delta_s, self.sigma_s ** 2)

            self.derivative_atLessThanMu.append(np.true_divide(tmp, -1 + self.delta_s * np.true_divide(2*self.alpha*(self.s_long_term_mean - current_s), self.sigma_s ** 2)))
        
        self.value_function[:self.half_implement_S] = (self.valueFunction_at_mu - np.cumsum(self.delta_s * np.array(self.derivative_atLessThanMu)))[::-1]
    
class Poisson_OU_asymptotic_fullBackwards_average(Poisson_OU_asymptotic):
    def run_greaterThanMu(self):
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S+i] + 0.5 * self.delta_s
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atGreaterThanMu[-1]
            a = self.gamma * self.delta_s * 0.25
            b = 0.5 * (self.gamma * old_value * self.delta_s \
                       - np.true_divide(2 * self.alpha, self.sigma_s ** 2) * (self.s_long_term_mean - current_s) * self.delta_s - 2)
            c = old_value + self.delta_s * (0.25 * self.gamma * old_value ** 2 \
                                            - np.true_divide(2 * self.alpha, self.sigma_s ** 2) * 0.5 * (self.s_long_term_mean - current_s) * old_value\
                                            - np.true_divide(2 * self.M, self.sigma_s ** 2) *  (np.exp(-self.kappa * (self.s_long_term_mean - current_s)) + np.exp(self.kappa * (self.s_long_term_mean - current_s)))\
                                            + np.true_divide(2 * self.ergodic_const, self.sigma_s ** 2))
            if b**2 < 4*a*c:
                print current_s
                print a, b, c
                print old_value
                raise Exception("determinant is negative. Maybe try smaller delta_s.")
            
            self.derivative_atGreaterThanMu.append(np.true_divide( - b - np.sqrt(b**2 - 4*a*c), 2 * a))
        
        self.value_function[self.half_implement_S+1:] = self.valueFunction_at_mu + np.cumsum(self.delta_s * np.array(self.derivative_atGreaterThanMu))

    def run_lessThanMu(self):
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S  - i] - 0.5 * self.delta_s
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atLessThanMu[-1]
            a = self.gamma * self.delta_s * 0.25
            b = 0.5 * (self.gamma * old_value * self.delta_s \
                       - np.true_divide(2 * self.alpha, self.sigma_s ** 2) * (self.s_long_term_mean - current_s) * self.delta_s + 2)
            
            c = -old_value + self.delta_s * (0.25 * self.gamma * old_value ** 2 \
                                            - np.true_divide(2 * self.alpha, self.sigma_s ** 2) * 0.5 * (self.s_long_term_mean - current_s) * old_value\
                                            - np.true_divide(2 * self.M, self.sigma_s ** 2) *  (np.exp(-self.kappa * (self.s_long_term_mean - current_s)) + np.exp(self.kappa * (self.s_long_term_mean - current_s)))\
                                            + np.true_divide(2 * self.ergodic_const, self.sigma_s ** 2))
    
            if b**2 < 4*a*c:
                print current_s
                print a, b, c
                print old_value
                raise Exception("determinant is negative. Maybe try smaller delta_s. ")
            
            self.derivative_atLessThanMu.append(np.true_divide( - b + np.sqrt(b**2 - 4*a*c), 2 * a))

        self.value_function[:self.half_implement_S] = (self.valueFunction_at_mu - np.cumsum(self.delta_s * np.array(self.derivative_atLessThanMu)))[::-1]

class Poisson_OU_asymptotic_RK4(Poisson_OU_asymptotic):
    
    def run_greaterThanMu(self):
        
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S + i]
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atGreaterThanMu[-1]
            k1 = self.updatePart_helper(old_value, current_s)
            k2 = self.updatePart_helper(old_value + 0.5 * k1 * self.delta_s, current_s + 0.5 * self.delta_s)
            k3 = self.updatePart_helper(old_value + 0.5 * k2 * self.delta_s, current_s + 0.5 * self.delta_s)
            k4 = self.updatePart_helper(old_value + k3 * self.delta_s, current_s + self.delta_s)
            self.derivative_atGreaterThanMu.append( old_value + np.true_divide(1, 6) * (k1 + 2 * k2 + 2 * k3 + k4) * self.delta_s)
        self.value_function[self.half_implement_S+1:] = self.valueFunction_at_mu + np.cumsum(self.delta_s * np.array(self.derivative_atGreaterThanMu))
    def run_lessThanMu(self):
        for i in xrange(self.half_implement_S):
            current_s = self.implement_s_space[self.half_implement_S - i]
            old_value = self.derivative_at_mu if i == 0 else self.derivative_atLessThanMu[-1]
            k1 = self.updatePart_helper(old_value, current_s)
            k2 = self.updatePart_helper(old_value - 0.5 * k1 * self.delta_s, current_s - 0.5 * self.delta_s)
            k3 = self.updatePart_helper(old_value - 0.5 * k2 * self.delta_s, current_s - 0.5 * self.delta_s)
            k4 = self.updatePart_helper(old_value - k3 * self.delta_s, current_s - self.delta_s)
            self.derivative_atLessThanMu.append( old_value - np.true_divide(1, 6) * (k1 + 2 * k2 + 2 * k3 + k4) * self.delta_s)
        self.value_function[:self.half_implement_S] = (self.valueFunction_at_mu - np.cumsum(self.delta_s * np.array(self.derivative_atLessThanMu)))[::-1]



    
    
    
    
    
    
    
    
    
    
    
                
        
         
    
        