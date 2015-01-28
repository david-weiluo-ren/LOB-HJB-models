'''
Created on Jan 27, 2015

@author: weiluo
'''
import numpy as np
from _pyio import __metaclass__
import abc
from scipy.sparse.linalg import spsolve
from scipy import sparse
class AbstractLOB(object):
    '''
    An abstract model as parents for all LOB models, Brownian Motion or Poisson, exp utility or linear with penalty, with
    price impact or not.
    
    The user are interested in results in [-half_I, half_I], while in the implementation, the boundary point may be bad. So we 
    introduce "implementation_I" which extends the user_defined I a little. 
    '''

    __metaclass__ = abc.ABCMeta
    
    @abc.abstractproperty
    def q_space(self):
        pass
    
    @abc.abstractproperty
    def delta_q(self):
        pass
    
    @abc.abstractproperty
    def half_I(self):
        pass
   
    @abc.abstractproperty
    def resutlt(self):
        pass
    
    @abc.abstractproperty
    def a_control(self):
        pass
    
    @abc.abstractproperty
    def b_control(self):
        pass
    
    def __init__(self, gamma = 0.1, A = 1, kappa = 0.1, beta = 0.02, N = 20, halfI = 2000, sigma_s = 0.05,\
                 q_0 = 0, x_0 = 3.0, s_0 = 5.0,  delta_t = 0.01, verbose = False, num_time_step = 100):
        self.gamma = gamma
        self.A = A
        self.kappa = kappa
        self.beta = beta
        self.N = N
        self.sigma_s = sigma_s
       
        self.delta_t = delta_t
        self.I = 2 * self.half_I + 1
        self.num_time_step = num_time_step
        self.T = self.delta_t * self.num_time_step
        self.step_index = 0
        self.v_init = self.terminal_condition()
        
        self._result = []
        self._a_control = []
        self._b_control = []
    
        
    @abc.abstractmethod
    def terminal_condition(self):
        pass
    
    @abc.abstractmethod
    def one_step_back(self, v_curr):
        pass
    
    def run(self, K = None):
        if K is None:
            K = self.num_time_step

        K = K - len(self._result) if K > len(self._result) else 0
        v_curr = self.v_init if len(self.result) == 0 else self._result.pop()
           
        for i in xrange(K):
            self._result.append(v_curr)
            v_curr = self.oneStepBack(v_curr)
            self.step_index += 1
            
class AbstractImplicitLOB(AbstractLOB):
    
    __metaclass__ = abc.ABCMeta
    
    @abc.abstractproperty
    def implement_I(self):
        pass
    
    
    def __init__(self, gamma = 0.1, A = 1, kappa = 0.1, beta = 0.02, N = 20, halfI = 2000, sigma_s = 0.05,\
                 q_0 = 0, x_0 = 3.0, s_0 = 5.0,  delta_t = 0.01, verbose = False, num_time_step = 100,\
                 iter_max = 200, new_weight = 0.1, abs_threshold = 10**(-4), rlt_threshold = 10**(-2)):
        super.__init__(gamma = gamma, A = A, kappa = kappa, beta = beta, N = N, halfI = halfI, sigma_s = sigma_s,\
                       q_0 = q_0, x_0 = x_0, s_0 = s_0, delta_t = delta_t, verbose = verbose, num_time_step = num_time_step)
        self.iter_max = iter_max
        self.new_weight = new_weight
        self.abs_threshold = abs_threshold
        self.rlt_threshold = rlt_threshold
    
    def one_step_back(self, v_curr):
        v_tmp = v_curr
        iter_count = 0
        while True:
            curr_control = self.feedback_control(v_tmp)
            v_new = self.one_iteration(v_curr, curr_control)
            if self.close_enough(v_new, v_curr):
                if self.verbose:
                    print "the {}th iteration converges in {} iterations".format(self.step_index, iter_count),
                return_control= self.feedback_control(v_new)
                self._a_control.append(return_control[0])
                self._b_control.append(return_control[1])
                return v_new
            v_tmp = self.new_weight * v_new + (1 - self.new_weight) * v_curr
            iter_count += 1
            if iter_count > self.iter_max:
                raise Exception('iteration cannot converge!')
    
    @abc.abstractmethod
    def feedback_control(self, v):
        pass
    
    def close_enough(self, v_new, v_curr):
        zero_index = np.abs(v_curr) < self.abs_threshold
        non_zero_index = ~zero_index
        
        if np.in1d(True, zero_index) and \
            np.max(np.abs(v_new[zero_index])) > 2 * self.abs_threshold:
            return False
        
        return np.max(np.abs((v_new[non_zero_index] - v_curr[non_zero_index])/v_curr[non_zero_index]))\
                < self.rlt_threshold
    
    
    def oneIter(self, v_curr, curr_control):
        co_left = self.coef_at_minus_one(v_curr, curr_control)
        co_right = self.coef_at_plus_one(v_curr, curr_control)
        co_mid = self.coef_at_curr(v_curr, curr_control)
        eq_right= self.equation_right(v_curr, curr_control)
        
        data = [co_left, co_mid, co_right]   #mind the sign here.
        diags = [-1, 0, 1]
        co_matrix = sparse.spdiags(data, diags, self.implement_I, self.implement_I, format = 'csc')
        return spsolve(co_matrix, eq_right)
        
    @abc.abstractmethod
    def coef_at_minus_one(self, v_curr, curr_control):
        pass
        
    @abc.abstractmethod
    def coef_at_plus_one(self, v_curr, curr_control):
        pass
        
    @abc.abstractmethod
    def coef_at_curr(self, v_curr, curr_control):
        pass
    
    @abc.abstractmethod
    def equation_right(self, v_curr, curr_control):
        
    
            
        