'''
Created on Jan 27, 2015

@author: weiluo
'''
import numpy as np
from _pyio import __metaclass__
import abc
from scipy.sparse.linalg import spsolve
from scipy import sparse
from abc import abstractmethod
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
    
   
    @property
    def result(self):
        return self._data_helper(self._index_result_2darray)
    
    @property
    def a_control(self):
        return self._data_helper(self._index_a_control_2darray)
    
    @property
    def b_control(self):
        return self._data_helper(self._index_b_control_2darray)
    
    def _data_helper(self, data_index, front_extend_space=None, behind_extend_space=None):
        if len(self._data[data_index]) == 0:
            return self._data[data_index]
        else:
            return self.user_friendly_list_of_array(self._data[data_index], data_index, front_extend_space, behind_extend_space)

        
    @property
    def _result(self):
        return self._data[self._index_result_2darray]
    @property
    def _a_control(self):
        return self._data[self._index_a_control_2darray]
    
    @property
    def _b_control(self):
        return self._data[self._index_b_control_2darray]
    
    @property
    def extend_space(self):
        return self._extend_space
       
    def __init__(self, gamma = 0.1, A = 1, kappa = 0.1, beta = 0.02, N = 20, half_I = 2000, sigma_s = 0.05,\
                 q_0 = 0, x_0 = 3.0, s_0 = 5.0,  delta_t = 0.01, verbose = False, num_time_step = 100, extend_space = 2,\
                 expIntensity_lower_tolerance_power = 2):
        self.gamma = gamma
        self.A = A
        self.kappa = kappa
        self.beta = beta
        self.N = N
        self.sigma_s = sigma_s
        self._half_I = half_I
        self.delta_t = delta_t
        self.I = 2 * self.half_I + 1
        self._extend_space = extend_space
        self.verbose = verbose
        self._delta_q = None
        self._q_space = None
        self.compute_q_space()
        self.implement_I = self.I + 2 * self.extend_space
        self.implement_q_space = np.hstack((-np.arange(self.extend_space, 0, -1) * self.delta_q + self.q_space[0],\
                                             self.q_space, \
                                             np.arange(1, self.extend_space+1) * self.delta_q + self.q_space[-1]))
        
        self.num_time_step = num_time_step
        self.T = self.delta_t * self.num_time_step
        self.step_index = 0
        self.v_init = self.terminal_condition()
        
        self._data = [[], [], []]
        
        self._cache = [None, None, None]
        self._index_result_2darray = 0
        self._index_a_control_2darray = 1
        self._index_b_control_2darray = 2
        self.expIntensity_lower_tolerance_power = expIntensity_lower_tolerance_power
        self.control_lower_bound = 0

        #A*exp(-kappa*control_upper_bound) < 10**(-expIntensity_lower_tolerance_power)
        self.control_upper_bound = (np.log(self.A) + np.log(10)*self.expIntensity_lower_tolerance_power)/self.kappa

    def get_extend_space(self):
        return self.__extend_space


    def set_extend_space(self, value):
        self.__extend_space = value


    def del_extend_space(self):
        del self.__extend_space

        
    @abstractmethod
    def compute_q_space(self):
        pass
    
    def check_user_friendly_helper(self, front_extend_space=None, behind_extend_space=None, totalLen = None):
        if front_extend_space is None:
            front_extend_space = self.extend_space
        if behind_extend_space is None:
            behind_extend_space = self.extend_space
        if totalLen is None:
            totalLen = self.I
        return [front_extend_space, behind_extend_space, totalLen]
    def user_friendly_index(self, front_extend_space=None, behind_extend_space=None, totalLen = None):
        front_extend_space, behind_extend_space, totalLen = self.check_user_friendly_helper(front_extend_space, behind_extend_space, totalLen)
        return np.arange(0, totalLen)[front_extend_space: -behind_extend_space]
    
    def user_friendly_array(self, arr, front_extend_space, behind_extend_space):
        return arr[self.user_friendly_index(front_extend_space, behind_extend_space, len(arr))] 
         
    def user_friendly_list_of_array(self, list_of_array, cache_index=None, front_extend_space=None, behind_extend_space=None, call_func=None):
        friendly_index = self.user_friendly_index(front_extend_space, behind_extend_space, len(list_of_array[0]))
        if cache_index is not None:
            if self._cache[cache_index] is None:
                self._cache[cache_index] = np.vstack(list_of_array)
           
            elif self._cache[cache_index].shape[0] < len(list_of_array):
                self._cache[cache_index] = np.vstack((self._cache[cache_index],\
                                                     list_of_array[self._cache[cache_index].shape[0]:]))
            if call_func is not None:
                return call_func(self._cache[cache_index]).T[friendly_index].T
            else:
                return self._cache[cache_index].T[friendly_index].T
        else:
            tmp_2darray = np.vstack(list_of_array)
            if call_func is not None:
                tmp_2darray = call_func(tmp_2darray)
            return tmp_2darray.T[friendly_index].T
    @abc.abstractmethod
    def terminal_condition(self):
        pass
    
    @abc.abstractmethod
    def one_step_back(self, v_curr):
        """
        To compute the value function and optimal feedback controls.
        """
        pass
    
    def run(self, K = None, use_cache = False):
        if not use_cache:
            self._data = [[], [], []]  #If no cache should be used, everything is restarted, all previous data are lost.
        
        if K is None:
            K = self.num_time_step
        if use_cache:
            K = K - len(self._result) if K > len(self._result) else 0
        v_curr = self.v_init if not use_cache or len(self.result) == 0 else self._result.pop()
           
        for i in xrange(K):
            self._result.append(v_curr)
           
            v_curr = self.one_step_back(v_curr)
            self.step_index += 1
             
class AbstractImplicitLOB(AbstractLOB):
    
    __metaclass__ = abc.ABCMeta
    
    
    
    
    def __init__(self, iter_max = 200, new_weight = 0.1, abs_threshold = 10**(-4), rlt_threshold = 10**(-2),\
                 *args, **kwargs):
        super(AbstractImplicitLOB, self).__init__(*args, **kwargs)
        self.iter_max = iter_max
        self.new_weight = new_weight
        self.abs_threshold = abs_threshold
        self.rlt_threshold = rlt_threshold
        self.index_for_debug = 0
    
    def one_step_back(self, v_curr):
        self.index_for_debug += 1
        AbstractLOB.one_step_back(self, v_curr)
        v_tmp = v_curr
        iter_count = 0
        while True:
            curr_control = self.feedback_control(v_tmp)
            v_new = self.one_iteration(v_curr, curr_control)
            if self.close_enough(v_new, v_tmp):
                if self.verbose:
                    print "the {}th iteration converges in {} iterations".format(self.step_index, iter_count),
                return_control= self.feedback_control(v_new)
                self._a_control.append(return_control[0])
                self._b_control.append(return_control[1])
                return v_new
            v_tmp = self.new_weight * v_new + (1 - self.new_weight) * v_tmp
            
            iter_count += 1
            if iter_count > self.iter_max:
                raise Exception('iteration cannot converge!')
    
    @abc.abstractmethod
    def feedback_control(self, v):
        """
        should return [optimal_a, optimal_b]
        """
        pass
    
    def close_enough(self, v_new, v_curr):
        return np.allclose(v_curr, v_new, self.rlt_threshold, self.abs_threshold)\
            and np.allclose(v_new, v_curr, self.rlt_threshold, self.abs_threshold)
        
    def close_enough_old(self, v_new, v_curr):
        zero_index = np.abs(v_curr) < self.abs_threshold
        non_zero_index = ~zero_index
        
        if np.in1d(True, zero_index) and \
            np.max(np.abs(v_new[zero_index])) > 2 * self.abs_threshold:
            return False
        
        return np.max(np.abs((v_new[non_zero_index] - v_curr[non_zero_index])/v_curr[non_zero_index]))\
                < self.rlt_threshold
    
    
    def one_iteration(self, v_curr, curr_control):
       
        
        
        eq_right, co_matrix = self.linear_system(v_curr, curr_control)
        x = spsolve(co_matrix, eq_right)
        return spsolve(co_matrix, eq_right)
    
    @abstractmethod
    def linear_system(self, v_curr, curr_control):
        """
        should return [eq_right, co_matrix]
        
        eq_right should be an array of length implement_I.
        
        co_matrix should be a implement_I times implement_I sparse matrix.
        
        """  
        pass