'''
Created on Jan 28, 2015

@author: weiluo
'''
import numpy as np
from abstractLOB import AbstractImplicitLOB
import abc
from _pyio import __metaclass__
class AbstractImplicitLOB_NeumannBC(AbstractImplicitLOB):
    '''
    Neumann Boundary Condition.
    Provide helper function for coef_at_minus_one,
    coef_at_plus_one, coef_at_curr, equation_right
    '''
    __metaclass__ = abc.ABCMeta

    def __init__(self, *args, **kwargs):
        super(AbstractImplicitLOB_NeumannBC, self).__init__(*args, **kwargs)
        
    def coef_at_minus_one_helper(self, short_arr, positive=True):
        """
        short_arr should be the coef at minus for 
        space = implement_q_space[1],...,implement_q_space[-2]
        
        It will return an array with length implement_I. The 
        return array will append the input two 1 at the end.
        """
        append_sign = 1 if positive else -1
        return np.append(short_arr, np.array([1,1]) * append_sign)
    
    def coef_at_plus_one_helper(self, short_arr, positive=True):
        append_sign = 1 if positive else -1
        return np.hstack((np.array([1,1]) * append_sign, short_arr))
    
    def coef_at_curr_helper(self, short_arr, positive=True):
        append_sign = 1 if positive else -1
        return np.hstack((append_sign, short_arr, append_sign))
    
    def equation_right_helper(self, short_arr):
        return np.hstack((0, short_arr, 0))
    
    
    
    
    