'''
Created on Feb 2, 2015

@author: weiluo
'''
import numpy as np
from scipy import sparse

class ImplicitLOB_NeumannBC(object):
    '''
    helper functions almost the same as AbstractImplicitLOB_NeumannBC in 
    abstractDerivedLOB.py
    '''


    def __init__(self, linear_system_helper, implement_I):
        self.linear_system_helper = linear_system_helper
        self.implement_I = implement_I
    @classmethod
    def coef_at_minus_one_helper(cls, short_arr, positive=True):
        """
        short_arr should be the coef at minus for 
        space = implement_q_space[1],...,implement_q_space[-2]
        
        It will return an array with length implement_I. The 
        return array will append the input two 1 at the end.
        """
        append_sign = 1 if positive else -1
        return np.append(short_arr, np.array([1,1]) * append_sign)
    @classmethod
    def coef_at_plus_one_helper(cls, short_arr, positive=True):
        append_sign = 1 if positive else -1
        return np.hstack((np.array([1,1]) * append_sign, short_arr))
    @classmethod
    def coef_at_curr_helper(cls, short_arr, positive=True):
        append_sign = 1 if positive else -1
        return np.hstack((append_sign, short_arr, append_sign))
    @classmethod
    def equation_right_helper(cls, short_arr):
        return np.hstack((0, short_arr, 0))
    def linear_system(self, v_curr, curr_control):
        co_left, co_right, co_mid, eq_right = self.linear_system_helper(v_curr, curr_control)
        data = [co_left, co_mid, co_right]   #mind the sign here.
        diags = [-1, 0, 1]
        co_matrix = sparse.spdiags(data, diags, self.implement_I, self.implement_I, format = 'csc')
       
        return [eq_right, co_matrix]
 
 

class ImplicitLOB_sameSlopeBC(object):
    '''
    The boundary condition for HJB equation becomes 
    V_1 - V_0 = V_2 - V_1 and V_T - V_(T-1) = V_(T-1) - V_(T-2)
    which makes the co_matrix not a triangle matrix anymore    
    '''
   

    def __init__(self, linear_system_helper, implement_I):
        self.linear_system_helper = linear_system_helper
        self.implement_I = implement_I
    @classmethod
    def coef_at_minus_one_helper(cls, short_arr, positive=True):
        """
        short_arr should be the coef at minus for 
        space = implement_q_space[1],...,implement_q_space[-2]
        
        It will return an array with length implement_I. The 
        return array will append the input two 1 at the end.
        """
        append_sign = 1 if positive else -1
        return np.append(short_arr, np.array([2,2]) * append_sign)
    @classmethod
    def coef_at_plus_one_helper(cls, short_arr, positive=True):
        append_sign = 1 if positive else -1
        return np.hstack((np.array([2,2]) * append_sign, short_arr))
    @classmethod
    def coef_at_curr_helper(cls, short_arr, positive=True):
        append_sign = 1 if positive else -1
        return np.hstack((append_sign, short_arr, append_sign))
    @classmethod
    def equation_right_helper(cls, short_arr):
        return np.hstack((0, short_arr, 0))
    
    
   
   
    def linear_system(self, v_curr, curr_control):
        co_left, co_right, co_mid, eq_right = self.linear_system_helper(v_curr, curr_control)
       
        minus2_diag = np.zeros(self.implement_I)
        
        minus2_diag[-3] = 1
        plus2_diag = np.zeros(self.implement_I)
        plus2_diag[2] = 1        
        data = [minus2_diag, co_left, co_mid, co_right, plus2_diag]   #mind the sign here.
        diags = [-2, -1, 0, 1, 2]
        co_matrix = sparse.spdiags(data, diags, self.implement_I, self.implement_I)
       
        return [eq_right, co_matrix]
     
    
   
   
   
   
     