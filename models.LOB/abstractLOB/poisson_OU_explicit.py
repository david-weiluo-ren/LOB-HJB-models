'''
Created on Mar 2, 2015

@author: weiluo
'''

import abc
from abstract_OU_LOB import Abstract_OU_LOB
import numpy as np
from pylab import plot, show
class Poisson_explicit_OU_LOB(Abstract_OU_LOB):
    '''
    classdocs
    '''
    
   
    
        
    
    
    def terminal_condition_real(self):
        return np.outer(self.implement_s_space, self.implement_q_space)
    
    def feedback_control(self, v):
        v_s = np.true_divide(v[2:,:] - v[:-2, :], 2*self.delta_s)[:, 1:-1]
        
        v_q_backward = (v[:, :-2] - v[:, 1:-1])[1:-1, :]
        v_q_forward = (v[:,2:] - v[:, 1:-1])[1:-1, :]
        s_space_casted = np.outer(self.implement_s_space[1:-1], np.ones(len(self.implement_q_space) - 2)) 
        LARGE_NUM = 10
        
        
        v_s_reshape = np.reshape(v_s, (-1,))
        optimal_a_reshape = np.ones(len(v_s_reshape)) * LARGE_NUM
        optimal_b_reshape = np.ones(len(v_s_reshape)) * LARGE_NUM
        v_q_backward_reshape = np.reshape(v_q_backward, (-1,))
        v_q_forward_reshape = np.reshape(v_q_forward, (-1,))
        s_space_casted_reshape = np.reshape(s_space_casted, (-1,))

       
        #index_a = np.all([(1+self.gamma*self.beta*v_s_reshape>0) ], axis=0)
        index_a = 1+self.gamma*self.beta*v_s_reshape>0
        index_b = 1-self.gamma*self.beta*v_s_reshape>0

        try:    
            optimal_a_reshape[index_a] = np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)) - np.log(1+self.gamma*self.beta*v_s_reshape[index_a]))\
                    -(v_q_backward_reshape[index_a] + s_space_casted_reshape[index_a])
            optimal_b_reshape[index_b] = np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)) - np.log(1-self.gamma*self.beta*v_s_reshape[index_b ]))\
                    -(v_q_forward_reshape[index_b ] - s_space_casted_reshape[index_b ])
            return [np.reshape(optimal_a_reshape,np.shape(v_s)), np.reshape(optimal_b_reshape,np.shape(v_s))]
        except Exception as e:
            print e
            plot(v)
            show()
            
            raise Exception()
      
    
    def truncate_control_at_zero(self, arrays):
        return [self.truncate_at_zero(arr) for arr in arrays]
    
    def one_step_back(self, v_curr, step_index=None):
        v_new = np.ndarray((len(self.implement_s_space), len(self.implement_q_space)))
        
        optimal_a, optimal_b = self.truncate_control_at_zero (self.feedback_control(v_curr))
        
        v_ss = np.true_divide(v_curr[2:,:] + v_curr[:-2, :] - 2* v_curr[1:-1, :], self.delta_s**2)[:, 1:-1]
        v_s = np.true_divide(v_curr[2:,:] - v_curr[:-2, :], 2*self.delta_s)[:, 1:-1]
        s_space_casted = np.outer(self.implement_s_space[1:-1], np.ones(len(self.implement_q_space) - 2)) 
        '''
        v_new[1:-1, 1:-1] = v_curr[1:-1, 1:-1] + self.delta_t*((v_ss - self.gamma * v_s*v_s)*self.sigma_s**2*0.5\
                                + v_s * (self.alpha*(self.s_0 - s_space_casted) + self.beta*self.A*(np.exp(-self.kappa*optimal_a) - np.exp(-self.kappa*optimal_b)) )\
                                 + np.true_divide(self.A, self.gamma) * np.exp(-self.kappa*optimal_a)*(1 - np.true_divide(self.kappa, self.kappa+self.gamma)*(1+self.kappa*self.beta*v_s) )\
                                + np.true_divide(self.A, self.gamma) * np.exp(-self.kappa*optimal_b)*(1 - np.true_divide(self.kappa, self.kappa+self.gamma)*(1-self.kappa*self.beta*v_s) )\
                                 )
        '''
        np.seterr(all="raise")
        try:
            v_new[1:-1, 1:-1] = v_curr[1:-1, 1:-1] +self.delta_t*((v_ss - self.gamma * v_s*v_s)*self.sigma_s**2*0.5\
                                + v_s * self.alpha*(self.s_long_term_mean-s_space_casted)+\
                                np.true_divide(self.A, self.kappa+self.gamma)*(np.exp(-self.kappa*optimal_a)*(1+self.beta*self.gamma*v_s) + \
                                np.exp(-self.kappa*optimal_b)*(1-self.beta*self.gamma*v_s)  )                                             )
        
            v_new[0,1:-1] = v_new[1, 1:-1]
            v_new[-1, 1:-1] = v_new[-2, 1:-1]
            v_new[:, 0] = v_new[:, 1]
            v_new[:, -1] = v_new[:, -2]
            self._a_control.append(optimal_a)
            self._b_control.append(optimal_b)
            return v_new
        except Exception as e:
            print e
            for i in xrange(np.shape(optimal_a)[1]):
                plot(optimal_a[:, i])
           
            show()

    
    
    
    
    
    
    
    
    
class Poisson_explicit_OU_LOB_noTruncation(Poisson_explicit_OU_LOB):
    
    def truncate_control_at_zero(self, arrays):
        return arrays
    
    
    
class Poisson_explicit_OU_LOB_sameslope(Poisson_explicit_OU_LOB):
    def one_step_back(self, v_curr, step_index=None):
        v_new = np.ndarray((len(self.implement_s_space), len(self.implement_q_space)))
        
        optimal_a, optimal_b = self.truncate_control_at_zero (self.feedback_control(v_curr))
        
        v_ss = np.true_divide(v_curr[2:,:] + v_curr[:-2, :] - 2* v_curr[1:-1, :], self.delta_s**2)[:, 1:-1]
        v_s = np.true_divide(v_curr[2:,:] - v_curr[:-2, :], 2*self.delta_s)[:, 1:-1]
        s_space_casted = np.outer(self.implement_s_space[1:-1], np.ones(len(self.implement_q_space) - 2)) 
        '''
        v_new[1:-1, 1:-1] = v_curr[1:-1, 1:-1] + self.delta_t*((v_ss - self.gamma * v_s*v_s)*self.sigma_s**2*0.5\
                                + v_s * (self.alpha*(self.s_0 - s_space_casted) + self.beta*self.A*(np.exp(-self.kappa*optimal_a) - np.exp(-self.kappa*optimal_b)) )\
                                 + np.true_divide(self.A, self.gamma) * np.exp(-self.kappa*optimal_a)*(1 - np.true_divide(self.kappa, self.kappa+self.gamma)*(1+self.kappa*self.beta*v_s) )\
                                + np.true_divide(self.A, self.gamma) * np.exp(-self.kappa*optimal_b)*(1 - np.true_divide(self.kappa, self.kappa+self.gamma)*(1-self.kappa*self.beta*v_s) )\
                                 )
        '''
        np.seterr(all="raise")
        try:
            v_new[1:-1, 1:-1] = v_curr[1:-1, 1:-1] +self.delta_t*((v_ss - self.gamma * v_s*v_s)*self.sigma_s**2*0.5\
                                + v_s * self.alpha*(self.s_long_term_mean-s_space_casted)+\
                                np.true_divide(self.A, self.kappa+self.gamma)*(np.exp(-self.kappa*optimal_a)*(1+self.beta*self.gamma*v_s) + \
                                np.exp(-self.kappa*optimal_b)*(1-self.beta*self.gamma*v_s)  )                                             )
        
            v_new[0,1:-1] = 2*v_new[1, 1:-1] - v_new[2, 1:-1]
            v_new[-1, 1:-1] = 2*v_new[-2, 1:-1] - v_new[-3, 1:-1]
            v_new[:, 0] = 2 * v_new[:, 1] - v_new[:, 2]
            v_new[:, -1] = 2 * v_new[:, -2] - v_new[:, -3]
            self._a_control.append(optimal_a)
            self._b_control.append(optimal_b)
            return v_new
        except Exception as e:
            print e
            for i in xrange(np.shape(optimal_a)[1]):
                plot(optimal_a[:, i])
           
            show()

   
    
    