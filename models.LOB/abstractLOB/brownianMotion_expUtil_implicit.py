'''
Created on Jan 28, 2015

@author: weiluo
'''
from abstractDerivedLOB import AbstractImplicitLOB_NeumannBC
import numpy as np
class BrownianMotion_ExpUtil_Implicit_NeumannBC(AbstractImplicitLOB_NeumannBC):
    '''
    Basically "secondTry".
    '''


    @property
    def half_I(self):
        return self._half_I
    
    @property
    def q_space(self):
        return self._q_space
    
    @property
    def delta_q(self):
        return self._delta_q
    
    def compute_q_space(self):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).compute_q_space()
        self._q_space, self._delta_q  = np.linspace(-self.N, self.N, self.I, retstep = True)
        
    @property
    def a_control(self):
        return self._data_helper(self._index_a_control_2darray, self.extend_space -1, self.extend_space - 1)
    
    @property
    def b_control(self):
        return self._data_helper(self._index_b_control_2darray, self.extend_space - 1, self.extend_space - 1)
    
    def __init__(self, expIntensity_lower_tolerance_power = 2, tolerance_power = 2, *args, **kwargs):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).__init__(*args, **kwargs)
        self.expIntensity_lower_tolerance_power = expIntensity_lower_tolerance_power
        self.control_lower_bound = 0

        #A*exp(-kappa*control_upper_bound) < 10**(-expIntensity_lower_tolerance_power)
        self.control_upper_bound = (np.log(self.A) + np.log(10)*self.expIntensity_lower_tolerance_power)/self.kappa
        self.tolerance_power = tolerance_power

    def terminal_condition(self):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).terminal_condition()
        return -1 * np.ones(len(self.implement_q_space))
  
    def fixedData(self, v):
        delta_v_forward = np.true_divide(v[2:] - v[1:-1], self.delta_q)
        delta_v_backward = np.true_divide(v[1:-1] - v[:-2], self.delta_q)
        delta_v_second_order = np.true_divide(v[2:] - 2*v[1:-1] + v[:-2], self.delta_q ** 2)
        v_mid = v[1:-1]
        implement_q_space_mid = self.implement_q_space[1:-1]
        return [delta_v_forward, delta_v_backward, delta_v_second_order, v_mid, implement_q_space_mid]
    def F_a_noExp(self, delta_a, delta_v_forward, delta_v_backward, delta_v_second_order, v_mid, implement_q_space_mid):
        return 0.5 * self.gamma**2 * v_mid * delta_a**2 + (self.gamma * delta_v_forward - self.gamma * v_mid * (1 + self.beta * implement_q_space_mid)) * delta_a + (0.5 * delta_v_second_order - delta_v_backward)

    def F_b_noExp(self, delta_b, delta_v_forward, delta_v_backward, delta_v_second_order, v_mid, implement_q_space_mid):
        return 0.5 * self.gamma**2 * v_mid * delta_b**2 + (-self.gamma * delta_v_backward - self.gamma * v_mid * (1 - self.beta * implement_q_space_mid)) * delta_b + (0.5 * delta_v_second_order + delta_v_forward)

    
    def feedback_control(self, v):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).feedback_control(v)
        [delta_v_forward, delta_v_backward, delta_v_second_order, v_mid, implement_q_space_mid] = self.fixedData(v)
    
        Delta_a = self.kappa**2*(self.gamma * delta_v_forward - self.gamma * v[1:-1]*(1 + self.beta * implement_q_space_mid))**2 + self.gamma**4 * v[1:-1]**2\
         - 2 * self.kappa**2 * self.gamma**2 * v[1:-1] * (0.5 * delta_v_second_order - delta_v_backward)
        Delta_b = self.kappa**2*(self.gamma * delta_v_backward + self.gamma * v[1:-1]*(1 - self.beta * implement_q_space_mid))**2 + self.gamma**4 * v[1:-1]**2\
         - 2 * self.kappa**2 * self.gamma**2 * v[1:-1] * (0.5 * delta_v_second_order + delta_v_forward)
        
        optimal_a = np.zeros(self.implement_I - 2)
        optimal_b = np.zeros(self.implement_I - 2)


        optimal_a[Delta_a<0] = self.control_upper_bound
        optimal_b[Delta_b<0] = self.control_upper_bound

        index_positive_delta_a = Delta_a >= 0
        index_positive_delta_b = Delta_b >= 0
        
        optimal_a[index_positive_delta_a] = np.true_divide(1, self.kappa) + np.true_divide(1+self.beta * implement_q_space_mid[index_positive_delta_a], self.gamma)- np.true_divide(1, self.gamma)*np.true_divide(delta_v_forward, v_mid)[index_positive_delta_a] + np.true_divide(np.sqrt(Delta_a[index_positive_delta_a]), self.gamma**2 * self.kappa * v_mid[index_positive_delta_a])
       
        optimal_a[np.all([Delta_a >=0, optimal_a<0], axis=0)] = 0
       
        optimal_a[np.all([Delta_a >=0, self.F_a_noExp(optimal_a, delta_v_forward, delta_v_backward, delta_v_second_order\
                                                      , v_mid, implement_q_space_mid) <= -10**(-self.tolerance_power)], axis=0)] = self.control_upper_bound
       
        optimal_b[index_positive_delta_b] = np.true_divide(1, self.kappa) + np.true_divide(1 - self.beta * implement_q_space_mid[index_positive_delta_b], self.gamma) + np.true_divide(1, self.gamma)*np.true_divide(delta_v_backward, v_mid)[index_positive_delta_b] + np.true_divide(np.sqrt(Delta_b[index_positive_delta_b]), self.gamma**2 * self.kappa * v_mid[index_positive_delta_b])
        optimal_b[np.all([Delta_b >=0, optimal_b<0], axis=0)] = 0
        optimal_b[np.all([Delta_b >=0, self.F_b_noExp(optimal_b, delta_v_forward, delta_v_backward, delta_v_second_order\
                                                      , v_mid, implement_q_space_mid) <= -10**(-self.tolerance_power)], axis=0)] = self.control_upper_bound
      
        return [optimal_a, optimal_b]


    def coef_offset_helper(self, v_curr, curr_control):
        curr_a_star = curr_control[0]
        curr_b_star = curr_control[1]
        co_left = ( self.A * np.exp(- self.kappa * curr_a_star) + np.true_divide(self.A * np.exp(- self.kappa * curr_a_star)\
             + self.A * np.exp(- self.kappa * curr_b_star), 2*self.delta_q)+ self.gamma * self.A * np.exp(- self.kappa * curr_b_star) * curr_b_star) * self.delta_t
        co_right = ( self.A * np.exp(- self.kappa * curr_b_star) + np.true_divide(self.A * np.exp(- self.kappa * curr_a_star)\
             + self.A * np.exp(- self.kappa * curr_b_star), 2*self.delta_q)+ self.gamma *  self.A * np.exp(- self.kappa * curr_a_star) * curr_a_star) * self.delta_t
        return [co_left, co_right]

    
    def coef_at_minus_one(self, v_curr, curr_control):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).coef_at_minus_one(v_curr, curr_control)
        
        return -self.coef_at_minus_one_helper(self.coef_offset_helper(v_curr, curr_control)[0])
    
    def coef_at_plus_one(self, v_curr, curr_control):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).coef_at_plus_one(v_curr, curr_control)
        curr_a_star = curr_control[0]
        curr_b_star = curr_control[1]
        return -self.coef_at_plus_one_helper(self.coef_offset_helper(v_curr, curr_control)[1])
        
    def coef_at_curr(self, v_curr, curr_control):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).coef_at_curr(v_curr, curr_control)
        curr_a_star = curr_control[0]
        curr_b_star = curr_control[1]
        tmp = self.delta_q + self.delta_t * self.delta_q *( self.gamma *(curr_a_star * self.A * np.exp(- self.kappa * curr_a_star) \
            + curr_b_star * self.A * np.exp(- self.kappa * curr_b_star) ) + self.gamma*self.implement_q_space[1:-1] * self.beta * \
            (curr_a_star * self.A * np.exp(- self.kappa * curr_a_star) - curr_b_star * self.A * np.exp(- self.kappa * curr_b_star))\
              - 0.5 * self.gamma**2 * ( self.A * np.exp(- self.kappa * curr_a_star)  *curr_a_star * curr_a_star + self.A * np.exp(- self.kappa * curr_b_star) \
              *curr_b_star*curr_b_star) - 0.5*self.gamma**2*self.implement_q_space[1:-1]*self.implement_q_space[1:-1]*self.sigma_s**2 )
        co_left, co_right = self.coef_offset_helper(v_curr, curr_control)
        tmp -= (co_left + co_right)
        return self.coef_at_curr_helper(tmp)


    
    def equation_right(self, v_curr, curr_control):
        super(BrownianMotion_ExpUtil_Implicit_NeumannBC, self).equation_right(v_curr, curr_control)
        return self.equation_right_helper(v_curr[1:-1])





    