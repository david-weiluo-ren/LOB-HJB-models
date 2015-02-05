'''
Created on Feb 2, 2015

@author: weiluo
'''
from abstractLOB import AbstractImplicitLOB
import numpy as np
class BrownianMotion_ExpUtil_Implicit(AbstractImplicitLOB):
    '''
    Basically "secondTry".
    Default is using NeumannBC
    '''

    @property
    def extend_space(self):
        return self._extend_space * int(self.half_I / self.N)
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
        super(BrownianMotion_ExpUtil_Implicit, self).compute_q_space()
        self._q_space, self._delta_q  = np.linspace(-self.N, self.N, self.I, retstep = True)
        
    @property
    def a_control(self):
        return self._data_helper(self._index_a_control_2darray, self.extend_space -1, self.extend_space - 1)
    
    @property
    def b_control(self):
        return self._data_helper(self._index_b_control_2darray, self.extend_space - 1, self.extend_space - 1)
    
    def __init__(self, BC_class, linear_system_helper, tolerance_power = 2, *args, **kwargs):
        super(BrownianMotion_ExpUtil_Implicit, self).__init__(*args, **kwargs)
        
        self.BC = BC_class(linear_system_helper, self.implement_I)

        self.tolerance_power = tolerance_power
        
    def terminal_condition(self):
        super(BrownianMotion_ExpUtil_Implicit, self).terminal_condition()
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
        
        super(BrownianMotion_ExpUtil_Implicit, self).feedback_control(v)
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


    def linear_system_matrix_helper(self, v_curr, curr_control):
        curr_a_star = curr_control[0]
        curr_b_star = curr_control[1]
        co_left = ( self.A * np.exp(- self.kappa * curr_a_star) + np.true_divide(self.A * np.exp(- self.kappa * curr_a_star)\
             + self.A * np.exp(- self.kappa * curr_b_star), 2*self.delta_q)+ self.gamma * self.A * np.exp(- self.kappa * curr_b_star) * curr_b_star) * self.delta_t
        co_right = ( self.A * np.exp(- self.kappa * curr_b_star) + np.true_divide(self.A * np.exp(- self.kappa * curr_a_star)\
             + self.A * np.exp(- self.kappa * curr_b_star), 2*self.delta_q)+ self.gamma *  self.A * np.exp(- self.kappa * curr_a_star) * curr_a_star) * self.delta_t
        co_mid = co_left + co_right + self.delta_q + self.delta_t * self.delta_q *( self.gamma *(curr_a_star * self.A * np.exp(- self.kappa * curr_a_star) \
            + curr_b_star * self.A * np.exp(- self.kappa * curr_b_star) ) + self.gamma*self.implement_q_space[1:-1] * self.beta * \
            (curr_a_star * self.A * np.exp(- self.kappa * curr_a_star) - curr_b_star * self.A * np.exp(- self.kappa * curr_b_star))\
              - 0.5 * self.gamma**2 * ( self.A * np.exp(- self.kappa * curr_a_star)  *curr_a_star * curr_a_star + self.A * np.exp(- self.kappa * curr_b_star) \
              *curr_b_star*curr_b_star) - 0.5*self.gamma**2*self.implement_q_space[1:-1]*self.implement_q_space[1:-1]*self.sigma_s**2 )
        
       
        
        return [-self.BC.coef_at_minus_one_helper(co_left),  -self.BC.coef_at_plus_one_helper(co_right), self.BC.coef_at_curr_helper(co_mid)]
        
            
    def equation_right(self, v_curr, curr_control):
        
        return self.BC.equation_right_helper(v_curr[1:-1]*self.delta_q)
    
    def linear_system_helper(self, v_curr, curr_control):

        matrix_data = self.linear_system_matrix_helper(v_curr, curr_control)
        matrix_data.append(self.equation_right(v_curr, curr_control))
        return matrix_data
    
    
      
    
    def simulate_one_step_forward(self, index):
        super(BrownianMotion_ExpUtil_Implicit, self).simulate_one_step_forward(index)
        curr_control_a, curr_control_b = self.control_at_current_point(index, self.q[-1])
        drift_q_a = self.A * np.exp(-self.kappa * curr_control_a)
        drift_q_b = self.A * np.exp(-self.kappa * curr_control_b)
        delta_B_a =  np.sqrt(drift_q_a) * np.sqrt(self.delta_t)*np.random.normal(0,1,1)
        delta_B_b =  np.sqrt(drift_q_b) * np.sqrt(self.delta_t)*np.random.normal(0,1,1)

        delta_q_a  = drift_q_a * self.delta_t + delta_B_a
        delta_q_b  = drift_q_b * self.delta_t + delta_B_b

        delta_q = delta_q_b - delta_q_a
        delta_x = (self.s[-1] + curr_control_a) * delta_q_a - (self.s[-1] + curr_control_b) * delta_q_b
        delta_s = self.beta * self.A * (np.exp(-self.kappa * curr_control_a) * curr_control_a -np.exp(-self.kappa * curr_control_b) * curr_control_b ) * self.delta_t\
         + self.sigma_s * np.sqrt(self.delta_t) * np.random.normal(0,1,1)
        self.q.append(self.q[-1] + delta_q)
        self.q_a.append(self.q_a[-1] + delta_q_a)
        self.q_b.append(self.q_b[-1] + delta_q_b)

        self.x.append(self.x[-1] + delta_x)        
        self.s.append(self.s[-1] + delta_s)
        self.simulate_control_a.append(curr_control_a)
        self.simulate_control_b.append(curr_control_b)
    
    
    
    
    
    
    
    