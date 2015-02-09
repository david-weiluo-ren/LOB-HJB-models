'''
Created on Feb 2, 2015

@author: weiluo
'''
from abstractLOB import AbstractImplicitLOB
import numpy as np
from scipy.special import lambertw
from scipy.optimize import newton, fmin_tnc
import scipy as sp
import math
from pylab import plot, show
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
        delta_B_a =  np.sqrt(drift_q_a) * np.sqrt(self.delta_t)*np.random.normal(0,1,1)[0]
        delta_B_b =  np.sqrt(drift_q_b) * np.sqrt(self.delta_t)*np.random.normal(0,1,1)[0]

        delta_q_a  = drift_q_a * self.delta_t + delta_B_a
        delta_q_b  = drift_q_b * self.delta_t + delta_B_b

        delta_q = delta_q_b - delta_q_a
        delta_x = (self.s[-1] + curr_control_a) * delta_q_a - (self.s[-1] + curr_control_b) * delta_q_b
        delta_s = self.beta * self.A * (np.exp(-self.kappa * curr_control_a) * curr_control_a -np.exp(-self.kappa * curr_control_b) * curr_control_b ) * self.delta_t\
         + self.sigma_s * np.sqrt(self.delta_t) * np.random.normal(0,1,1)[0]
        self.q.append(self.q[-1] + delta_q)
        self.q_a.append(self.q_a[-1] + delta_q_a)
        self.q_b.append(self.q_b[-1] + delta_q_b)

        self.x.append(self.x[-1] + delta_x)        
        self.s.append(self.s[-1] + delta_s)
        self.simulate_control_a.append(curr_control_a)
        self.simulate_control_b.append(curr_control_b)
    
    
    
    
class Poisson_expUtil_implicit(AbstractImplicitLOB):
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
        super(Poisson_expUtil_implicit, self).terminal_condition()
        return np.ones(len(self.implement_q_space))

    def F_1(self, x, beta1, beta2, beta3):

        return 0 if x == self.control_upper_bound else self.A*np.exp(-self.kappa*x)*(-beta1 * x + beta2 * np.exp(-self.gamma * x) - beta3)

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
                
            def function_for_root_helper(beta_1, beta_2, beta_3, x):
                return -(self.gamma + self.kappa)*beta_2 * np.exp(-self.gamma*x) + self.kappa * beta_1 * x - beta_1 + self.kappa * beta_3

            def function_for_root_derivative_helper(beta_1, beta_2, beta_3, x):
                return self.gamma * (self.gamma + self.kappa) * beta_2 * np.exp(-self.gamma * x) + self.kappa * beta_1
            
            def function_to_minimize_helper(beta_1, beta_2, beta_3, x):
                return np.exp(-self.gamma*x)*(-beta_1 * x + beta_2 * np.exp(-self.kappa * x)\
                                       -beta_3)
            def a_function_for_root(x):
                return function_for_root_helper(a_beta1, a_beta2, a_beta3, x)
            def a_function_for_root_derivative(x):
                return function_for_root_derivative_helper(a_beta1, a_beta2, a_beta3, x)
            
            def b_function_for_root(x):
                return function_for_root_helper(b_beta1, b_beta2, b_beta3, x)
            def b_function_for_root_derivative(x):
                return function_for_root_derivative_helper(b_beta1, b_beta2, b_beta3, x)
            
            def a_function_to_minimize(x):
                return function_to_minimize_helper(a_beta1, a_beta2, a_beta3, x)
            def b_function_to_minimize(x):
                return function_to_minimize_helper(b_beta1, b_beta2, b_beta3, x)
               
            if v[i-1] < 0 or v[i]<0 or v[i+1]<0:
                plot(v)
                show()
                print self.index_for_debug
            a_beta1_zero_guess = np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)) + np.log(v[i-1]) - np.log(v[i]))
            b_beta1_zero_guess = np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)) + np.log(v[i+1]) - np.log(v[i]))

            
            #When exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1)) overflow. Try compute manully or another way to compute
            #without using the lambertW function.
            if a_beta1==0:
                a_curr[i] = a_beta1_zero_guess
                b_curr[i] = b_beta1_zero_guess
            else:
                
                if q > 0:
                    try:
                        a_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(a_beta3, a_beta1) + np.true_divide(1, self.gamma)\
                            * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*a_beta1)\
                                           * self.gamma * a_beta2 *np.exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1))))
                        if math.isnan(a_curr[i]):
                            raise(FloatingPointError)
                    except FloatingPointError:
                        a_curr[i] = newton(func = a_function_for_root, x0=a_beta1_zero_guess, fprime = a_function_for_root_derivative)
                       
                    try:
                        if np.true_divide(self.gamma + self.kappa, self.kappa*b_beta1)  * self.gamma * b_beta2 *np.exp(-np.true_divide(self.gamma*(b_beta1 - self.kappa*b_beta3),self.kappa*b_beta1)) <= - np.exp(-1):
                            b_curr[i] = self.control_upper_bound
                        else:
                            b_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(b_beta3, b_beta1) + \
                            np.true_divide(1, self.gamma) * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*b_beta1)  \
                                            * self.gamma * b_beta2 *np.exp(-np.true_divide(self.gamma*(b_beta1 - self.kappa*b_beta3),self.kappa*b_beta1)), -1))
                    
                        if math.isnan(b_curr[i]):
                            raise(FloatingPointError)
                    except FloatingPointError:
                        #b_curr[i] = newton(func = b_function_for_root, x0=b_beta1_zero_guess, fprime = b_function_for_root_derivative)
                        opt_result = fmin_tnc(func=b_function_to_minimize, x0=np.array([b_beta1_zero_guess]), \
                                fprime=lambda x: np.exp(-self.kappa*x)*b_function_for_root(x), \
                                bounds=[(self.control_lower_bound, self.control_upper_bound)], disp=0)
                       

                        b_curr[i] = opt_result[0][0]
                else:
                    try:
                        b_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(b_beta3, b_beta1)\
                         + np.true_divide(1, self.gamma) * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*b_beta1) \
                             * self.gamma * b_beta2 *np.exp(-np.true_divide(self.gamma*(b_beta1 - self.kappa*b_beta3),self.kappa*b_beta1))))
                        if math.isnan(b_curr[i]):
                            raise(FloatingPointError)
                    except FloatingPointError:
                        b_curr[i] = newton(func = b_function_for_root, x0=b_beta1_zero_guess, fprime = b_function_for_root_derivative)
                    
                    try:
                        if np.true_divide(self.gamma + self.kappa, self.kappa*a_beta1)  * self.gamma * a_beta2 *np.exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1)) <= - np.exp(-1):
                            a_curr[i] = self.control_upper_bound
                        else:
                            a_curr[i] = np.true_divide(1, self.kappa) - np.true_divide(a_beta3, a_beta1)\
                                                + np.true_divide(1, self.gamma) * sp.real(lambertw( np.true_divide(self.gamma + self.kappa, self.kappa*a_beta1) \
                                                * self.gamma * a_beta2 *np.exp(-np.true_divide(self.gamma*(a_beta1 - self.kappa*a_beta3),self.kappa*a_beta1)), -1))
                        if math.isnan(a_curr[i]):
                            raise(FloatingPointError)
                    except FloatingPointError:
                        #a_curr[i] = newton(func = a_function_for_root, x0=a_beta1_zero_guess, fprime = a_function_for_root_derivative)
                        opt_result = fmin_tnc(func=a_function_to_minimize, x0=np.array([a_beta1_zero_guess]), \
                                fprime=lambda x: np.exp(-self.kappa*x)*a_function_for_root(x),\
                                 bounds=[(self.control_lower_bound, self.control_upper_bound)],disp=0)
                        a_curr[i] = opt_result[0][0]
                        
            if self.truncated_at_zero:
                a_curr[i] = 0 if a_curr[i] < 0 else a_curr[i]
                b_curr[i] = 0 if b_curr[i] < 0 else b_curr[i]
            a_curr[i] = self.control_upper_bound if self.F_1(a_curr[i], a_beta1, a_beta2, a_beta3)>0 else a_curr[i]
            b_curr[i] = self.control_upper_bound if self.F_1(b_curr[i], b_beta1, b_beta2, b_beta3)>0 else b_curr[i]
        return [a_curr[1:-1], b_curr[1:-1]]
    
    def linear_system_helper(self, v_curr, curr_control):
        a_curr, b_curr = curr_control
        co_left = self.A * self.delta_t * np.exp(-(self.kappa + self.gamma) * a_curr)
        co_right = self.A * self.delta_t * np.exp(-(self.kappa + self.gamma) * b_curr)
        co_mid = 1 + self.delta_t * ( - 0.5 * self.sigma_s**2 * self.gamma**2 * self.implement_q_space[1:-1]**2 + self.gamma * self.beta * self.A * self.implement_q_space[1:-1] \
                                      * (np.exp(-self.kappa* a_curr) * a_curr - np.exp(-self.kappa * b_curr) * b_curr)\
                                       + self.A * (np.exp(-self.kappa* a_curr) + np.exp(-self.kappa* b_curr)))
        
        
        
        eq_right = v_curr.copy()
        eq_right[0] = 0
        eq_right[-1] = 0
        return [-self.BC.coef_at_minus_one_helper(co_left), -self.BC.coef_at_plus_one_helper(co_right), self.BC.coef_at_curr_helper(co_mid), eq_right]
        
    def run(self, K=None, use_cache=False):
        old_settings = np.seterr(all='raise')
        super(Poisson_expUtil_implicit, self).run( K=K, use_cache=use_cache)
        np.seterr(**old_settings)
        
    
    
    def simulate_one_step_forward(self, index):
        super(Poisson_expUtil_implicit, self).simulate_one_step_forward(index)
        curr_control_a, curr_control_b = self.control_at_current_point(index, self.q[-1])
       
        a_intensity = self.delta_t * self.A * np.exp(-self.kappa* curr_control_a)
        b_intensity = self.delta_t * self.A * np.exp(-self.kappa* curr_control_b)
        a_prob_0 = np.exp(-a_intensity)
        b_prob_0 = np.exp(-b_intensity)
        #Here we only want our intensity small enough that with extremely low probability that Poisson event could happen more than twice in a small time interval.
        random_a = np.random.random()
        random_b = np.random.random()
        delta_N_a = 0 if random_a < a_prob_0 else 1
        delta_N_b = 0 if random_b < b_prob_0 else 1
        a_prob_1 = np.exp(-a_intensity) * a_intensity
        b_prob_1 = np.exp(-b_intensity) * b_intensity
        if random_a > a_prob_0 + a_prob_1:
            print "too large A_intensity!", index
        if random_b > b_prob_0 + b_prob_1:
            print "too large B_intensity!", index
        
        delta_x = (self.s[-1] + curr_control_a) * delta_N_a - (self.s[-1] - curr_control_b) * delta_N_b
        delta_q = delta_N_b - delta_N_a
        delta_s = self.sigma_s*np.sqrt(self.delta_t)*np.random.normal(0,1,1) + self.delta_t * self.beta*(self.A* np.exp(-self.kappa * curr_control_a) * curr_control_a\
                         - self.A* np.exp(-self.kappa * curr_control_b) * curr_control_b)
        self.x.append(self.x[-1] + delta_x)
        self.q.append(self.q[-1] + delta_q)
        self.s.append(self.s[-1] + delta_s)
        self.simulate_control_a.append(curr_control_a)
        self.simulate_control_b.append(curr_control_b)        
        #self.a_intensity_simulate.append(a_intensity)
        #self.b_intensity_simulate.append(b_intensity)
            
    def __init__(self, BC_class, linear_system_helper, truncated_at_zero = True,   *args, **kwargs):
        super(Poisson_expUtil_implicit, self).__init__(*args, **kwargs)
        self.half_implement_I = self.half_I + self.extend_space   
        self.BC = BC_class(linear_system_helper, self.implement_I)
        self.truncated_at_zero = truncated_at_zero
    