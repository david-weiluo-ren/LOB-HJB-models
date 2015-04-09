'''
Created on Mar 4, 2015

@author: weiluo
'''
import abc
from abstractLOB import AbstractLOB
import numpy as np
class Abstract_OU_LOB(AbstractLOB):
    

    __metaclass__ = abc.ABCMeta
    def user_friendly_list_of_array(self, list_of_array, cache_index=None, front_extend_space=None, behind_extend_space=None, call_func=None):
        if front_extend_space is None:
            front_extend_space = self.extend_space
        if behind_extend_space is None:
            behind_extend_space = self.extend_space
        return [arr[front_extend_space:(-behind_extend_space), front_extend_space:(-behind_extend_space)] for arr in list_of_array]
    @property
    def extend_space(self):
        return int(self._extend_space)
    def compute_q_space(self):
        self._q_space  = np.linspace(-self.N, self.N, self.I)
        self._delta_q = 1
    
    @property
    def half_I(self):
        return int(self.N)
    
    @property
    def q_space(self):
        return self._q_space
    
    @property
    def delta_q(self):
        return self._delta_q

    def compute_s_space(self):
        self._s_space, self._delta_s  = np.linspace(self.s_long_term_mean-self.half_S, self.s_long_term_mean+self.half_S, self.I_S, retstep = True)
 
    @property
    def s_space(self):
        return self._s_space
    
    @property
    def delta_s(self):
        return self._delta_s
    def terminal_condition(self):
        pass
    @abc.abstractmethod
    def terminal_condition_real(self):
        pass
    def __init__(self, alpha = 5, half_S = 2.0, half_I_S=40, s_long_term_mean=None, *args, **kwargs):
        super(Abstract_OU_LOB, self).__init__(*args, **kwargs)

        self.alpha = alpha
        self._s_space = None
        self._delta_s = None
        self.half_S = half_S
        self.half_I_S = half_I_S
        self.I_S = 2*self.half_I_S + 1
        self.s_long_term_mean = s_long_term_mean if s_long_term_mean is not None else self.s_0

        self.compute_s_space()
        self.implement_s_space = np.hstack((-np.arange(self.extend_space, 0, -1) * self.delta_s + self.s_space[0],\
                                             self.s_space, \
                                             np.arange(1, self.extend_space+1) * self.delta_s + self.s_space[-1]))
        self.implement_S = len(self.implement_s_space)
        self.v_init = self.terminal_condition_real()
        self.s_drift_impact = [0]
        self.s_drift_OU = [0]
        self._a_price = []
        self._b_price = []
        self.simulate_price_a = []
        self.simulate_price_b = [] 
    def truncate_at_zero(self, arr):
        return np.maximum(arr, 0)    
    def q_to_index_for_simulate_control(self, q):
        if q > self.N or q < -self.N:
                print q, self.N
                self.failed_simulation += 1
                raise Exception("Too large inventory")
        curr_control_vector_length = np.shape(self._a_control[0])[1]
        return int(np.true_divide(q, self.delta_q)) + (curr_control_vector_length-1)/2
    
    def s_to_index_for_simulate_control(self, s):
        if s > self.s_long_term_mean+self.half_S or s < self.s_long_term_mean-self.half_S:
                print "overflow S =", s, self.S
                self.failed_simulation += 1
                raise Exception("Too large price")
        curr_control_vector_length = np.shape(self._a_control[0])[0]
        return int(np.true_divide(s-self.s_long_term_mean, self.delta_s)) + (curr_control_vector_length-1)/2
    
    def control_at_current_point(self, index, curr_q, curr_s):
        
        curr_control_a = self._a_control[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) ,self.q_to_index_for_simulate_control(curr_q)]
        curr_control_b = self._b_control[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) ,self.q_to_index_for_simulate_control(curr_q)]
        return [curr_control_a, curr_control_b] 
    def price_at_current_point(self, index, curr_q, curr_s):
        
        curr_price_a = self._a_price[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) ,self.q_to_index_for_simulate_control(curr_q)]
        curr_price_b = self._b_price[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) ,self.q_to_index_for_simulate_control(curr_q)]
        return [curr_price_a, curr_price_b]    
    
    def simulate_one_step_forward_helper(self, index, random_a, random_b, random_s):

        curr_control_a, curr_control_b = self.control_at_current_point(index, self.q[-1], self.s[-1])
        curr_price_a, curr_price_b = self.price_at_current_point(index, self.q[-1], self.s[-1])

        a_intensity = self.delta_t * self.A * np.exp(-self.kappa* curr_control_a)
        b_intensity = self.delta_t * self.A * np.exp(-self.kappa* curr_control_b)
        a_prob_0 = np.exp(-a_intensity)
        b_prob_0 = np.exp(-b_intensity)
        #Here we only want our intensity small enough that with extremely low probability that Poisson event could happen more than twice in a small time interval.
        
        delta_N_a = 0 if random_a <= a_prob_0 else 1
        delta_N_b = 0 if random_b <= b_prob_0 else 1
        a_prob_1 = np.exp(-a_intensity) * a_intensity
        b_prob_1 = np.exp(-b_intensity) * b_intensity
        if random_a > a_prob_0 + a_prob_1:
            print "too large A_intensity!", index
        if random_b > b_prob_0 + b_prob_1:
            print "too large B_intensity!", index
        
        delta_x = (self.s[-1] + curr_control_a) * delta_N_a - (self.s[-1] - curr_control_b) * delta_N_b
        delta_q = delta_N_b - delta_N_a
        delta_s_price_impact_part = self.delta_t * self.beta*(self.A* np.exp(-self.kappa * curr_control_a) \
                         - self.A* np.exp(-self.kappa * curr_control_b) )
        delta_s_OU_part = self.alpha*(self.s_long_term_mean-self.s[-1])*self.delta_t
        delta_s_drift_part = delta_s_price_impact_part + delta_s_OU_part
        delta_s = self.sigma_s*np.sqrt(self.delta_t) * random_s + delta_s_drift_part 
        self.x.append(self.x[-1] + delta_x)
        self.q.append(self.q[-1] + delta_q)
        self.s.append(self.s[-1] + delta_s)
        self.simulate_control_a.append(curr_control_a)
        self.simulate_control_b.append(curr_control_b) 
        self.simulate_price_a.append(curr_price_a)
        self.simulate_price_b.append(curr_price_b)         
        self.s_drift.append(self.s_drift[-1] + delta_s_drift_part) 
        self.s_drift_impact.append(self.s_drift_impact[-1] + delta_s_price_impact_part)
        #self.s_drift_OU.append(self.s_drift_OU[-1] +  delta_s_price_impact_part) 
        self.s_drift_OU.append(self.s_drift_OU[-1] + delta_s_OU_part)      
        #self.a_intensity_simulate.append(a_intensity)
        #self.b_intensity_simulate.append(b_intensity)    
        
    def simulate_one_step_forward(self, index):
        self.simulate_one_step_forward_helper(index, np.random.random(), np.random.random(), np.random.normal(0,1,1))
        
    def simulate_one_step_forward_use_givenRandom(self, index, randomSource):
        self.simulate_one_step_forward_helper(index, randomSource[0][index], randomSource[1][index], randomSource[2][index])
    def generate_random_source(self):
        return [np.random.random(self.num_time_step), np.random.random(self.num_time_step), np.random.normal(0, 1, self.num_time_step)]
    def init_forward_data(self, q_0 = None, x_0 = None, s_0 = None ):
        super(Abstract_OU_LOB, self).init_forward_data(q_0, x_0, s_0)
        self.simulate_price_a[:] = []
        self.simulate_price_b[:] = []
        
