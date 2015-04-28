'''
Created on Mar 4, 2015

@author: weiluo
'''
import numpy as np
from scipy.sparse.linalg import spsolve
from scipy import sparse
from abc import abstractmethod

from scipy.linalg import solve 
from pylab import plot, show
from numpy.linalg import matrix_rank
from abstract_OU_LOB import Abstract_OU_LOB
from abstractLOB import AbstractImplicitLOB

np.seterr(all='raise')



class Poisson_OU_implicit(Abstract_OU_LOB):
    
    def __init__(self, iter_max = 200, new_weight = 0.1, \
                 abs_threshold_power = -4, rlt_threshold_power = -3,\
                 use_sparse=True, *args,  **kwargs):
        super(Poisson_OU_implicit, self).__init__(*args, **kwargs)
        self.iter_max = iter_max
        self.new_weight = new_weight
        self.abs_threshold = 10**abs_threshold_power
        self.rlt_threshold = 10**rlt_threshold_power
        self.use_sparse = use_sparse
        
        
        self.simulate_price_a_test = []
        self.simulate_price_b_test = []
        self.valid_index = self.construct_valid_index()
    def construct_valid_index(self):
        result = np.array([True] * (self.implement_I * self.implement_S))
        result[:self.implement_S] = False
        result[-self.implement_S:] = False
        for i in xrange(1, self.implement_I-1):
            result[i*self.implement_S] = False
            result[(i+1)*self.implement_S - 1] = False
        return result

    '''
    Borrowed from AbstractImplicitLOB starts here
    
    I do think we should use a inner helper function to replace those code.
    '''
    def close_enough(self, v_new, v_curr):
        #return np.allclose(v_curr, v_new, self.rlt_threshold, self.abs_threshold)\
            #and np.allclose(v_new, v_curr, self.rlt_threshold, self.abs_threshold)
    
        return np.allclose(v_curr[self.valid_index], v_new[self.valid_index], self.rlt_threshold, self.abs_threshold)\
            and np.allclose(v_new[self.valid_index], v_curr[self.valid_index], self.rlt_threshold, self.abs_threshold)
    

    def one_iteration(self, v_curr, v_iter_old, curr_control, step_index):
        eq_right, co_matrix = self.linear_system(v_curr, v_iter_old, curr_control, step_index)
        if self.use_sparse:           
            return spsolve(co_matrix, eq_right)
        else:
            return solve(co_matrix.todense(), eq_right)
    
        
        
        
        
    def one_step_back(self, v_curr, step_index):
        
        return_control= self.exp_neg_feedback_control(v_curr)
        self._a_control.append(return_control[0])
        self._b_control.append(return_control[1])
        
        
        optimal_price = self.exp_neg_feedback_control(v_curr, True)
        self._a_price.append(optimal_price[0])
        self._b_price.append(optimal_price[1])

        
        v_tmp = v_curr
        iter_count = 0
        while True:
            curr_control = self.exp_neg_feedback_control(v_tmp)
            v_new = self.one_iteration(v_curr, v_tmp, curr_control, step_index)
            if self.close_enough(v_new, v_tmp):
                if self.verbose:
                    print "the {}th iteration converges in {} iterations".format(self.step_index, iter_count),
                
                if step_index % 500 == 0:
                    print step_index, iter_count
                return v_new
            v_tmp = self.new_weight * v_new + (1 - self.new_weight) * v_tmp
            
            iter_count += 1
            if iter_count > self.iter_max:
                print step_index
                for arr in self._a_control:
                    plot(arr)
                show()
                raise Exception('iteration cannot converge!')  

    '''
    Borrowed from AbstractImplicitLOB ends here
    '''

    def terminal_condition_real(self):
        return np.outer(self.implement_q_space, self.implement_s_space).reshape((1, -1))[0]
   
    def exp_neg_feedback_control(self, v , price=False):
        v_s_forward = np.zeros(len(v))
        v_s_backward =  np.zeros(len(v))
        v_q_forward = np.ones(len(v))
        v_q_backward = np.ones(len(v))
        
        
        v_s_forward[ : -1] = np.true_divide(v[1 : ] - v[ : -1], self.delta_s)
        v_s_backward[1 : ] = np.true_divide(v[1 : ] - v[ : -1], self.delta_s)

        v_q_forward[ : -self.implement_S] = v[self.implement_S : ] -  v[ : -self.implement_S] 
        v_q_backward[self.implement_S : ] = v[self.implement_S : ] -   v[ : -self.implement_S]
       
        implement_s_space_casted = np.tile(self.implement_s_space, self.implement_I)
        LARGE_NUM = 100
        a_critical_value = 1 + self.beta * self.gamma * v_s_forward
        b_critical_value = 1 - self.beta * self.gamma * v_s_backward
        a_critical_value[:self.implement_S] = 0
        a_critical_value[-self.implement_S:] = 0
        b_critical_value[:self.implement_S] = 0
        b_critical_value[-self.implement_S:] = 0
        for i in xrange(1, self.implement_I-1):
            a_critical_value[i*self.implement_S] = 0
            a_critical_value[(i+1) * self.implement_S - 1] = 0
            b_critical_value[i*self.implement_S] = 0
            b_critical_value[(i+1) * self.implement_S - 1] = 0
            

        
        
        if not price:
            exp_neg_optimal_a = np.zeros(len(v))
            exp_neg_optimal_b = np.zeros(len(v)) 
            
            exp_neg_optimal_a[a_critical_value>0] = (1+self.gamma/self.kappa)**(-1.0/self.gamma)\
                    * ((a_critical_value[a_critical_value>0])**(1.0/self.gamma))\
                    * np.exp(implement_s_space_casted[a_critical_value>0] - v_q_backward[a_critical_value>0])
            exp_neg_optimal_b[b_critical_value>0] = (1+self.gamma/self.kappa)**(-1.0/self.gamma)\
                    * ((b_critical_value[b_critical_value>0])**(1.0/self.gamma))\
                    * np.exp(-implement_s_space_casted[b_critical_value>0] + v_q_forward[b_critical_value>0])
            
            return [exp_neg_optimal_a, exp_neg_optimal_b]
        else :
            price_a = np.ones(len(v)) * LARGE_NUM
            price_b = -np.ones(len(v)) * LARGE_NUM
            price_a[a_critical_value > 0] = 1.0/self.gamma * np.log(1 + self.gamma/self.kappa)\
             - 1.0/self.gamma * np.log(a_critical_value[a_critical_value>0] )\
             +  v_q_backward[a_critical_value>0]
            price_b[b_critical_value > 0] = - 1.0/self.gamma * np.log(1 + self.gamma/self.kappa)\
            +  1.0/self.gamma * np.log(b_critical_value[b_critical_value>0] )\
            + v_q_forward[b_critical_value>0]
            return [price_a, price_b]
            
            
            
            
            
            
            
    def linear_system(self, v_curr, v_iter_old, curr_exp_neg_control, step_index):
        a_curr_exp_neg, b_curr_exp_neg = curr_exp_neg_control
        totalLength = self.implement_I * self.implement_S
        eq_right = v_curr.copy()
        eq_right[1:-1] += - 0.5 * self.sigma_s**2 * self.gamma * self.delta_t * ((v_iter_old[2:] - v_iter_old[:-2])/(2*self.delta_s))**2\
            + self.A*self.delta_t/(self.kappa + self.gamma) * ((a_curr_exp_neg[1:-1]) ** self.kappa + (b_curr_exp_neg[1:-1]) ** self.kappa)
        eq_right[:self.implement_S] = 0
        eq_right[-self.implement_S:] = 0
        for i in xrange(1, self.implement_I-1):
            eq_right[i*self.implement_S] = 0
            eq_right[(i+1) * self.implement_S - 1] = 0

        
        
        diagBlock_diagnal = np.ones(totalLength)
        diagBlock_upper = np.zeros(totalLength)
        diagBlock_lower = np.zeros(totalLength)
        for i in xrange(1, self.implement_I-1):
            q_value = self.implement_q_space[i]
            q_sign = 1 if q_value > 0 else -1
            q_negativeIndicator = 1 if q_value < 0 else 0
            q_positiveIndicator = 1 if q_value > 0 else 0
            s_array = self.implement_s_space[1:-1]
            s_array_relSign = np.ones(len(s_array))
            s_array_relSign[s_array>self.s_long_term_mean] = -1
            
            s_array_relLessThanMean = np.ones(len(s_array))
            s_array_relLessThanMean[s_array >= self.s_long_term_mean] = 0
            s_array_relGreaterThanMean = np.ones(len(s_array))
            s_array_relGreaterThanMean[s_array <= self.s_long_term_mean] = 0

            
            try:
                diagBlock_diagnal[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] = \
                1 + (self.sigma_s / self.delta_s)**2 * self.delta_t\
                 + self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relSign\
                 + self.A / (self.kappa + self.gamma) * self.delta_t / self.delta_s * self.beta * self.gamma\
                * (a_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa\
                    + b_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa)
                
                
                
            except Exception as e:
                print e
                            
                raise Exception()    
               
            diagBlock_upper[(i * self.implement_S + 1)] = -1
            diagBlock_upper[(i * self.implement_S + 2) : ((i + 1) * self.implement_S)] = \
            -0.5 * (self.sigma_s / self.delta_s)**2 * self.delta_t\
            -self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relLessThanMean\
            - self.A / (self.kappa + self.gamma) * self.delta_t / self.delta_s * self.beta * self.gamma\
             * a_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa
            
                  
            diagBlock_lower[((i + 1) * self.implement_S)-2] = -1
            diagBlock_lower[i * self.implement_S : ((i + 1) * self.implement_S - 2)] = \
            -0.5 * (self.sigma_s / self.delta_s)**2 * self.delta_t\
            +self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relGreaterThanMean\
             - self.A / (self.kappa + self.gamma) * self.delta_t / self.delta_s * self.beta * self.gamma\
             *b_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa

      
        
        upperBlock_diagonal = np.zeros(totalLength)
        upperBlock_diagonal[:2 * self.implement_S] = -1
        
        
        lowerBlock_diagonal = np.zeros(totalLength)
        lowerBlock_diagonal[(-2*self.implement_S):] = -1
        
    
    
                    
        matrix_data = [ lowerBlock_diagonal, diagBlock_lower, diagBlock_diagnal, diagBlock_upper, upperBlock_diagonal]
        matrix_offset = [-1*self.implement_S, -1,0,1, self.implement_S]
    
        co_matrix = sparse.spdiags(matrix_data, matrix_offset, totalLength, totalLength)
        tmp = co_matrix.todense()
        #print step_index,
        #print tmp.shape, matrix_rank(tmp)
        #print tmp
        #print "right: ", eq_right
        return [eq_right, co_matrix]



    def q_to_index_for_simulate_control(self, q):
        if q > self.N or q < -self.N:
                print q, self.N
                self.failed_simulation += 1
                raise Exception("Too large inventory")
       
        return int(np.true_divide(q, self.delta_q)) + (self.implement_I-1)/2
    
    def s_to_index_for_simulate_control(self, s):
        if s > self.s_long_term_mean+self.half_S or s < self.s_long_term_mean-self.half_S:
                print "overflow S =", s, self.S
                self.failed_simulation += 1
                raise Exception("Too large price")
        return int(np.true_divide(s-self.s_long_term_mean, self.delta_s)) + (self.implement_S-1)/2




    def control_at_current_point(self, index, curr_q, curr_s):
        
        curr_control_a = self._a_control[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        curr_control_b = self._b_control[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        return [curr_control_a, curr_control_b]


    def price_at_current_point(self, index, curr_q, curr_s):
        
        curr_price_a = self._a_price[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        curr_price_b = self._b_price[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        return [curr_price_a, curr_price_b]


    def simulate_one_step_forward_helper(self, index, random_a, random_b, random_s):
        curr_control_a, curr_control_b = self.control_at_current_point(index, self.q[-1], self.s[-1])
        curr_price_a, curr_price_b = self.price_at_current_point(index, self.q[-1], self.s[-1])
        LARGW_NUM=100
        
        a_spread = -np.log(curr_control_a) if curr_control_a > 0 else 100
        b_spread = -np.log(curr_control_b) if curr_control_b > 0 else 100
        a_intensity = self.delta_t * self.A * curr_control_a ** self.kappa
        b_intensity = self.delta_t * self.A * curr_control_b ** self.kappa
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
        
        
        
        
        delta_x = (self.s[-1] + a_spread) * delta_N_a - (self.s[-1] - b_spread) * delta_N_b
        delta_q = delta_N_b - delta_N_a
        delta_s_price_impact_part = self.delta_t * self.beta*(self.A* curr_control_a ** self.kappa \
                         - self.A* curr_control_b ** self.kappa )

        delta_s_OU_part = self.alpha*(self.s_long_term_mean-self.s[-1])*self.delta_t
        
        delta_s_drift_part = delta_s_price_impact_part + delta_s_OU_part
        delta_s = self.sigma_s*np.sqrt(self.delta_t) * random_s +delta_s_drift_part 
        self.x.append(self.x[-1] + delta_x)
        self.q.append(self.q[-1] + delta_q)
        self.s.append(self.s[-1] + delta_s)
        self.simulate_control_a.append(a_spread)
        self.simulate_control_b.append(b_spread) 
        self.simulate_price_a.append(curr_price_a)
        self.simulate_price_b.append(curr_price_b)
        self.simulate_price_a_test.append(self.s[-1] + a_spread)
        self.simulate_price_b_test.append(self.s[-1] - b_spread)


        #self.s_drift.append(self.s_drift[-1] +  delta_s_price_impact_part) 
        self.s_drift.append(self.s_drift[-1] + delta_s_drift_part) 
        self.s_drift_impact.append(self.s_drift_impact[-1] + delta_s_price_impact_part)
        
        self.s_drift_OU.append(self.s_drift_OU[-1] + delta_s_OU_part)    

    def init_forward_data(self, q_0 = None, x_0 = None, s_0 = None ):
        super(Poisson_OU_implicit, self).init_forward_data(q_0, x_0, s_0)
        self.simulate_price_a = []
        self.simulate_price_b = []
        self.simulate_price_a_test = []
        self.simulate_price_b_test = []

class Poisson_OU_implicit_truncateControlAtZero(Poisson_OU_implicit):
    
    def exp_neg_feedback_control(self, v , price=False):
        if price:
            price_a, price_b = super(Poisson_OU_implicit_truncateControlAtZero, self).exp_neg_feedback_control(v, price)
            implement_s_space_casted = np.tile(self.implement_s_space, self.implement_I)
            price_a = np.maximum(price_a, implement_s_space_casted)
            price_b = np.minimum(price_b, implement_s_space_casted)
            return [price_a, price_b]
        else:
            exp_neg_optimal = super(Poisson_OU_implicit_truncateControlAtZero, self).exp_neg_feedback_control(v, price)
            return np.minimum(1, exp_neg_optimal)


























class Poisson_OU_implicit_explode(Abstract_OU_LOB):
    '''
    classdocs
    '''

    def __init__(self, iter_max = 200, new_weight = 0.1, \
                 abs_threshold_power = -4, rlt_threshold_power = -2,\
                 use_sparse=True, *args,  **kwargs):
        super(Poisson_OU_implicit_explode, self).__init__(*args, **kwargs)
        self.iter_max = iter_max
        self.new_weight = new_weight
        self.abs_threshold = 10**abs_threshold_power
        self.rlt_threshold = 10**rlt_threshold_power
        self.use_sparse = use_sparse

    '''
    Borrowed from AbstractImplicitLOB starts here
    
    I do think we should use a inner helper function to replace those code.
    '''
    def close_enough(self, v_new, v_curr):
        return np.allclose(v_curr, v_new, self.rlt_threshold, self.abs_threshold)\
            and np.allclose(v_new, v_curr, self.rlt_threshold, self.abs_threshold)
    


    def one_iteration(self, v_curr, curr_control, step_index):
        eq_right, co_matrix = self.linear_system(v_curr, curr_control, step_index)
        if self.use_sparse:           
            return spsolve(co_matrix, eq_right)
        else:
            return solve(co_matrix.todense(), eq_right)
   
    def one_step_back(self, v_curr, step_index):
        
        
        v_tmp = v_curr
        iter_count = 0
        while True:
            curr_control = self.exp_neg_feedback_control(v_tmp)
            v_new = self.one_iteration(v_curr, curr_control, step_index)
            if self.close_enough(v_new, v_tmp):
                if self.verbose:
                    print "the {}th iteration converges in {} iterations".format(self.step_index, iter_count),
                return_control= self.exp_neg_feedback_control(v_new)
                self._a_control.append(return_control[0])
                self._b_control.append(return_control[1])
                print step_index, iter_count
                return v_new
            v_tmp = self.new_weight * v_new + (1 - self.new_weight) * v_tmp
            
            iter_count += 1
            if iter_count > self.iter_max:
                print step_index
                for arr in self.a_control:
                    plot(arr)
                show()
                raise Exception('iteration cannot converge!')  

    '''
    Borrowed from AbstractImplicitLOB ends here
    '''

    def twoDimCoordToOneDim(self, i_q, j_s):
        return i_q * self.implement_S + j_s

    def oneDimCoordToQ_index(self, loc):
        return loc/self.implement_S
    
    
    def oneDimCoordToQ_value(self, loc):
        return self.implement_q_space[self.oneDimCoordToQ_index(loc)]
    
    def linear_system(self, v_curr, curr_exp_neg_control, step_index):
        a_curr_exp_neg, b_curr_exp_neg = curr_exp_neg_control
        totalLength = self.implement_I * self.implement_S
        eq_right = v_curr.copy() 
        eq_right[:self.implement_S] = 0
        eq_right[-self.implement_S:] = 0
        for i in xrange(1, self.implement_I-1):
            eq_right[i*self.implement_S] = 0
            eq_right[(i+1) * self.implement_S - 1] = 0

        
        
        diagBlock_diagnal = np.ones(totalLength)
        diagBlock_upper = np.zeros(totalLength)
        diagBlock_lower = np.zeros(totalLength)
        for i in xrange(1, self.implement_I-1):
            q_value = self.implement_q_space[i]
            q_sign = 1 if q_value > 0 else -1
            q_negativeIndicator = 1 if q_value < 0 else 0
            q_positiveIndicator = 1 if q_value > 0 else 0
            s_array = self.implement_s_space[1:-1]
            s_array_relSign = np.ones(len(s_array))
            s_array_relSign[s_array>self.s_long_term_mean] = -1
            
            s_array_relLessThanMean = np.ones(len(s_array))
            s_array_relLessThanMean[s_array >= self.s_long_term_mean] = 0
            s_array_relGreaterThanMean = np.ones(len(s_array))
            s_array_relGreaterThanMean[s_array <= self.s_long_term_mean] = 0

            
            try:
                diagBlock_diagnal[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] = \
                1 - 0.5 * (self.sigma_s * self.gamma * q_value) ** 2 * self.delta_t\
                + 0.5 * self.sigma_s ** 2 * self.delta_t * (2 * self.gamma * q_value / self.delta_s * q_sign + 2/(self.delta_s ** 2))\
                + self.alpha * (self.s_long_term_mean - s_array) * self.delta_t \
                * (self.gamma * q_value + 1 / self.delta_s * s_array_relSign)\
                +self.beta * self.gamma * q_value *self.delta_t * self.A * (a_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa\
                                            - b_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa)\
                + (1+self.beta/self.delta_s)*self.A * self.delta_t * (a_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa \
                                            +  b_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa )
            except Exception as e:
                print e
                            
                raise Exception()    
               
            diagBlock_upper[(i * self.implement_S + 1)] = -1
            diagBlock_upper[(i * self.implement_S + 2) : ((i + 1) * self.implement_S)] = \
            self.delta_t / self.delta_s * 2 * self.gamma * q_value * (self.sigma_s)**2 *0.5 * q_negativeIndicator\
            - (self.sigma_s)**2 *0.5 * self.delta_t/(self.delta_s ** 2)\
             - self.delta_t/self.delta_s * self.alpha * (self.s_long_term_mean - s_array) * s_array_relLessThanMean\
             -self.beta * self.A *  a_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa * self.delta_t/self.delta_s
             
            diagBlock_lower[((i + 1) * self.implement_S)-2] = -1
            diagBlock_lower[i * self.implement_S : ((i + 1) * self.implement_S - 2)] = \
             -self.delta_t / self.delta_s * 2 * self.gamma * q_value * (self.sigma_s)**2 *0.5 * q_positiveIndicator\
             - (self.sigma_s)**2 *0.5 * self.delta_t/(self.delta_s ** 2)\
            + self.delta_t/self.delta_s * self.alpha * (self.s_long_term_mean - s_array) * s_array_relGreaterThanMean\
            -self.beta * self.A *  b_curr_exp_neg[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] ** self.kappa * self.delta_t/self.delta_s


        
        upperBlock_diagonal = -1 * np.ones(totalLength)
        upperBlock_diagonal[2 * self.implement_S:] = -self.A * self.delta_t * (b_curr_exp_neg[self.implement_S:(-self.implement_S)]) ** (self.kappa+self.gamma)
        for i in xrange(2, self.implement_I):
            upperBlock_diagonal[i * self.implement_S] = 0
            upperBlock_diagonal[(i+1) * self.implement_S - 1] = 0
            
        
        lowerBlock_diagonal = -1 * np.ones(totalLength)
        lowerBlock_diagonal[:(-2*self.implement_S)] = -self.A * self.delta_t * (a_curr_exp_neg[self.implement_S:(-self.implement_S)]) ** (self.kappa + self.gamma)
        for i in xrange(0, self.implement_I-2):
            lowerBlock_diagonal[i * self.implement_S] = 0
            lowerBlock_diagonal[(i+1) * self.implement_S - 1] = 0
    
    
    
                    
        matrix_data = [ lowerBlock_diagonal, diagBlock_lower, diagBlock_diagnal, diagBlock_upper, upperBlock_diagonal]
        matrix_offset = [-1*self.implement_S, -1,0,1, self.implement_S]
    
        co_matrix = sparse.spdiags(matrix_data, matrix_offset, totalLength, totalLength)
        tmp = co_matrix.todense()
        #print step_index,
        #print tmp.shape, matrix_rank(tmp)
        #print tmp
        #print "right: ", eq_right
        return [eq_right, co_matrix]

    
    
    
    
    def exp_neg_feedback_control(self, v):
        v_s_forward = np.zeros(len(v))
        v_s_backward =  np.zeros(len(v))
        v_q_forward = np.ones(len(v))
        v_q_backward = np.ones(len(v))
        
        
        v_s_forward[ : -1] = np.true_divide(v[1 : ] - v[ : -1], self.delta_s)
        v_s_backward[1 : ] = np.true_divide(v[1 : ] - v[ : -1], self.delta_s)

        v_q_forward[ : -self.implement_S] = np.true_divide(v[self.implement_S : ] ,   v[ : -self.implement_S]) 
        v_q_backward[self.implement_S : ] = np.true_divide(v[self.implement_S : ] ,   v[ : -self.implement_S])
       
        implement_q_space_casted = np.repeat(self.implement_q_space, self.implement_S)
        LARGE_NUM = 100
        
        
        exp_neg_optimal_a = np.zeros(len(v))
        exp_neg_optimal_b = np.zeros(len(v)) 
        a_critical_value = 1 + self.beta * self.gamma * implement_q_space_casted\
         - self.beta * np.true_divide(v_s_forward, v)
        b_critical_value = 1 - self.beta * self.gamma * implement_q_space_casted\
         - self.beta * np.true_divide(v_s_backward, v)
        
        exp_neg_optimal_a[a_critical_value>0] = (1+self.gamma/self.kappa)**(-1.0/self.gamma)\
                * ((a_critical_value[a_critical_value>0])**(1.0/self.gamma))\
                * (np.true_divide(1, v_q_backward[a_critical_value>0])**(-1.0/self.gamma))
        exp_neg_optimal_b[b_critical_value>0] = (1+self.gamma/self.kappa)**(-1.0/self.gamma)\
                * ((b_critical_value[b_critical_value>0])**(1.0/self.gamma))\
                *(v_q_forward[b_critical_value>0] ** (-1.0/self.gamma))
        
        return [exp_neg_optimal_a, exp_neg_optimal_b]
    
    
    
    def terminal_condition_real(self):
        return -1*np.ones(self.implement_I * self.implement_S)
   
    
    
    