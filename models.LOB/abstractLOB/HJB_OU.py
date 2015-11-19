import numpy as np
from scipy.sparse.linalg import spsolve
from scipy import sparse
from scipy.linalg import solve 
from collections import namedtuple
from scipy.stats import norm
from scipy import linalg
from numpy import indices

class HJB_OU_solver(object):
    '''
    Numerical solver of the HJB equation from an HFT model with 
    an OU reference price.
    One use case is
    
    import HJB_OU
    obj = HJB_OU.HJB_OU_solver()
    obj.run()
    for i in xrange(1, obj.implement_I - 1):
        plot( obj._a_price[-1][i * obj.implement_S + 1 : (i + 1) * obj.implement_S])
    '''

    """
    Parameters:
    (1) gamma, A, kappa, sigma_s, alpha, s_long_term_mean, lambda_tilde 
        are the ones in the model
    (2) half_I: Q in the model, namely half length of the q-space.
    (3) half_S: S in the model, namely half length of the s-space
    (4) half_I_S: Number of grid points between s_long_term_mean 
        and "s_long_term_mean + half_S" in the discretized s-space
    (5) delta_t, num_time_step: T = delta_t * num_time_step
    (6) gueant_boundary: True if using Gueant's boundary condition, 
        False if using zero second-derivative boundary condition.
    """
    def __init__(self, gamma = 1.0, A = 10, kappa = 1.5, 
                 sigma_s = 3.0, alpha = 5.0, s_long_term_mean=5.0, lambda_tilde = 0, 
                 half_I = 10, half_S = 3.0, half_I_S=300, delta_t = 0.001,
                 x_0 = 0, q_0 = 0, s_0 = None,
                 num_time_step = 10000, extend_space = 0, boundary_factor = 0, quadratic_boundary_factor = 0,
                 iter_max = 2000, new_weight = 0.1, abs_threshold_power = -4, rlt_threshold_power = -3, 
                 data_storing_jump_size = -1, OU_step = 0, delta_t_factor = 2, diagonal_factor = 0, sign_factor = 1,
                 verbose = False, use_sparse=True, gueant_boundary = False, raise_overflow=False, record_time_lower_bound=-1):
        
        """
        Parameters in the model
        """
        self.gamma = gamma
        self.A = A
        self.kappa = kappa
        self.sigma_s = sigma_s # sigma in the model
        self.alpha = alpha
        self.s_long_term_mean = s_long_term_mean # mu in the model 
        self.lambda_tilde = lambda_tilde #lambda in the model

        """
        choose boundary condition: 
        (1) Gueant's bounded inventory condition
        (2) The one I used: zero second derivative at the boundaries.
        """
        self.linear_system = self.linear_system_gueant if gueant_boundary else self.linear_system_zero_2nd_derivative

        self.quadratic_boundary_factor = quadratic_boundary_factor

        self.data_storing_jump_size = data_storing_jump_size
        """
        Number of time steps we will compute. Here T = \Delta t * num_time_step.
        """
        self.num_time_step = num_time_step 
        self.OU_step = OU_step
        self.delta_t_factor = delta_t_factor
        self.diagonal_factor = diagonal_factor
        self.sign_factor = sign_factor
        self.record_time_lower_bound = record_time_lower_bound
        """
        Compute the q-space.
        
        We will extend the q-space specified by the user. 
        For instance, if the q-space considered by the user is [-Q, Q], 
        then here the implement_q_space = [-Q - extend_space, Q + extend_space]
        
        The extended part is expected to be affected by the boundary condition.
        """
        self.extend_space = int(extend_space)
        self.half_I = int(half_I)
        self.delta_t = delta_t
        self.I = 2 * self.half_I + 1
        self.q_space  = np.linspace(-half_I, half_I, self.I)
        if self.extend_space > 0:
            self.implement_I = self.I + 2 * self.extend_space
            self.implement_q_space = np.hstack((-np.arange(self.extend_space, 0, -1)  + self.q_space[0],\
                                             self.q_space, \
                                             np.arange(1, self.extend_space+1)  + self.q_space[-1]))
        else:
            self.implement_I = self.I
            self.implement_q_space = self.q_space
        
        """
        Compute the s-space.
        We extend it similarly to the case of q-space.
        """        
        self.half_S = half_S
        self.half_I_S = half_I_S
        self.I_S = 2*self.half_I_S + 1
        self.s_space, self.delta_s  = np.linspace(self.s_long_term_mean-self.half_S, self.s_long_term_mean+self.half_S, self.I_S, retstep = True)
        
        if self.extend_space > 0:   
            self.implement_s_space = np.hstack((-np.arange(self.extend_space, 0, -1) * self.delta_s + self.s_space[0],
                                             self.s_space, 
                                             np.arange(1, self.extend_space+1) * self.delta_s + self.s_space[-1]))
        else:
            self.implement_s_space = self.s_space
            
        self.implement_S = len(self.implement_s_space)

        """
        Use the terminal condition of HJB equation to compute
        the v function at initial time (use \tau instead of t).
        """
        self.v_init = self.terminal_condition()
        self.v_init_PC = self.terminal_condition_PC()
        """
        Parameter used in iterations.
        """
        self.iter_max = iter_max
        self.new_weight = new_weight
        
        if abs_threshold_power > 0:
            raise Exception("the power of absolute threshold should be negative")
        if rlt_threshold_power > 0:
            raise Exception("the power of relative threshold should be negative")
        
        self.abs_threshold = 10**abs_threshold_power 
        self.rlt_threshold = 10**rlt_threshold_power
        self.use_sparse = use_sparse
        self.verbose = verbose
        self.boundary_factor = boundary_factor

        """
        The index of elements that are used when determining
        whether to stop the iteration. Basically the elements
        that are not in the extended space.
        """
        self.valid_index = self.construct_valid_index()

        """
        The data we store during the computation.
        (1) exp(-\delta^{a*})
        (2) exp(-\delta^{b*})
        (3) p^{a*}
        (4) p^{b*}
        (5) v(\tau, q, s)
        """
        self._a_exp_neg_control = []
        self._b_exp_neg_control = []
        self._a_price = []
        self._b_price = []
        self.value_function = []
        self.value_function_PC = []
        self.step_index = 0
        
        self.simulate_price_a = []
        self.simulate_price_b = []
        self.q_0 = q_0
        self.x_0 = x_0
        self.s_0 = self.s_long_term_mean if s_0 is None else s_0
        self.simulate_control_a = []
        self.simulate_control_b = []
        self.simulate_price_a_test = []
        self.simulate_price_b_test = []
        self.q = []
        self.q_a = [0]
        self.q_b = [0]
        self.x = []
        self.s = []
        self.s_drift = [0]
        self.failed_simulation = 0
        
        
        self.OU_transition_prob = self.construct_OU_transition_prob()
        self.exp_M_matrix =self.construct_exp_M_matrix()
        if raise_overflow:
            np.seterr(all='raise')
        else:
            np.seterr(all='warn')
        try:
            self.exp_kappa_SQ_matrix =  np.exp(-self.kappa * np.outer(self.implement_q_space, self.implement_s_space).reshape((1, -1))[0])
            self.exp_gamma_SQ_matrix =  np.exp(-self.gamma * np.outer(self.implement_q_space, self.implement_s_space).reshape((1, -1))[0])
        except:
            raise Exception("exp_SQ_matrix overflow / underflow")
        finally:
            np.seterr(all='warn')
            
    def construct_exp_M_matrix(self):
        np.seterr(all='raise')
        #C_array = np.ones(self.implement_I) * self.A * np.true_divide(self.gamma, self.kappa + self.gamma)\
        #    * (1 + np.true_divide(self.gamma, self.kappa)) ** (- np.true_divide(self.kappa, self.gamma)) * self.delta_t
        
        C_array = self.sign_factor * self.delta_t * np.ones(self.implement_I) * self.A\
            * (1 + np.true_divide(self.gamma, self.kappa)) ** (-1 - np.true_divide(self.kappa, self.gamma))
        diag_array = -self.delta_t * self.diagonal_factor * 0.5 * self.sigma_s ** 2 * self.gamma * self.kappa \
            * self.implement_q_space ** 2 
        m_matrix = sparse.spdiags([C_array, diag_array, C_array],
                      [-1, 0, 1], 
                      self.implement_I,
                      self.implement_I, 
                     format='csc')
        np.seterr(all='warn')
        
        return linalg.expm(m_matrix).todense()
    
    """
    Method used to run the solver. The only method
    that needs to be called after the object is constructed.
    After calling this method, the data are stored.
    """    
    def run(self, K = None):
        self._a_exp_neg_control[:] = []
        self._b_exp_neg_control[:] = []
        self._a_price[:] = []
        self._b_price[:] = []
        self.value_function = []
        self.step_index = 0

        if K is None:
            K = self.num_time_step
        
        v_curr = self.v_init 
           
        for i in xrange(K):
            if self.data_storing_jump_size < 0 or (
                self.data_storing_jump_size > 0 and i % self.data_storing_jump_size == 0):
                self.value_function.append(v_curr.copy())  
            
                    
            v_curr = self.one_step_back(v_curr, i)
            self.step_index += 1 
    
    def run_PC(self, value_function=None, start_index=0, K = None, exact=False, parallel=False, normalization=False):
        self._a_exp_neg_control[:] = []
        self._b_exp_neg_control[:] = []
        self._a_price[:] = []
        self._b_price[:] = []
        self.value_function_PC = []
        self.step_index = 0


        if K is None:
            K = self.num_time_step
        if value_function is None:
            v_curr = self.v_init_PC
        else:
            v_curr = value_function

           
        for i in xrange(start_index, K):
            if self.data_storing_jump_size < 0 or (
                self.data_storing_jump_size > 0 and i % self.data_storing_jump_size == 0):
                 
                if self.record_time_lower_bound <  0 or (
                   i > self.record_time_lower_bound):
                     
                    self.value_function_PC.append(v_curr.copy())  
            
            if exact:
                if parallel:
                    v_curr = self.one_step_back_PC_exact_parallel(v_curr, i)
                else:
                    v_curr = self.one_step_back_PC_exact(v_curr, i, normalization)
            else:      
                v_curr = self.one_step_back_PC(v_curr, i)
            self.step_index += 1 
    def run_PC_log(self, K = None, opt=False, iter=False):
        self._a_exp_neg_control[:] = []
        self._b_exp_neg_control[:] = []
        self._a_price[:] = []
        self._b_price[:] = []
        self.value_function[:] = []
        self.step_index = 0

        if K is None:
            K = self.num_time_step
        
        v_curr = self.v_init
           
        for i in xrange(K):
            if self.data_storing_jump_size < 0 or (
                self.data_storing_jump_size > 0 and i % self.data_storing_jump_size == 0):
                if self.record_time_lower_bound < 0 or i > self.record_time_lower_bound:
                    self.value_function.append(v_curr.copy())  
            
            if opt:       
                v_curr = self.one_step_back_PC_log_opt(v_curr, i)
            elif iter:
                v_curr = self.one_step_back_PC_log_iter(v_curr, i)
            else:
                v_curr = self.one_step_back_PC_log(v_curr, i)
            self.step_index += 1 
    
    def run__OU_PC_log_hybrid(self, value_function=None, start_index=0, OU_step=None, K=None, iter_flag=False, opt_flag=False):
        self._a_exp_neg_control[:] = []
        self._b_exp_neg_control[:] = []
        self._a_price[:] = []
        self._b_price[:] = []
        self.value_function[:] = []
        self.step_index = 0
        
        
        if K is None:
            K = self.num_time_step
        if OU_step is None:
            OU_step = self.OU_step
        if OU_step < 0:
            OU_step = K
        if value_function is None:
            v_curr = self.v_init
        else:
            v_curr = value_function
            
        for i in xrange(start_index, OU_step):
            if self.data_storing_jump_size < 0 or (
                self.data_storing_jump_size > 0 and i % self.data_storing_jump_size == 0):
                if self.record_time_lower_bound < 0 or i > self.record_time_lower_bound:
                    self.value_function.append(v_curr.copy())  
            
                    
            v_curr = self.one_step_back(v_curr, i)
            self.step_index += 1 
        
        for i in xrange(OU_step, K):
            if self.data_storing_jump_size < 0 or (
                self.data_storing_jump_size > 0 and i % self.data_storing_jump_size == 0):
                self.value_function.append(v_curr.copy())  
            
            if iter_flag:
                v_curr = self.one_step_back_PC_log_iter(v_curr, i)
            elif opt_flag:
                v_curr = self.one_step_back_PC_log_opt(v_curr, i)
            else:
                v_curr = self.one_step_back_PC_log(v_curr, i)
            self.step_index += 1 

    """
    Terminal Condition of HJB equation.
    """
    def terminal_condition(self):
        return np.outer(self.implement_q_space, self.implement_s_space).reshape((1, -1))[0] - \
            self.lambda_tilde * np.outer(self.implement_q_space**2, np.ones(self.implement_S)).reshape((1, -1))[0]
    def terminal_condition_PC(self):
        return np.exp( -self.gamma * (np.outer(self.implement_q_space, self.implement_s_space).reshape((1, -1))[0] - \
            self.lambda_tilde * np.outer(self.implement_q_space**2, np.ones(self.implement_S)).reshape((1, -1))[0]))
   
    

    def construct_OU_transition_prob(self):
        result = []
        s_space_modified_forward = self.implement_s_space.copy() + 0.5 * self.delta_s
        s_space_modified_forward = np.insert(s_space_modified_forward, 0, -np.infty)
        s_space_modified_forward[-1] = np.infty
        
        for s in self.implement_s_space:
            weight = np.exp(-self.alpha * 0.5 * self.delta_t * self.delta_t_factor) 
            loc = np.average([s, self.s_long_term_mean], weights=[weight, 1 - weight])
            scale = np.sqrt(np.true_divide(self.sigma_s ** 2, 2 * self.alpha) * (1 - np.exp(-2 * self.alpha * self.delta_t * self.delta_t_factor * 0.5)))
            cdfs = norm.cdf(s_space_modified_forward, loc=loc, scale=scale)
            result.append(cdfs[1:] - cdfs[:-1])
        return np.array(result)
    def construct_OU_transition_prob_tmp(self):
        weight = np.exp(-self.alpha * 0.5 * self.delta_t)
        OUVar = np.true_divide(self.sigma_s ** 2, 2 * self.alpha) * (1 - np.exp(-self.delta_t * self.alpha))
        OUMean = self.s_space * weight + self.s_long_term_mean * (1 - weight)
        
        s_space_modified_forward = self.s_space.copy() + 0.5 * self.delta_s
        s_space_modified_forward = np.insert(s_space_modified_forward, 0, -np.infty)
        s_space_modified_forward[-1] = np.infty
        
        result = norm.cdf(s_space_modified_forward.reshape((-1, 1)), loc=OUMean.reshape((1, -1)), scale=np.sqrt(OUVar))
        return  result[1:, :] - result[:-1, :]
    
    
    
    def one_step_back_PC(self, v_curr_PC, step_index):
        exp_neg_control = self.exp_neg_feedback_control_PC(v_curr_PC)
        optimal_price = self.optimal_price_PC(v_curr_PC)
        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            self._a_exp_neg_control.append(exp_neg_control[0])
            self._b_exp_neg_control.append(exp_neg_control[1])
            self._a_price.append(optimal_price[0])
            self._b_price.append(optimal_price[1])
        if step_index % 500 == 0:
            print step_index
        v_intermediate_PC = np.exp(-self.A * np.true_divide(self.gamma, self.gamma + self.kappa) * 
                                   (exp_neg_control[0] ** self.kappa + exp_neg_control[1] ** self.kappa) * 0.5 * self.delta_t * self.delta_t_factor) * v_curr_PC
        
        v_result_PC_matrix = np.dot(self.OU_transition_prob, v_intermediate_PC.reshape((self.implement_I, -1)).T)
        return v_result_PC_matrix.T.reshape((1, -1))[0]
    
    def one_step_back_PC_exact(self, v_curr_PC, step_index, normalization=False): #Following Tzu-Wei's new method
        np.seterr(all='raise')
 
        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            if self.record_time_lower_bound <  0 or (
                   step_index > self.record_time_lower_bound):
                exp_neg_control = self.exp_neg_feedback_control_PC(v_curr_PC)
                optimal_price = self.optimal_price_PC(v_curr_PC)
                self._a_exp_neg_control.append(exp_neg_control[0])
                self._b_exp_neg_control.append(exp_neg_control[1])
                self._a_price.append(optimal_price[0])
                self._b_price.append(optimal_price[1])
        if step_index % 500 == 0:
            print step_index
        if normalization:
            v_curr_PC = v_curr_PC / np.median(v_curr_PC)
        w_curr = (self.exp_kappa_SQ_matrix * v_curr_PC ** (-np.true_divide(self.kappa, self.gamma))).reshape((self.implement_I, -1))
        #print np.max(v_curr_PC), np.min(v_curr_PC)
        w_intermediate = np.asarray(np.dot(self.exp_M_matrix, w_curr)).reshape(1, -1)[0]
        
        v_intermediate_PC = self.exp_gamma_SQ_matrix * w_intermediate ** (-np.true_divide(self.gamma, self.kappa))
        
       
        v_result_PC_matrix = np.dot(self.OU_transition_prob, v_intermediate_PC.reshape((self.implement_I, -1)).T)
        
        np.seterr(all='warn')
        return v_result_PC_matrix.T.reshape((1, -1))[0]
    
    def one_step_back_PC_exact_parallel(self, v_curr_PC, step_index): #Following Tzu-Wei's new method
        np.seterr(all='raise')

        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            exp_neg_control = self.exp_neg_feedback_control_PC(v_curr_PC)
            optimal_price = self.optimal_price_PC(v_curr_PC)
            self._a_exp_neg_control.append(exp_neg_control[0])
            self._b_exp_neg_control.append(exp_neg_control[1])
            self._a_price.append(optimal_price[0])
            self._b_price.append(optimal_price[1])
        if step_index % 500 == 0:
            print step_index
        
        w_curr = (self.exp_kappa_SQ_matrix * v_curr_PC ** (-np.true_divide(self.kappa, self.gamma))).reshape((self.implement_I, -1))
        w_intermediate = np.asarray(np.dot(self.exp_M_matrix, w_curr)).reshape(1, -1)[0]
        
        v_intermediate_PC_q = self.exp_gamma_SQ_matrix * w_intermediate ** (-np.true_divide(self.gamma, self.kappa))
        
        v_intermediate_PC_s = np.dot(self.OU_transition_prob, v_curr_PC.reshape((self.implement_I, -1)).T).T.reshape((1, -1))[0]
        
        np.seterr(all='warn')
        return v_intermediate_PC_q + v_intermediate_PC_s - v_curr_PC

    '''
    def one_step_back_PC_log_exact(self, v_curr, step_index):
        np.seterr(all='raise')
        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            exp_neg_control = self.exp_neg_feedback_control(v_curr)
            optimal_price = self.optimal_price(v_curr)
            self._a_exp_neg_control.append(exp_neg_control[0])
            self._b_exp_neg_control.append(exp_neg_control[1])
            self._a_price.append(optimal_price[0])
            self._b_price.append(optimal_price[1])
        if step_index % 500 == 0:
            print step_index
        
        w_curr = np.exp(-self.kappa * (v_curr + np.outer(self.implement_q_space, self.implement_s_space).reshape((1, -1))[0]))
        w_intermediate = np.dot(self.exp_M_matrix, w_curr.reshape(self.implement_I, -1))
        exp_v_intermediate = np.exp(-self.kappa * np.outer(self.implement_q_space, self.implement_s_space).reshape((1, -1))[0]) \
            * (w_intermediate) ** (-np.true_divide(self.gamma, self.kappa))
        
        v_result_PC_matrix = np.dot(self.OU_transition_prob, exp_v_intermediate.reshape((self.implement_I, -1)).T)
        v_result_PC_matrix_flatten = v_result_PC_matrix.T.reshape((1, -1))[0]
        
        np.seterr(all='warn')

        return -np.true_divide(1, self.gamma) * np.log(v_result_PC_matrix_flatten) 
    '''
    
    def one_step_back_PC_log(self, v_curr, step_index, useMean=True):
        exp_neg_control = self.exp_neg_feedback_control(v_curr)
        optimal_price = self.optimal_price(v_curr)
        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            self._a_exp_neg_control.append(exp_neg_control[0])
            self._b_exp_neg_control.append(exp_neg_control[1])
            self._a_price.append(optimal_price[0])
            self._b_price.append(optimal_price[1])
        if step_index % 500 == 0:
            print step_index
        v_intermediate =  np.true_divide(self.A, self.gamma + self.kappa) * \
                                   (exp_neg_control[0] ** self.kappa + exp_neg_control[1] ** self.kappa) * 0.5 * self.delta_t * self.delta_t_factor + v_curr
        if useMean:
            K = np.mean(v_intermediate)
        else:
            K = np.min(v_intermediate)
        exp_v_intermediate_normalized = np.exp(-self.gamma * (v_intermediate - K))
        v_result_PC_matrix = np.dot(self.OU_transition_prob, exp_v_intermediate_normalized.reshape((self.implement_I, -1)).T)
        v_result_PC_matrix_flatten = v_result_PC_matrix.T.reshape((1, -1))[0]
        return -np.true_divide(1, self.gamma) * np.log(v_result_PC_matrix_flatten) + K
    
    def one_step_back_PC_log_iter(self, v_curr, step_index):
        if self.verbose:
            print step_index    
        exp_neg_control = self.exp_neg_feedback_control(v_curr)
        optimal_price = self.optimal_price(v_curr)
        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            self._a_exp_neg_control.append(exp_neg_control[0])
            self._b_exp_neg_control.append(exp_neg_control[1])
            self._a_price.append(optimal_price[0])
            self._b_price.append(optimal_price[1])

        
        v_tmp = v_curr.copy()
        itr = 0
        while(itr <= self.iter_max):
            exp_neg_control = self.exp_neg_feedback_control(v_tmp)

            v_intermediate =  np.true_divide(self.A, self.gamma + self.kappa) * \
                                   (exp_neg_control[0] ** self.kappa + exp_neg_control[1] ** self.kappa) * 0.5 * self.delta_t * self.delta_t_factor + v_curr
            K = np.min(v_intermediate)
            exp_v_intermediate_normalized = np.exp(-self.gamma * (v_intermediate - K))
            v_result_PC_matrix = np.dot(self.OU_transition_prob, exp_v_intermediate_normalized.reshape((self.implement_I, -1)).T)
            v_result_PC_matrix_flatten = v_result_PC_matrix.T.reshape((1, -1))[0]
            v_new = -np.true_divide(1, self.gamma) * np.log(v_result_PC_matrix_flatten) + K
            if self.close_enough(v_new, v_tmp):
                if step_index % 500 == 0:
                    print step_index, itr
                return v_new
            v_tmp = self.new_weight * v_new + (1 - self.new_weight) * v_tmp
            itr += 1
        
        
        print step_index
        raise Exception('iteration cannot converge!')  
        
    
    def one_step_back_PC_log_opt(self, v_curr, step_index):
        exp_neg_control = self.exp_neg_feedback_control(v_curr)
        optimal_price = self.optimal_price(v_curr)
        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            self._a_exp_neg_control.append(exp_neg_control[0])
            self._b_exp_neg_control.append(exp_neg_control[1])
            self._a_price.append(optimal_price[0])
            self._b_price.append(optimal_price[1])
        if step_index % 500 == 0:
            print step_index
        v_intermediate =  np.true_divide(self.A, self.gamma + self.kappa) * \
                                   (exp_neg_control[0] ** self.kappa + exp_neg_control[1] ** self.kappa) * 0.5 * self.delta_t + v_curr
        
        v_intermediate_reshaped = v_intermediate.reshape((self.implement_I, -1))
        result = []

        for prob_row in self.OU_transition_prob:

            indices = prob_row > 10 ** (-7)
            prob_row[~indices] = 0
            v_intermediate_reshaped_copy = v_intermediate_reshaped.copy()
            v_intermediate_reshaped_copy[:, ~indices] = 0
            #K_array = v_intermediate_reshaped_copy.sum(axis=1) / sum(indices)
            K_array = np.min(v_intermediate_reshaped_copy[:, indices])
            v_intermediate_reshaped_copy[:, indices] = v_intermediate_reshaped_copy[:, indices] - K_array.reshape((-1, 1))            
            exp_v_intermediate_normalized = np.exp(-self.gamma * v_intermediate_reshaped_copy)
            v_result_PC_matrix = np.dot(prob_row, exp_v_intermediate_normalized.T)
            result.append( -np.true_divide(1, self.gamma) * np.log(v_result_PC_matrix) + np.array(K_array))
            
        v_result_PC_matrix_flatten = np.array(result).T.reshape((1, -1))[0]
        return v_result_PC_matrix_flatten
    
    """
    Use iteration to solve for the v function at next temporal position.
    """
    def one_step_back(self, v_curr, step_index):
        
        exp_neg_control= self.exp_neg_feedback_control(v_curr)
        optimal_price = self.optimal_price(v_curr)
        
        if self.data_storing_jump_size < 0 or (
            self.data_storing_jump_size > 0 and step_index % self.data_storing_jump_size == 0):
            if self.record_time_lower_bound < 0 or step_index > self.record_time_lower_bound:
                self._a_exp_neg_control.append(exp_neg_control[0])
                self._b_exp_neg_control.append(exp_neg_control[1])
                self._a_price.append(optimal_price[0])
                self._b_price.append(optimal_price[1])

        v_tmp = v_curr 
        iter_count = 0
        while True:
            curr_exp_neg_control = self.exp_neg_feedback_control(v_tmp)
            v_new = self.one_iteration(v_curr, v_tmp, curr_exp_neg_control)
            if self.close_enough(v_new, v_tmp):
                if self.verbose:
                    print "the {}th iteration converges in {} iterations".format(self.step_index, iter_count),
                
                if step_index % 500 == 0:
                    print step_index, iter_count
                return v_new
            
            """
            v_tmp is the temporary value of v function used to computed the 
            expression on the RHS of discrete equation. It would be updated
            in every iteration until convergence.
            """
            v_tmp = self.new_weight * v_new + (1 - self.new_weight) * v_tmp
            
            iter_count += 1
            if iter_count > self.iter_max:
                print step_index
                raise Exception('iteration cannot converge!')  
           
    """
    One iteration for computing the v function. Basically solving a linear system.
    """
    def one_iteration(self, v_curr, v_iter_old, curr_exp_neg_control):
        eq_right, co_matrix = self.linear_system(v_curr, v_iter_old, curr_exp_neg_control) 

        if self.use_sparse:           
            return spsolve(co_matrix, eq_right)
        else:
            return solve(co_matrix.todense(), eq_right)  
        
        
    def exp_neg_feedback_control_PC(self, v_PC):
        
        return self.exp_neg_feedback_control(-np.true_divide(1, self.gamma) * np.log(v_PC))
    
    def optimal_price_PC(self, v_PC):
        return self.optimal_price(-np.true_divide(1, self.gamma) * np.log(v_PC))
    
    """
    Given v function at current time, compute 
    (exp^{-\delta^{a*}}, exp^{-\delta^{b*}})
    """
    def exp_neg_feedback_control(self, v):
        v_q_forward = np.ones(len(v))
        v_q_backward = np.ones(len(v))
        v_q_forward[ : -self.implement_S] = v[self.implement_S : ] - v[ : -self.implement_S] 
        v_q_backward[self.implement_S : ] = v[self.implement_S : ] - v[ : -self.implement_S]
       
        implement_s_space_casted = np.tile(self.implement_s_space, self.implement_I)
     
        
        exp_neg_optimal_a = (1+self.gamma/self.kappa)**(-1.0/self.gamma)\
                * np.exp(implement_s_space_casted - v_q_backward)
        exp_neg_optimal_a[: self.implement_S] = 0 if (self.boundary_factor == 0) else self.boundary_factor * \
            exp_neg_optimal_a[self.implement_S : (2 * self.implement_S)] / exp_neg_optimal_a[(2 * self.implement_S) : (3 * self.implement_S)] * \
            exp_neg_optimal_a[self.implement_S : (2 * self.implement_S)]
        exp_neg_optimal_b = (1+self.gamma/self.kappa)**(-1.0/self.gamma)\
                * np.exp(-implement_s_space_casted + v_q_forward)
        exp_neg_optimal_b[-self.implement_S:] = 0 if (self.boundary_factor == 0) else self.boundary_factor * \
            exp_neg_optimal_b[(-2 * self.implement_S) : (-self.implement_S)] / exp_neg_optimal_b[(-3 * self.implement_S) : (-2 * self.implement_S)] * \
            exp_neg_optimal_b[(-2 * self.implement_S) : (-self.implement_S)]
        return [exp_neg_optimal_a, exp_neg_optimal_b]
    
    """
    Given v function at current time, compute 
    (p^{a*}, p^{b*})
    """
    def optimal_price(self, v):
        v_q_forward = np.ones(len(v))
        v_q_backward = np.ones(len(v))
        v_q_forward[ : -self.implement_S] = v[self.implement_S : ] -  v[ : -self.implement_S] 
        v_q_backward[self.implement_S : ] = v[self.implement_S : ] -   v[ : -self.implement_S]
       
        LARGE_NUM = 100

        price_a = 1.0/self.gamma * np.log(1 + self.gamma/self.kappa) +  v_q_backward
        price_a[: self.implement_S] = LARGE_NUM
        price_b = - 1.0/self.gamma * np.log(1 + self.gamma/self.kappa) + v_q_forward
        price_b[-self.implement_S:] = LARGE_NUM

        return [price_a, price_b]            
                      
            
           
    """
    Compute the coefficient matrix on the LHS and the expression
    on the RHS of the discrete equation using the zero second-derivative
    boundary condition.
    """       
    def linear_system_zero_2nd_derivative(self, v_curr, v_iter_tmp, curr_exp_neg_control):
        a_curr_exp_neg, b_curr_exp_neg = curr_exp_neg_control
        totalLength = self.implement_I * self.implement_S
        
        """
        RHS of discrete equation.
        The positions corresponding to the boundary points 
        are set to be 0.
        """
        eq_right = v_curr.copy()
        eq_right[1:-1] += - 0.5 * self.sigma_s ** 2 * self.gamma * self.delta_t * ((v_iter_tmp[2:] - v_iter_tmp[:-2]) / (2 * self.delta_s)) ** 2\
            + self.A * self.delta_t / (self.kappa + self.gamma) * ((a_curr_exp_neg[1:-1]) ** self.kappa + (b_curr_exp_neg[1:-1]) ** self.kappa)
        eq_right[:self.implement_S] = self.quadratic_boundary_factor
        eq_right[-self.implement_S:] = self.quadratic_boundary_factor
        for i in xrange(1, self.implement_I-1):
            eq_right[i * self.implement_S] = 0
            eq_right[(i+1) * self.implement_S - 1] = 0

        
        """
        diagBlock_diagnal corresponds to the coefficients of v^{n+1}_{q, j},
        diagBlock_upper corresponds to the coefficients of v^{n+1}_{q, j + 1},
        diagBlock_lower corresponds to the coefficients of v^{n+1}_{q, j - 1}.
        If not considering the boundary condition, then the coeffficent matrix is 
        a tri-diagonal matrix, and those three arrays correspond to the non-zero
        diagonals in that tri-diagonal matrix.        
        """
        diagBlock_diagnal = np.ones(totalLength)
        diagBlock_upper = np.zeros(totalLength)
        diagBlock_lower = np.zeros(totalLength)
        
        """
        diagBlock_upper2 is the diagonal above diagBlock_upper,
        diagBlock_lower2 is the diagonal below diagBlock_lower.
        They are used to construct the boundary conditions for the
        boundaries in the s-space.
        """
        diagBlock_upper2 = np.zeros(totalLength)
        diagBlock_lower2 = np.zeros(totalLength)
        
        
        for i in xrange(1, self.implement_I-1):
            
            s_array = self.implement_s_space[1:-1] # Non-boundary points in the s-space
            
            
            """
            The arrays of indicator functions
            s_array_relSign: 1_{s_long_term_mean > s} - 1_{s_long_term_mean < s}
            s_array_relGreaterThanMean: 1_{s_long_term_mean < s}
            s_array_relLessThanMean: 1_{s_long_term_mean > s}
            """
            s_array_relSign = np.ones(len(s_array))
            s_array_relSign[s_array>self.s_long_term_mean] = -1 
            s_array_relLessThanMean = np.ones(len(s_array))
            s_array_relLessThanMean[s_array >= self.s_long_term_mean] = 0
            s_array_relGreaterThanMean = np.ones(len(s_array))
            s_array_relGreaterThanMean[s_array <= self.s_long_term_mean] = 0

            
            try:
                diagBlock_diagnal[(i * self.implement_S + 1) : ((i + 1) * self.implement_S - 1)] = \
                1 + (self.sigma_s / self.delta_s)**2 * self.delta_t\
                 + self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relSign                
                
            except Exception as e:
                print e
                            
                raise Exception()    
               
            diagBlock_upper[(i * self.implement_S + 1)] = -2
            diagBlock_upper[(i * self.implement_S + 2) : ((i + 1) * self.implement_S)] = \
                -0.5 * (self.sigma_s / self.delta_s)**2 * self.delta_t\
                -self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relLessThanMean            
            
            diagBlock_upper2[(i * self.implement_S + 2)] = 1
                  
            diagBlock_lower[((i + 1) * self.implement_S)-2] = -2
            diagBlock_lower[i * self.implement_S : ((i + 1) * self.implement_S - 2)] = \
                -0.5 * (self.sigma_s / self.delta_s)**2 * self.delta_t\
                +self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relGreaterThanMean

            diagBlock_lower2[((i + 1) * self.implement_S)-3] = 1
      
        
        
        """
        upperBlock_diagonal, upperBlock_diagonal2, lowerBlock_diagonal, lowerBlock_diagonal2
        are the ones used to construct boundary conditions for boundary points in q-space. 
        """
        upperBlock_diagonal = np.zeros(totalLength)
        upperBlock_diagonal[:2 * self.implement_S] = -2
        
        upperBlock_diagonal2 = np.zeros(totalLength)
        upperBlock_diagonal2[:3 * self.implement_S] = 1
        
        
        lowerBlock_diagonal = np.zeros(totalLength)
        lowerBlock_diagonal[(-2*self.implement_S):] = -2
        
        lowerBlock_diagonal2 = np.zeros(totalLength)
        lowerBlock_diagonal2[(-3*self.implement_S):] = 1
                    
        
        """
        Collect all the diagonals constructed above and create the coefficient matrix.
        """
        matrix_data = [lowerBlock_diagonal2, lowerBlock_diagonal, diagBlock_lower2, diagBlock_lower, diagBlock_diagnal,\
                        diagBlock_upper, diagBlock_upper2, upperBlock_diagonal, upperBlock_diagonal2]
        matrix_offset = [-2 * self.implement_S, -1*self.implement_S, -2, -1,0,1, 2, self.implement_S, 2 * self.implement_S]
    
        co_matrix = sparse.spdiags(matrix_data, matrix_offset, totalLength, totalLength, format='csc')

        return [eq_right, co_matrix]  
    
    
    """
    Compute the coefficient matrix on the LHS and the expression
    on the RHS of the discrete equation using Gueant's boundary condition.
    """ 
    def linear_system_gueant(self, v_curr, v_iter_tmp, curr_exp_neg_control):
        a_curr_exp_neg, b_curr_exp_neg = curr_exp_neg_control
        totalLength = self.implement_I * self.implement_S
        eq_right = v_curr.copy()
        eq_right[1:-1] += - 0.5 * self.sigma_s**2 * self.gamma * self.delta_t * ((v_iter_tmp[2:] - v_iter_tmp[:-2])/(2*self.delta_s))**2\
            + self.A*self.delta_t/(self.kappa + self.gamma) * ((a_curr_exp_neg[1:-1]) ** self.kappa + (b_curr_exp_neg[1:-1]) ** self.kappa)
        for i in xrange(self.implement_I):
            eq_right[i*self.implement_S] = 0
            eq_right[(i+1) * self.implement_S - 1] = 0

        
        
        diagBlock_diagnal = np.ones(totalLength)
        diagBlock_upper = np.zeros(totalLength)
        diagBlock_lower = np.zeros(totalLength)
        diagBlock_upper2 = np.zeros(totalLength)
        diagBlock_lower2 = np.zeros(totalLength)
        
        for i in xrange(self.implement_I):
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
                 + self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relSign

            except Exception as e:
                print e
                            
                raise Exception()    
               
            
            diagBlock_upper[(i * self.implement_S + 1)] = -2
            diagBlock_upper[(i * self.implement_S + 2) : ((i + 1) * self.implement_S)] = \
                -0.5 * (self.sigma_s / self.delta_s)**2 * self.delta_t\
                -self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relLessThanMean            
            
            diagBlock_upper2[(i * self.implement_S + 2)] = 1
                  
            diagBlock_lower[((i + 1) * self.implement_S)-2] = -2
            diagBlock_lower[i * self.implement_S : ((i + 1) * self.implement_S - 2)] = \
                -0.5 * (self.sigma_s / self.delta_s)**2 * self.delta_t\
                +self.alpha * (self.s_long_term_mean - s_array) * self.delta_t / self.delta_s * s_array_relGreaterThanMean

            diagBlock_lower2[((i + 1) * self.implement_S)-3] = 1
            
            
                    
        matrix_data = [diagBlock_lower2, diagBlock_lower, diagBlock_diagnal,\
                        diagBlock_upper, diagBlock_upper2]
        matrix_offset = [ -2, -1,0,1, 2]

    
        co_matrix = sparse.spdiags(matrix_data, matrix_offset, totalLength, totalLength, format='csc')
        return [eq_right, co_matrix]  
    
    def construct_valid_index(self):
        result = np.array([True] * (self.implement_I * self.implement_S))
        result[:self.implement_S] = False
        result[-self.implement_S:] = False
        for i in xrange(1, self.implement_I-1):
            result[i*self.implement_S] = False
            result[(i+1)*self.implement_S - 1] = False
        return result
    
    def close_enough(self, v_new, v_curr):
        return np.allclose(v_curr[self.valid_index], v_new[self.valid_index], self.rlt_threshold, self.abs_threshold)\
            and np.allclose(v_new[self.valid_index], v_curr[self.valid_index], self.rlt_threshold, self.abs_threshold)
    def simulate_one_step_forward(self, index):
        self.simulate_one_step_forward_helper(index, np.random.random(), np.random.random(), np.random.normal(0,1,1))
        
    def simulate_one_step_forward_use_givenRandom(self, index, randomSource):
        self.simulate_one_step_forward_helper(index, randomSource[0][index], randomSource[1][index], randomSource[2][index])
    def generate_random_source(self):
        return [np.random.random(self.num_time_step), np.random.random(self.num_time_step), np.random.normal(0, 1, self.num_time_step)]
   
    def simulate_forward(self, K = None, q_0 = None, x_0 = None, s_0 = None, randomSource = None):
        if self.data_storing_jump_size > 0:
            raise Exception("The value function data is not complete because data_storing_jump_size > 0.")
        self.init_forward_data(q_0, x_0, s_0)
        SimulationResult = namedtuple("SimulationResult",
                                      ["success",
                                       "control_a",
                                       "control_b",
                                       "q",
                                       "q_a",
                                       "q_b",
                                       "x",
                                       "s"])

        for index in xrange(K if (K is not None) else self.num_time_step):
            try:
                if randomSource is None:
                    self.simulate_one_step_forward(index)
                else:
                    self.simulate_one_step_forward_use_givenRandom(index, randomSource)
            except Exception, e:
                print e
                print "exit current simulation"
                return SimulationResult(False, self.simulate_control_a, self.simulate_control_b,\
                self.q, self.q_a, self.q_b, self.x, self.s)
            
        return SimulationResult(True, self.simulate_control_a, self.simulate_control_b,\
                self.q, self.q_a, self.q_b, self.x, self.s)
    
    def init_forward_data(self, q_0 = None, x_0 = None, s_0 = None ):
        self.simulate_price_a[:] = []
        self.simulate_price_b[:] = []
        q_0 = self.q_0 if q_0 is None else q_0
        x_0 = self.x_0 if x_0 is None else x_0
        s_0 = self.s_0 if s_0 is None else s_0
        self.simulate_control_a[:] = []
        self.simulate_control_b[:] = []
        self.simulate_price_a_test[:] = []
        self.simulate_price_b_test[:] = []
        self.q[:] = [q_0]
        self.q_a[:] = [0]
        self.q_b[:] = [0]
        self.x[:] = [x_0]
        self.s[:] = [s_0]
        self.s_drift[:]=[0]
        
    def control_at_current_point(self, index, curr_q, curr_s):
        
        curr_exp_neg_control_a = self._a_exp_neg_control[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        curr_exp_neg_control_b = self._b_exp_neg_control[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        return [curr_exp_neg_control_a, curr_exp_neg_control_b]

    def q_to_index_for_simulate_control(self, q):
        if q > self.half_I or q < -self.half_I:
                print q, self.half_I
                self.failed_simulation += 1
                raise Exception("Too large inventory")
       
        return int(q) + (self.implement_I - 1) / 2
    
    def s_to_index_for_simulate_control(self, s):
        if s > self.s_long_term_mean+self.half_S or s < self.s_long_term_mean-self.half_S:
                print "overflow S =", s, self.S
                self.failed_simulation += 1
                raise Exception("Too large price")
        return int(np.true_divide(s-self.s_long_term_mean, self.delta_s)) + (self.implement_S-1)/2


    def price_at_current_point(self, index, curr_q, curr_s):
        
        curr_price_a = self._a_price[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        curr_price_b = self._b_price[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) + self.q_to_index_for_simulate_control(curr_q) * self.implement_S]
        return [curr_price_a, curr_price_b]
    
    def simulate_one_step_forward_helper(self, index, random_a, random_b, random_s):
        curr_exp_neg_control_a, curr_exp_neg_control_b = self.control_at_current_point(index, self.q[-1], self.s[-1])
        curr_price_a, curr_price_b = self.price_at_current_point(index, self.q[-1], self.s[-1])
        
        a_spread = -np.log(curr_exp_neg_control_a) if curr_exp_neg_control_a > 0 else 100
        b_spread = -np.log(curr_exp_neg_control_b) if curr_exp_neg_control_b > 0 else 100
        a_intensity = self.delta_t * self.A * curr_exp_neg_control_a ** self.kappa
        b_intensity = self.delta_t * self.A * curr_exp_neg_control_b ** self.kappa
        a_prob_0 = np.exp(-a_intensity)
        b_prob_0 = np.exp(-b_intensity)
        #Here we only want our intensity small enough that with extremely low probability that Poisson event could happen more than twice in a small time interval.
        
        delta_N_a = 0 if random_a <= a_prob_0 else 1
        delta_N_b = 0 if random_b <= b_prob_0 else 1
        a_prob_1 = np.exp(-a_intensity) * a_intensity
        b_prob_1 = np.exp(-b_intensity) * b_intensity
        if random_a > a_prob_0 + a_prob_1:
            print "too large A_intensity!", index, random_a, a_prob_0 + a_prob_1
        if random_b > b_prob_0 + b_prob_1:
            print "too large B_intensity!", index, random_b, b_prob_0 + b_prob_1
        
        delta_x = (curr_price_a) * delta_N_a - (curr_price_b) * delta_N_b
        delta_q = delta_N_b - delta_N_a
        delta_s_OU_part = self.alpha * (self.s_long_term_mean - self.s[-1]) * self.delta_t
        delta_s = self.sigma_s * np.sqrt(self.delta_t) * random_s +delta_s_OU_part
        
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
        self.s_drift.append(self.s_drift[-1] + delta_s_OU_part) 
        
