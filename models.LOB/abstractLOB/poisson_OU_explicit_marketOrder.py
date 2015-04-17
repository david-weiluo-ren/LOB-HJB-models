'''
Created on Apr 16, 2015

@author: weiluo
'''
from abstract_OU_LOB import Abstract_OU_LOB
import numpy as np
from pylab import plot, show
class Poisson_OU_explicit_marketOrder(Abstract_OU_LOB):
    '''
    classdocs
    '''
    
    def compute_q_space(self):
        self._q_space, self._delta_q  = np.linspace(-self.N, self.N, self.I, retstep = True)
        
    @property
    def half_I(self):
        return self._half_I

    def __init__(self, l =0.01, LARGE_NUM=100, *args, **kwargs):
        super(Poisson_OU_explicit_marketOrder, self).__init__(*args, **kwargs)
        self.steps_for_one_share = self.half_I/self.N
        self.l = l
        self._xi = []
        self.simulate_xi = []
        self.simulate_accumulate_xi = [0]
        self.simulate_price_a_withLargeNum = []
        self.simulate_price_b_withLargeNum = []

        self.LARGE_NUM = LARGE_NUM
    def init_forward_data(self, q_0 = None, x_0 = None, s_0 = None ):
        super(Poisson_OU_explicit_marketOrder, self).init_forward_data(q_0, x_0, s_0)
        self.simulate_xi[:] = []
        self.simulate_accumulate_xi[:] = [0]
        self.simulate_price_a_withLargeNum[:] = []
        self.simulate_price_b_withLargeNum[:] = []
    def terminal_condition_real(self):
        return np.outer(self.implement_s_space, self.implement_q_space)
  
    def feedback_control(self, v, price=False):
        
        #If not price, return exp(-\delta^a*)*I_{\xi>=0} and exp(-delta^b*)*I_{xi<=0}
        
    
        
        #v_q_center = np.true_divide((v[1:-1, 2*self.steps_for_one_share:] - v[1:-1, :-2*self.steps_for_one_share]), 2 * self.delta_q)
        v_q_center = np.true_divide((v[1:-1, self.steps_for_one_share+1:-(self.steps_for_one_share-1)] - v[1:-1, (self.steps_for_one_share-1):-(self.steps_for_one_share+1)]), 2 * self.delta_q)
        v_q_center_reshape = np.reshape(v_q_center, (-1,))
        v_q_backward_one_share = v[1 : -1, : -2*self.steps_for_one_share] - v[ 1:-1, self.steps_for_one_share:-self.steps_for_one_share]
        v_q_forward_one_share = v[ 1:-1 , 2*self.steps_for_one_share : ] - v[1:-1, self.steps_for_one_share:-self.steps_for_one_share]
        v_q_backward_one_share_reshape = np.reshape(v_q_backward_one_share, (-1,))
        v_q_forward_one_share_reshape = np.reshape(v_q_forward_one_share, (-1,))

        
        s_space_casted = np.outer(self.implement_s_space, np.ones(len(self.implement_q_space) - 2 * self.steps_for_one_share))[1:-1,:]
        s_space_casted_reshape = np.reshape(s_space_casted, (-1,))
        sup_xi_func_positiveXi = np.zeros(len(s_space_casted_reshape))
        sup_xi_func_positiveXi[v_q_center_reshape > s_space_casted_reshape] = self.gamma * \
        np.true_divide((v_q_center_reshape[v_q_center_reshape > s_space_casted_reshape] - s_space_casted_reshape[v_q_center_reshape > s_space_casted_reshape]) ** 2, 4 * self.l)
        sup_xi_func_negativeXi = np.zeros(len(s_space_casted_reshape))
        sup_xi_func_negativeXi[v_q_center_reshape < s_space_casted_reshape] = self.gamma * \
        np.true_divide((v_q_center_reshape[v_q_center_reshape < s_space_casted_reshape] - s_space_casted_reshape[v_q_center_reshape < s_space_casted_reshape]) ** 2, 4 * self.l)
        F_a_optimal = (np.true_divide(self.kappa, self.kappa + self.gamma)) ** np.true_divide(self.kappa, self.gamma)\
         * np.exp(self.kappa * (s_space_casted_reshape + v_q_backward_one_share_reshape))
        F_b_optimal = (np.true_divide(self.kappa, self.kappa + self.gamma)) ** np.true_divide(self.kappa, self.gamma)\
         * np.exp(self.kappa * ( -s_space_casted_reshape + v_q_forward_one_share_reshape))
         
        optimal_xi = np.zeros(len(s_space_casted_reshape))
        bidLimitAskLimit_checkAskLimit_indice = np.all([v_q_center_reshape < s_space_casted_reshape, F_a_optimal > sup_xi_func_negativeXi], axis=0)
        bidLimitAskLimit_checkBidLimit_indice =  np.all([v_q_center_reshape > s_space_casted_reshape, F_b_optimal > sup_xi_func_positiveXi], axis=0)     
        bidLimitAskLimit_indice = np.any([bidLimitAskLimit_checkAskLimit_indice, bidLimitAskLimit_checkBidLimit_indice, (v_q_center_reshape == s_space_casted_reshape)], axis=0)   
        
        bidLimitAskMarket_indice = np.all([v_q_center_reshape < s_space_casted_reshape, F_a_optimal <= sup_xi_func_negativeXi], axis=0)
        bidMarketAskLimit_indice = np.all([v_q_center_reshape > s_space_casted_reshape, F_b_optimal <= sup_xi_func_positiveXi], axis=0)
        
        marketOrder_indice = np.any([bidLimitAskMarket_indice, bidMarketAskLimit_indice], axis=0)
        askLimitOrder_indice = np.any([bidLimitAskLimit_indice, bidMarketAskLimit_indice], axis=0)
        bidLimitOrder_indice = np.any([bidLimitAskLimit_indice, bidLimitAskMarket_indice], axis=0)
        
        optimal_xi[marketOrder_indice] = np.true_divide(v_q_center_reshape[marketOrder_indice] - s_space_casted_reshape[marketOrder_indice], 2 * self.l)
        F_xi_optimal = np.zeros(len(s_space_casted_reshape))
        F_xi_optimal[marketOrder_indice] = self.gamma * np.true_divide((v_q_center_reshape[marketOrder_indice] - s_space_casted_reshape[marketOrder_indice]) ** 2, 4 * self.l)
        F_a_optimal[~askLimitOrder_indice] = 0
        F_b_optimal[~bidLimitOrder_indice] = 0

        try:    
            optimal_a_reshape = np.zeros(len(s_space_casted_reshape))
            optimal_b_reshape = np.zeros(len(s_space_casted_reshape))
            if not price:
                optimal_a_reshape[askLimitOrder_indice] = np.true_divide(self.kappa, self.kappa + self.gamma)**(np.true_divide(1, self.gamma))\
                *np.exp(s_space_casted_reshape[askLimitOrder_indice] + v_q_backward_one_share_reshape[askLimitOrder_indice])
                optimal_b_reshape[bidLimitOrder_indice] = np.true_divide(self.kappa, self.kappa + self.gamma)**(np.true_divide(1, self.gamma))\
                *np.exp(-s_space_casted_reshape[bidLimitOrder_indice] + v_q_forward_one_share_reshape[bidLimitOrder_indice])
                return [np.reshape(optimal_a_reshape,np.shape(v_q_center)), np.reshape(optimal_b_reshape,np.shape(v_q_center)),\
                         np.reshape(optimal_xi, np.shape(v_q_center)), np.reshape(F_xi_optimal, np.shape(v_q_center)),\
                          np.reshape(F_a_optimal, np.shape(v_q_center)), np.reshape(F_b_optimal, np.shape(v_q_center))]
            else:
                optimal_a_reshape = np.ones(len(s_space_casted_reshape)) * self.LARGE_NUM
                optimal_b_reshape = -np.ones(len(s_space_casted_reshape)) * self.LARGE_NUM 
                optimal_a_reshape[askLimitOrder_indice] = np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)))\
                    -(v_q_backward_one_share_reshape[askLimitOrder_indice])
                optimal_b_reshape[bidLimitOrder_indice] = -1 * (np.true_divide(1, self.gamma)*(np.log(1+np.true_divide(self.gamma, self.kappa)))\
                    -(v_q_forward_one_share_reshape[bidLimitOrder_indice]))
                
                return [np.reshape(optimal_a_reshape,np.shape(v_q_center)), np.reshape(optimal_b_reshape,np.shape(v_q_center))]
        except Exception as e:
            print e
            plot(v)
            show()
            
            raise Exception()
        
            
    def one_step_back(self, v_curr, step_index=None):
        #print step_index
        v_new = np.ndarray((len(self.implement_s_space), len(self.implement_q_space)))
        
        exp_neg_optimal_a, exp_neg_optimal_b, optimal_xi, F_xi_optimal, F_a_optimal, F_b_optimal = self.feedback_control(v_curr)
        price_a, price_b= self.feedback_control(v_curr, price=True)
        

        
        v_ss = np.true_divide(v_curr[2:,:] + v_curr[:-2, :] - 2* v_curr[1:-1, :], self.delta_s**2)[:, self.steps_for_one_share:-self.steps_for_one_share]
        v_s = np.true_divide(v_curr[2:,:] - v_curr[:-2, :], 2*self.delta_s)[:, self.steps_for_one_share:-self.steps_for_one_share]
        
        s_space_casted = np.outer(self.implement_s_space, np.ones(len(self.implement_q_space) - 2 * self.steps_for_one_share))[1:-1, :]
        


        np.seterr(all="raise")
        
        try:
           
            v_new[1:-1, self.steps_for_one_share:-self.steps_for_one_share] = v_curr[1:-1, self.steps_for_one_share:-self.steps_for_one_share] +self.delta_t*((v_ss - self.gamma * v_s*v_s)*self.sigma_s**2*0.5\
                                + v_s * self.alpha*(self.s_long_term_mean-s_space_casted))
                               
        
            v_new[0,self.steps_for_one_share:-self.steps_for_one_share] = v_new[1, self.steps_for_one_share:-self.steps_for_one_share]
            v_new[-1, self.steps_for_one_share:-self.steps_for_one_share] = v_new[-2, self.steps_for_one_share:-self.steps_for_one_share]
            v_new[:, :self.steps_for_one_share] = np.outer(v_new[:, self.steps_for_one_share], np.ones(self.steps_for_one_share))
            v_new[:, -self.steps_for_one_share:] = np.outer(v_new[:, -self.steps_for_one_share - 1], np.ones(self.steps_for_one_share))
            
            
            self._a_control.append(exp_neg_optimal_a)
            self._b_control.append(exp_neg_optimal_b)
            self._a_price.append(price_a)
            self._b_price.append(price_b)
            self._xi.append(optimal_xi)
            return v_new
        except Exception as e:
            print e
            for i in xrange(np.shape(exp_neg_optimal_a)[1]):
                plot(exp_neg_optimal_a[:, i])
           
            show()
        
    def xi_at_current_point(self, index, curr_q, curr_s):
        
        curr_xi = self._xi[-1*(index+1)][self.s_to_index_for_simulate_control(curr_s) ,self.q_to_index_for_simulate_control(curr_q)]        
        return curr_xi 
    def simulate_one_step_forward_helper(self, index, random_a, random_b, random_s):

        curr_control_a, curr_control_b = self.control_at_current_point(index, self.q[-1], self.s[-1])
        curr_price_a, curr_price_b = self.price_at_current_point(index, self.q[-1], self.s[-1])
        curr_xi = self.xi_at_current_point(index, self.q[-1], self.s[-1])
        a_intensity = self.delta_t * self.A * curr_control_a**self.kappa if curr_xi >= 0 else 0
        b_intensity = self.delta_t * self.A * curr_control_b**self.kappa if curr_xi <= 0 else 0
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
        
        delta_x = (self.s[-1] + curr_control_a) * delta_N_a - (self.s[-1] - curr_control_b) * delta_N_b - curr_xi * (self.s[-1] + self.l * curr_xi)
        delta_q = delta_N_b - delta_N_a + curr_xi * self.delta_t
        delta_s_price_impact_part = 0
        #delta_s_price_impact_part = self.delta_t * self.beta*(self.A* np.exp(-self.kappa * curr_control_a) \
        #                 - self.A* np.exp(-self.kappa * curr_control_b) )
        delta_s_OU_part = self.alpha*(self.s_long_term_mean-self.s[-1])*self.delta_t
        delta_s_drift_part = delta_s_price_impact_part + delta_s_OU_part
        delta_s = self.sigma_s*np.sqrt(self.delta_t) * random_s + delta_s_drift_part 
        self.x.append(self.x[-1] + delta_x)
        self.q.append(self.q[-1] + delta_q)
        self.s.append(self.s[-1] + delta_s)
        self.simulate_control_a.append(curr_control_a)
        self.simulate_control_b.append(curr_control_b) 
        self.simulate_price_a_withLargeNum.append(curr_price_a)
        self.simulate_price_b_withLargeNum.append(curr_price_b)
        
        self.simulate_price_a.append(curr_price_a if curr_price_a != self.LARGE_NUM else self.s[-1])
        self.simulate_price_b.append(curr_price_b if curr_price_b != -self.LARGE_NUM else self.s[-1])

                 
        self.s_drift.append(self.s_drift[-1] + delta_s_drift_part) 
        self.s_drift_impact.append(self.s_drift_impact[-1] + delta_s_price_impact_part)
        #self.s_drift_OU.append(self.s_drift_OU[-1] +  delta_s_price_impact_part) 
        self.s_drift_OU.append(self.s_drift_OU[-1] + delta_s_OU_part)      
        #self.a_intensity_simulate.append(a_intensity)
        #self.b_intensity_simulate.append(b_intensity)
    
        self.simulate_xi.append(curr_xi)
        self.simulate_accumulate_xi.append(self.simulate_accumulate_xi[-1] + curr_xi*self.delta_t ) 
    
    
        