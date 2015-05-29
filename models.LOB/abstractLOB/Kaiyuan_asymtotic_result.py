'''
Created on May 12, 2015

@author: weiluo
'''
import numpy as np
from pylab import plot, show, legend, xlabel, ylabel
from poisson_OU_implicit import Poisson_OU_implicit
#tau = T - t
def zeta_0(obj, tau_index):
    tau = tau_index * obj.delta_t
    result = tau * np.true_divide(2 * obj.kappa * obj.A, (obj.kappa + obj.gamma) * obj.gamma) * (np.true_divide(obj.gamma, obj.kappa)\
                                                                                           -np.log(1 + np.true_divide(obj.gamma, obj.kappa)))
    result -= np.true_divide(obj.sigma_s ** 2 * obj.kappa * obj.A * obj.gamma, 2 * (obj.kappa + obj.gamma)) * \
            np.true_divide( np.exp(-2*obj.alpha * tau) - 1 + 2 * obj.alpha * tau, obj.alpha**2)
    return result
    
    
def zeta_1(obj, tau_index):
    tau = tau_index * obj.delta_t
    return obj.s_long_term_mean * (1 - np.exp(-obj.alpha * tau))

def zeta_2(obj, tau_index):
    tau = tau_index * obj.delta_t
    return np.true_divide(obj.sigma_s ** 2 * obj.gamma * (1 - np.exp(-2*obj.alpha * tau)), 2 * obj.alpha) 
def approx_value_function(obj, tau_index, q):
    return zeta_0(obj, tau_index) + q * zeta_1(obj, tau_index) - 0.5*q**2 * zeta_2(obj, tau_index)

def cmp_value_func_diff(obj, tau_index, subplot=None):
    zeta_q = []
    for i in xrange(obj.implement_I):
        zeta_q.append(obj._result[-tau_index]\
                  [int((i + 0.5) * obj.implement_S)])
    #plot(zeta_q, "o")

    #plot(approx(obj, 1, np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2)), "ro")
    if subplot is None:
        plot(np.array(zeta_q) - np.array(approx_value_function(obj, tau_index, np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2)) ), "o")
    else:
        subplot.plot(np.array(zeta_q) - np.array(approx_value_function(obj, tau_index, np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2)) ), "o")
  
    
def cmp_value_func(obj, tau_index, subplot=None):
    if tau_index < 0:
        tau_index = obj.num_time_step - tau_index
    zeta_q = []
    for i in xrange(obj.implement_I):
        zeta_q.append(obj._result[tau_index]\
                  [int((i + 0.5) * obj.implement_S)])
    if subplot is None:
        plot(np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2), zeta_q, "o", label="numerical result")

        plot(np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2), approx_value_function(obj, tau_index, np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2)), "ro",\
    label="Kaiyuan's approximation")
        legend(loc="best")
        xlabel("inventory")
        ylabel("value function v")
    else:
        subplot.plot(np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2), zeta_q, "o", label="numerical result")

        subplot.plot(np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2), approx_value_function(obj, tau_index, np.arange(-(obj.implement_I-1)/2, (obj.implement_I+1)/2)), "ro",\
    label="Kaiyuan's approximation")
        subplot.legend(loc="best")
        subplot.set_xlabel("inventory")
        subplot.set_ylabel("value function v")
    
def approx_a_price(obj, tau_index, q):
    constant = np.log(1+np.true_divide(obj.gamma, obj.kappa))/obj.gamma 
    
    return constant + zeta_1(obj, tau_index) - (q-0.5) * zeta_2(obj, tau_index)

def aPrice_dot(obj, k):
    result = []
    for i in xrange(0,obj.implement_I):
        
        result.append(obj._a_price[k][int((i+0.5) * obj.implement_S)])
    return result

def show_aPrice(obj, k , l = 10):
    for i in xrange(2,int(2 * obj.N )-1):
        plot(obj.implement_s_space[l+1:-l], obj._a_price[k][i * obj.implement_S+l+1:(i+1)*obj.implement_S-l])
    xlabel("ask price")
    ylabel("reference price")
myObj_largeAlpha_largeA_noBeta_midKappa_longTermMean3_normalTime = Poisson_OU_implicit(A=40, s_0=3, kappa=0.1, sigma_s=6, num_time_step=1000,\
                                  half_S=3, delta_t=0.001, N=8.0, half_I_S=240,\
                                  beta=-0, q_0=0.0, alpha=40,gamma=1.0, s_long_term_mean=3.0,\
                                   new_weight=0.02, iter_max=2000,\
                                  abs_threshold_power = -6, rlt_threshold_power = -5)
