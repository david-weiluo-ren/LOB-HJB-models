'''
Created on Apr 21, 2015

@author: weiluo
'''
from abstractLOB.poisson_OU_explicit_marketOrder import Poisson_OU_explicit_marketOrder
from remoteExeFiles.SimulateSingle import mean_var_cmp, dumpData
if __name__ == '__main__':
    objNormalAlpha = Poisson_OU_explicit_marketOrder(A=20, s_0=2.0, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=4.0, delta_t=0.0001, N=10.0, half_I_S=60,\
                                  beta=0, q_0=0.0, alpha=10,gamma=1.0, s_long_term_mean=2.0, l=0.05,\
                                  half_I=100)
    objNormalAlpha.run()
    objNormalAlpha_lessRiskAversion_largeN = Poisson_OU_explicit_marketOrder(A=20, s_0=2.0, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=4.0, delta_t=0.0001, N=10.0, half_I_S=60,\
                                  beta=0, q_0=0.0, alpha=10,gamma=0.3, s_long_term_mean=2.0, l=0.05,\
                                  half_I=100)
    objNormalAlpha_lessRiskAversion_largeN.run()
    objNormalAlpha_greaterRiskAversion = Poisson_OU_explicit_marketOrder(A=20, s_0=2.0, kappa=1.5, sigma_s=3, num_time_step=10000,\
                                  half_S=4.0, delta_t=0.0001, N=10.0, half_I_S=60,\
                                  beta=0, q_0=0.0, alpha=10,gamma=3, s_long_term_mean=2.0, l=0.05,\
                                  half_I=100)
    objNormalAlpha_greaterRiskAversion.run()
    
    result = []
    simulate_num = 10000
    q_0 = 3
    
    result.append(mean_var_cmp(objNormalAlpha, simulate_num, q_0))
    result.append(mean_var_cmp(objNormalAlpha_lessRiskAversion_largeN, simulate_num, q_0))
    result.append(mean_var_cmp(objNormalAlpha_greaterRiskAversion, simulate_num, q_0))
    dumpData(result, "marketOrderDiffRiskAversion")

