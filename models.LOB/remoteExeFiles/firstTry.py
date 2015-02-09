'''
Created on Feb 6, 2015

@author: weiluo
'''
from abstractLOB.poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC
import numpy as np
import pickle

def saveObj1():
    params = {'gamma' : 1, 'sigma_s': 0.01, 'kappa': 0.3, 'A':0.05, 'num_time_step':10000, 'iter_max':500, 'N':7}
    
    myObj_largeT_beta_point1_smallSigma = Poisson_expUtil_implicit_NeumannBC(beta = 0.1, **params)
    myObj_largeT_beta_point1_smallSigma.run()
    myObj_largeT_beta_0_smallSigma = Poisson_expUtil_implicit_NeumannBC(beta = 0.1, **params)
    myObj_largeT_beta_0_smallSigma.run()
    myObj_largeT_beta_0_smallSigma_s_std = []
    myObj_largeT_beta_0_smallSigma_q_std = []
    myObj_largeT_beta_point1_smallSigma_s_std = []
    myObj_largeT_beta_point1_smallSigma_q_std = []
    for i in xrange(1000):
        tmp_result = myObj_largeT_beta_0_smallSigma.simulate_forward()
        if tmp_result[0]:
            myObj_largeT_beta_0_smallSigma_s_std.append(np.std(myObj_largeT_beta_0_smallSigma.s))
            myObj_largeT_beta_0_smallSigma_q_std.append(np.std(myObj_largeT_beta_0_smallSigma.q))
        
        tmp_result = myObj_largeT_beta_point1_smallSigma.simulate_forward()
        if tmp_result[0]:
            myObj_largeT_beta_point1_smallSigma_s_std.append(np.std(myObj_largeT_beta_point1_smallSigma.s))
            myObj_largeT_beta_point1_smallSigma_q_std.append(np.std(myObj_largeT_beta_point1_smallSigma.q))
    fileHandler = open("gmm_1_sgm01_kpp3_A05_N_7_beta1", 'r')
    
    pickle.dump([ myObj_largeT_beta_0_smallSigma_s_std, myObj_largeT_beta_0_smallSigma_q_std,\
    myObj_largeT_beta_point1_smallSigma_s_std, myObj_largeT_beta_point1_smallSigma_q_std], fileHandler)
if __name__ == '__main__':
    saveObj1()