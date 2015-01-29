'''
Created on Jan 28, 2015

@author: weiluo
'''
import unittest
from abstractDerivedLOB import AbstractImplicitLOB_NeumannBC
import numpy as np
class Test_AbstractImplicitLOB_NeumannBC(AbstractImplicitLOB_NeumannBC):
    
    def compute_q_space(self):
        self._q_space, self.delta_q  = np.linspace(-self.N, self.N, self.I, retstep = True)
    def __init__(self, *args, **kwargs):
        super(Test_AbstractImplicitLOB_NeumannBC, self).__init__(*args, **kwargs)
      
    
    def coef_at_curr(self, v_curr, curr_control):
        AbstractImplicitLOB_NeumannBC.coef_at_curr(self, v_curr, curr_control)
        return self.coef_at_curr_helper(np.zeros(self.implement_I-2))
    
    def coef_at_minus_one(self, v_curr, curr_control):
        AbstractImplicitLOB_NeumannBC.coef_at_minus_one(self, v_curr, curr_control)
        return -self.coef_at_curr_helper(np.ones(self.implement_I-2))
    
    def coef_at_plus_one(self, v_curr, curr_control):
        AbstractImplicitLOB_NeumannBC.coef_at_plus_one(self, v_curr, curr_control)
        return -self.coef_at_plus_one_helper(np.ones(self.implement_I-2))
    
    def equation_right(self, v_curr, curr_control):
        AbstractImplicitLOB_NeumannBC.equation_right(self, v_curr, curr_control)
        return self.equation_right_helper(np.zeros(self.implement_I-2))
    
    def feedback_control(self, v):
        AbstractImplicitLOB_NeumannBC.feedback_control(self, v)
    @property
    def half_I(self):
        return self.N
    
    @property
    def q_space(self):
        return self._q_space
    
    @property
    def delta_q(self):
        return 1
    def terminal_condition(self):
        return -1 * np.ones(len(self.implement_q_space))
  
    
class Test(unittest.TestCase):


    def test_construction(self):
        myObj = Test_AbstractImplicitLOB_NeumannBC(N=1, extend_space=1)
        print myObj.oneIter(None, None)

    def test_matrix(self):
        myObj = Test_AbstractImplicitLOB_NeumannBC(N=1, extend_space=0)
        myObj.coef_at_curr = lambda self,x: np.ones(myObj.implement_I)
        print myObj.oneIter(None, None)
    def test_matrix2(self):
        myObj = Test_AbstractImplicitLOB_NeumannBC(N=1, extend_space=1)
        myObj.coef_at_curr = lambda self,x: np.arange(1,myObj.implement_I+1)
        
        myObj.equation_right = lambda self,x: \
            np.hstack((0,np.arange(1,myObj.implement_I-1),0))
        print myObj.oneIter(None, None)
       
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()