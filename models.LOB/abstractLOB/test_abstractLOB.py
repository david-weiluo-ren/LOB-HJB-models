'''
Created on Jan 27, 2015

@author: weiluo
'''
import unittest
#from abstractLOB import AbstractLOB
import abstractLOB
import numpy as np

class SimpleLOB_Implementation(abstractLOB.AbstractLOB):
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
        self._q_space, self._delta_q  = np.linspace(-self.N, self.N, self.I, retstep = True)
    def __init__(self, *args, **kwargs):
        super(SimpleLOB_Implementation, self).__init__(*args, **kwargs)
         
    def terminal_condition(self):
        return -1 * np.ones(len(self.implement_q_space))
  
    
   
    def one_step_back(self, v_curr):
        abstractLOB.AbstractLOB.one_step_back(self, v_curr)
        return v_curr

class BrownianMotionExpUtilImplicit(SimpleLOB_Implementation):
    pass
        
class Test(unittest.TestCase):


    def test_simple_implementation(self):
        myObj = SimpleLOB_Implementation()
        print myObj.delta_q
        assert(myObj.implement_I == myObj.I+4)
        print myObj.implement_q_space[:2], myObj.implement_q_space[-2:]
        print myObj.q_space[:10]
        print myObj.q_space[-10:]
        print myObj.v_init
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()