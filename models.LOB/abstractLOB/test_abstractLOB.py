'''
Created on Jan 27, 2015

@author: weiluo
'''
import unittest
from abstractLOB import AbstractImplicitLOB, AbstractLOB

import numpy as np

class SimpleLOB_Implementation(AbstractLOB):
    @property
    def half_I(self):
        return self._half_I
    
    @property
    def q_space(self):
        return self._q_space
    
    @property
    def delta_q(self):
        return self._delta_q

    def __init__(self, *args, **kwargs):
        super(BrownianMotionExpUtilImplicit, self).__init__(*args, **kwargs)
        self._q_space, self._delta_q  = np.linspace(-self.N, self.N, self.I, retstep = True)
        
    def terminal_condition(self):
        return -1 * np.ones(len(self.implement_q_space))
  




class BrownianMotionExpUtilImplicit(SimpleLOB_Implementation):
    pass
        
class Test(unittest.TestCase):


    def testName(self):
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()