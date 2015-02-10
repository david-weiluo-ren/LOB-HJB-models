'''
Created on Feb 8, 2015

@author: weiluo
'''
import sys
import pickle
import argparse
from abstractLOB.brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC, BrownianMotion_ExpUtil_Implicit_sameSlopeBC
from abstractLOB.poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC, Poisson_expUtil_implicit_sameSlopeBC
                                                
class basicReader(object):
    @classmethod
    def argParser(cls):
        parser = argparse.ArgumentParser()
        _type = "BM" #It's either "BM" or "Poisson"
        _BC = "Neumann" #Currently support "Neumann" and "sameSlope"
        _gamma = None
        _kappa = None
        _beta = None
        _N = None
        _half_I = None
        _sigma_s = None
        _q_0 = None
        _x_0 = None
        _s_0 = None
        _delta_t = None
        _verbose = None
        _num_time_step = None
        _extend_space = None
        _expIntensity_lower_tolerance_power = None  
          
        parser.add_argument('-type', type = str, default = _type, nargs='?', \
                            help="The type of model, either 'BM' or 'Poisson'")
        parser.add_argument('-BC', type = str, default = _BC, nargs = '?',\
                            help = "The type of boundary condition. \
                            Currently supporting 'Neumann' and 'sameSlope'")
        parser.add_argument('-A', type = float, nargs='?',\
                            help = "A")
        parser.add_argument('-gamma', type = float, nargs='?',\
                            help = "gamma")
        parser.add_argument('-kappa', type = float, nargs='?',\
                            help = "kappa")
        parser.add_argument('-beta', type = float, nargs='?',\
                            help = "beta")
        parser.add_argument('-N', type = float, nargs='?',\
                            help = "N")
        parser.add_argument('-half_I', type = int, nargs='?',\
                            help = "half_I")
        parser.add_argument('-sigma_s', type = float, nargs='?',\
                            help = "sigma_s")
        parser.add_argument('-q_0', type = float, nargs='?',\
                            help = "q_0")
        parser.add_argument('-x_0', type = float, nargs='?',\
                            help = "x_0")
        parser.add_argument('-s', type = float, nargs='?',\
                            help = "s_0")
        parser.add_argument('-delta_t', type = float, nargs='?',\
                            help = "delta_t")
        parser.add_argument('-verbose', type = bool, nargs='?',\
                            help = "verbose option")
        parser.add_argument('-num_time_step', type = int, nargs='?',\
                            help = "num_time_step")
        parser.add_argument('-extend_space', type = float, nargs='?',\
                            help = "extend_space with the same unit of N")
        parser.add_argument('-expIntensity_lower_tolerance_power', type = float, nargs='?',\
                            help = "expIntensity_lower_tolerance_power")
        
        
        
        return parser
    
    @classmethod
    def parserToArgsDict(cls, parser):
        opts = vars(parser.parse_args())
        if 'cmd' in opts:
            command = opts.pop('cmd')
        options = {k : opts[k] for k in opts if opts[k] is not None}
        return options



class ImplicitMethodReader(basicReader):
    
    def __init__(self):
        self.all_models = {\
                      ('BM', 'NEUMANN'): BrownianMotion_ExpUtil_Implicit_NeumannBC,\
                      ('BM', 'SAMESLOPE'): BrownianMotion_ExpUtil_Implicit_sameSlopeBC,\
                      ('POISSON', 'NEUMANN'): Poisson_expUtil_implicit_NeumannBC,\
                      ('POISSON', 'SAMESLOPE'): Poisson_expUtil_implicit_sameSlopeBC\
                      }    
    
    @classmethod
    def argParser(cls):
        parser = super(ImplicitMethodReader, cls).argParser()
        parser.add_argument('-iter_max', type = int, nargs='?',\
                            help = "the maximal number of iteration at each time")
        parser.add_argument('-new_weight', type = float, nargs = '?',\
                            help = "The weight on the new number when doing iteration")
        parser.add_argument('-abs_threshold_power', type = float, nargs = '?',\
                            help = "The power for the absolute threshold when deciding whether\
                            to terminate the iteration. The absolute threshold would be 10**abs_threshold_power,\
                            so usually such argument should be negative.")
        
        parser.add_argument('-rlt_threshold_power', type = float, nargs = '?',\
                            help = "The power for the relative threshold when deciding whether\
                            to terminate the iteration. The relative threshold would be 10**rlt_threshold_power,\
                            so usually such argument should be negative.")
        parser.add_argument('-use_sparse', type = bool, nargs = '?',\
                            help = 'If true, then the program will use sparse matrix when solving the linear system.')
      
        return parser
    
    def constructModel(self):
        options = self.parserToArgsDict(self.argParser())
        return self.constructModelFromOptions_helper(options)
    
    def constructModelFromOptions_helper(self, options):
        options_local = options.copy()
        model_type = options_local['type'].upper()
        options_local.pop('type')
        BC_type = options_local['BC'].upper()
        options_local.pop('BC')
        return self.all_models[(model_type, BC_type)](**options_local)
       
def constructObj_wrapper():
    myReader = ImplicitMethodReader()
    myObj = myReader.constructModel()
    return myObj
if __name__ == '__main__':
    pass 
    
    
    
    
    