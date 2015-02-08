'''
Created on Feb 8, 2015

@author: weiluo
'''
import sys
import pickle
import argparse
from abstractLOB.brownianMotion_expUtil_implicit import BrownianMotion_ExpUtil_Implicit_NeumannBC, BrownianMotion_ExpUtil_Implicit_sameSlopeBC
from abstractLOB.poisson_expUtil_implicit import Poisson_expUtil_implicit_NeumannBC, Poisson_expUtil_implicit_sameSlopeBC
                                                        

def basicParser_forAbstractLOB():
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
    
    
def parserToArgsDict(parser):
    opts = vars(parser.parse_args())
    if 'cmd' in opts:
        command = opts.pop('cmd')
    options = {k : opts[k] for k in opts if opts[k] is not None}
    return options
    
class BasicReader:
    all_models = {\
                      ('BM', 'NEUMANN'): BrownianMotion_ExpUtil_Implicit_NeumannBC,\
                      ('BM', 'SAMESLOPE'): BrownianMotion_ExpUtil_Implicit_sameSlopeBC,\
                      ('POISSON', 'NEUMANN'): Poisson_expUtil_implicit_NeumannBC,\
                      ('POISSON', 'SAMESLOPE'): Poisson_expUtil_implicit_sameSlopeBC\
                      }    
    
    def constructModel(self):
        options = parserToArgsDict(basicParser_forAbstractLOB())
        model_type = options['type'].upper()
        options.pop('type')
        BC_type = options['BC'].upper()
        options.pop('BC')
        return self.all_models[(model_type, BC_type)](options)
        
       
def constructObj_wrapper():
    myReader = BasicReader()
    myObj = myReader.constructModel()
    return myObj
if __name__ == '__main__':
    pass 
    
    
    
    
    