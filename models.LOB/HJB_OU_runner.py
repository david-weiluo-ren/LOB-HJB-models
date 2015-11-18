from abstractLOB.HJB_OU import HJB_OU_solver
import pickle
import abstractLOB
import argparse
import re
from sys import argv
import numpy as np
'''
Created on Sep 9, 2015

@author: weiluo
'''



def dumpData(data, fileName=None):
    if fileName is None:
        fileHandler = open(data[0], 'w')
    else:
        fileHandler = open(fileName, 'w')
    pickle.dump(data, fileHandler)

def run_HJB_OU_solver(PC_exact=False):
    data, sample_stepSize  = save_OU_sampleValueFunction_helper(PC_exact)
    dump_data = [data[0], data[1]]
    myObjImplicit_no_truncation = data[2]
    myObjImplicit_no_truncation_unrun = data[3]
    sample_result = myObjImplicit_no_truncation.value_function[::sample_stepSize]
    sample_a_price = myObjImplicit_no_truncation._a_price[::sample_stepSize]
    sample_b_price = myObjImplicit_no_truncation._b_price[::sample_stepSize]
    dump_data.append([sample_stepSize, myObjImplicit_no_truncation_unrun, sample_result, sample_a_price, sample_b_price])
    dumpData(dump_data)
    
def save_OU_sampleValueFunction_helper(PC_exact):
    options_forImplicit, fileName, sample_stepSize, normalization\
    = prepareOptions_forSaveSampleValueFunction()
    fileName += '_obj'
    data = [fileName[:200], options_forImplicit]
    myObjImplicit_no_truncation = HJB_OU_solver(**options_forImplicit)
    if PC_exact:
        myObjImplicit_no_truncation.run_PC(exact=True, normalization=normalization)
    else:
        myObjImplicit_no_truncation.run__OU_PC_log_hybrid()
    myObjImplicit_no_truncation_unrun = HJB_OU_solver(**options_forImplicit)
    myObjImplicit_no_truncation_unrun.linear_system = None
    if PC_exact and (not normalization):
        myObjImplicit_no_truncation.value_function[:] = -1/ myObjImplicit_no_truncation_unrun.gamma * np.log(myObjImplicit_no_truncation.value_function_PC)
    data.append(myObjImplicit_no_truncation)
    data.append(myObjImplicit_no_truncation_unrun)
    return data, sample_stepSize

def prepareOptions_forSaveSampleValueFunction():
    parser = argparse.ArgumentParser()
    _sample_stepSize = 100
    _dump_dir = ""

    parser.add_argument('-A', type = float, nargs='?',\
                        help = "A")
    parser.add_argument('-gamma', type = float, nargs='?',\
                        help = "gamma")
    parser.add_argument('-kappa', type = float, nargs='?',\
                        help = "kappa")
    parser.add_argument('-half_I', type = int, nargs='?',\
                        help = "half_I")
    parser.add_argument('-sigma_s', type = float, nargs='?',\
                        help = "sigma_s")
    parser.add_argument('-delta_t', type = float, nargs='?',\
                        help = "delta_t")
    parser.add_argument('-verbose', type = bool, nargs='?',\
                        help = "verbose option")
    parser.add_argument('-num_time_step', type = int, nargs='?',\
                        help = "num_time_step")
    parser.add_argument('-extend_space', type = float, nargs='?',\
                        help = "extend_space with the same unit of N")
    parser.add_argument('-use_sparse', type = bool, nargs = '?',\
                        help = 'If true, then the program will use sparse matrix when solving the linear system.')

    parser.add_argument('-sample_stepSize', type=int, nargs='?',default = _sample_stepSize,\
                        help="the step size for sampling the value function")
    parser.add_argument('-alpha', type = float, nargs='?',\
                        help="mean-reverting rate")
    parser.add_argument('-lambda_tilde', type = float, nargs='?',\
                        help="left inventory penalty")
    parser.add_argument('-half_S', type=float, nargs="?",\
                        help="radius of space we consider around s_0")
    
    parser.add_argument('-half_I_S', type=int, nargs='?',\
                        help="number of grid points of half of s space")
    parser.add_argument('-s_long_term_mean', type=float, nargs='?',\
                        help='long term mean of the OU price process')
    parser.add_argument('-new_weight', type=float, nargs='?',\
                        help='weight of new iteration result in the iteration')
    parser.add_argument('-abs_threshold_power', type=float, nargs='?',\
                        help="the power of absolute threshold for determining when to stop the iteration")
    parser.add_argument('-rlt_threshold_power', type=float, nargs='?',\
                        help="the power of relative threshold for determining when to stop the iteration")
    parser.add_argument('-iter_max', type=float, nargs='?',\
                        help="max number of iteration")
    parser.add_argument('-gueant_boundary', type=bool, nargs='?',\
                        help="whether to use gueant's boundary")
    parser.add_argument("-dump_dir", type = str, default = _dump_dir,\
                        nargs = '?', help="the directory containing the dumped objs")
    parser.add_argument("-boundary_factor", type = float, \
                        nargs = '?', help="factor used in the boundary of the zero 2nd-derivative boundary method.")
    parser.add_argument("-quadratic_boundary_factor", type = float, \
                        nargs = '?', help="quadratic factor used in the boundary on the second derivative of value function")
    parser.add_argument("-data_storing_jump_size", type = float, \
                        nargs = '?', help="jump size for storing data.")
    parser.add_argument("-OU_step", type=int, \
                       nargs='?', help="number of OU in hybrid run")

    parser.add_argument("-normalization", type=bool, default=False, nargs='?', help='normalization for run_PC')
    parser.add_argument("-record_time_lower_bound", type=int, nargs='?')

    parser.add_argument("PC_exact", nargs='?')
    
    

    options = parserToArgsDict(parser)
    directory = options['dump_dir']
    options.pop('dump_dir')
    if 'PC_exact' in options:
        options.pop('PC_exact')
    normalization = True
    if 'normalization' in options:
        normalization = options['normalization']
        options.pop('normalization')
    fileName = constructFileName(options, directory)

    sample_stepSize = options["sample_stepSize"]
    options.pop("sample_stepSize")
    return [options, fileName, sample_stepSize, normalization]



def parserToArgsDict(parser):
    opts = vars(parser.parse_args())
    options = {k : opts[k] for k in opts if opts[k] is not None}
    return options

def constructFileName(options, directory):
    print directory
    if len(options)==0:
        return "defaultParams"
    return directory + re.sub( r'[:,]',"_", re.sub(r'[\'{} ]', "", str(options)))

if __name__ == '__main__':
    if len(argv) > 1 and argv[1] == 'PC_exact':
        print 'Use PC exact and normalization'
        run_HJB_OU_solver(PC_exact=True)
    else:
        run_HJB_OU_solver()
    


