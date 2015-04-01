'''
Created on Mar 5, 2015

@author: weiluo
'''
from remoteExeFiles.SaveObj_helpers import basicReader
from remoteExeFiles.SimulateSingle import constructFileName, prepareOptionsHelper,\
summary_mean_var_helper,prepareOptionsHelper2, dumpData
from remoteExeFiles.SimulateComparison import simulateImplicitComparison
from abstractLOB.poisson_OU_explicit import Poisson_explicit_OU_LOB
from abstractLOB.poisson_OU_implicit import Poisson_OU_implicit, Poisson_OU_implicit_truncateControlAtZero

def prepareOptions():
    myReader = basicReader()
    parser = prepareOptionsHelper(myReader)
    
    parser.add_argument('-alpha', type = float, nargs='?',\
                        help="mean-reverting rate")
    
    parser.add_argument('-half_S', type=float, nargs="?",\
                        help="radius of space we consider around s_0")
    
    parser.add_argument('-half_I_S', type=int, nargs='?',\
                        help="number of grid points of half of s space")
    parser.add_argument('-s_long_term_mean', type=float, nargs='?',\
                        help='long term mean of the OU price process')
    
    options = myReader.parserToArgsDict(parser)
    options.pop("type")
    options.pop("BC")
    return prepareOptionsHelper2(options)

def prepareOptions_forSameRandomness():
    myReader = basicReader()
    parser = prepareOptionsHelper(myReader)
    
    parser.add_argument('-alpha', type = float, nargs='?',\
                        help="mean-reverting rate")
    
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
    
    options = myReader.parserToArgsDict(parser)
    options.pop("type")
    options.pop("BC")
    return prepareOptionsHelper3(options)
def prepareOptionsHelper3(options):
    directory = options['dump_dir']
    options.pop('dump_dir')

    fileName = constructFileName(options, directory)
    simulate_num = options['simulate_num']
    options.pop('simulate_num')
    
    random_q_0  =  options['random_q_0'].upper()
    options.pop("random_q_0")
    options_forImplicit = options.copy()
    if "new_weight" in options:
        options.pop('new_weight')
    if "abs_threshold_power" in options:
        options.pop('abs_threshold_power')
    if 'rlt_threshold_power' in options:
        options.pop('rlt_threshold_power')
    options_forExplicit = options.copy()
    return [options_forImplicit, options_forExplicit,  simulate_num, fileName, random_q_0]

def summary_mean_var(options,simulate_num,fileName, randomOpt = False):
    myObj = Poisson_explicit_OU_LOB(**options)
    myObj.run()    
    print "done with run"
    return summary_mean_var_helper(myObj, simulate_num, options, fileName, randomOpt, False)

def simulateImplicitComparison_OU():
    options,  simulate_num, fileName, random_q_0 = prepareOptions()
    fileName += '_comparison'
    random_q_0_opt = False if random_q_0.upper()=="FALSE" else True
    if 'beta' in options and options['beta']==0.0:
        dumpData(summary_mean_var(options,  simulate_num, fileName, random_q_0_opt))
        return
    data = [fileName]
    nonZeroBetaOptions = options.copy()
    zeroBetaOptions = options.copy()
    data.append(summary_mean_var(nonZeroBetaOptions,  simulate_num, fileName, random_q_0_opt))
    zeroBetaOptions['beta'] = 0.0
    data.append(summary_mean_var(zeroBetaOptions,  simulate_num, fileName, random_q_0_opt))
    
    dumpData(data)
    
def simulateImplicitComparison_OU_obj():
    options,  simulate_num, fileName, random_q_0 = prepareOptions()
    fileName += '_comparison'
    random_q_0_opt = False if random_q_0.upper()=="FALSE" else True
    if 'beta' in options and options['beta']==0.0:
        myObj = Poisson_OU_implicit(**options)
        myObj.run()  
        myObj.simulate_forward()
        dumpData([fileName, myObj])
        return
    data = [fileName]
    nonZeroBetaOptions = options.copy()
    zeroBetaOptions = options.copy()
    myObj = Poisson_OU_implicit(**nonZeroBetaOptions)
    myObj.run()  
    myObj.simulate_forward()
    
    data.append(myObj)
    zeroBetaOptions['beta'] = 0.0
    myObj = Poisson_OU_implicit(**zeroBetaOptions)
    myObj.run()  
    myObj.simulate_forward()
    data.append(myObj)
    
    dumpData(data)
    
def simulateComparison_OU_sameRandomness():
    options_forImplicit, options_forExplicit,  simulate_num, fileName, random_q_0 = prepareOptions_forSameRandomness()
    fileName += '_comparison'
    random_q_0_opt = False if random_q_0.upper()=="FALSE" else True
   
    data = [fileName]
    
    myObjImplicit = Poisson_OU_implicit(**options_forImplicit)
    myObjExplicit = Poisson_explicit_OU_LOB(**options_forExplicit)

    myObjImplicit.run()
    myObjExplicit.run()
    for _ in xrange(simulate_num):
        randomSource = myObjExplicit.generate_random_source()
        myObjExplicit.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        myObjImplicit.simulate_forward(useGivenRandom=True, randomSource=randomSource)
        tmpdata_explicit = [myObjExplicit.s, myObjExplicit.q, myObjExplicit.x, myObjExplicit.simulate_control_a, myObjExplicit.simulate_control_b]
        tmpdata_implicit = [myObjImplicit.s, myObjImplicit.q, myObjImplicit.x, myObjImplicit.simulate_control_a, myObjImplicit.simulate_control_b, \
                            myObjImplicit.simulate_price_a, myObjImplicit.simulate_price_b]
        data.append([tmpdata_explicit, tmpdata_implicit])

    
    dumpData(data)