'''
Created on Feb 8, 2015

@author: weiluo
'''

from remoteExeFiles.SaveObj_helpers import ImplicitMethodReader
import pickle
import re, random
import numpy as np
from pylab import plot
from Cython.Shadow import NULL
def constructFileName(options, directory):
    print directory
    if len(options)==0:
        return "defaultParams"
    return directory + re.sub( r'[:,]',"_", re.sub(r'[\'{} ]', "", str(options)))


def prepareOptions(_myReader=None):
    
    myReader = ImplicitMethodReader() if _myReader is None else _myReader
    
    parser = prepareOptionsHelper(myReader)
    options = myReader.parserToArgsDict(parser)
    return prepareOptionsHelper2(options)

def prepareOptionsHelper(myReader):
    parser = myReader.argParser()
    _simulate_num = 1000
    parser.add_argument('-simulate_num', type = int, default = _simulate_num,\
                        nargs = '?', help="number of trajectories to simulate")
    
    _dump_dir = ""
    parser.add_argument("-dump_dir", type = str, default = _dump_dir,\
                        nargs = '?', help="the directory containing the dumped objs")
    
    _random_q_0 = "FALSE"
    parser.add_argument("-random_q_0", type = str, default = _random_q_0,\
                        nargs = '?', help="whether or not to use random q_0")
    
    
    return parser
    

def prepareOptionsHelper2(options):

    directory = options['dump_dir']
    options.pop('dump_dir')

    fileName = constructFileName(options, directory)
    simulate_num = options['simulate_num']
    options.pop('simulate_num')
    
    random_q_0  =  options['random_q_0'].upper()
    options.pop("random_q_0")
    return [options,  simulate_num, fileName, random_q_0]


def summary_mean_var(options,simulate_num,fileName, randomOpt = False):
    myReader = ImplicitMethodReader()
    myObj= myReader.constructModelFromOptions_helper(options)

     
    myObj.run()
    print "done with run"
    return summary_mean_var_helper(myObj, simulate_num, options, fileName, randomOpt)
def mean_var_cmp(myObj, simulate_num=1000, q_0=0):
    
    
    mean_data = [0,0,0,0,0,0,0]  #[simulate_control_a_mean, simulate_control_b_mean, simulate_s_mean, simulate_q_mean]
    squared_data = [0, 0, 0, 0,0,0,0] #[simulate_control_a_squared, simulate_control_b_squared, simulate_s_squared, simulate_q_squared]
    successful_simulate_num = 0
    for _ in xrange(simulate_num):
        tmp_result = myObj.simulate_forward(q_0=q_0)

        if tmp_result[0]:
            successful_simulate_num += 1
            
            mean_data[0] += np.asarray(myObj.simulate_control_a)
            mean_data[1] += np.asarray(myObj.simulate_control_b)
            mean_data[2] += np.asarray(myObj.s)
            mean_data[3] += np.asarray(myObj.q)
            
            try:
                mean_data[4] += np.asarray(myObj.s_drift)
            except:
                pass
            try: 
                mean_data[5] += np.asarray(myObj.simulate_price_a)
                mean_data[6] += np.asarray(myObj.simulate_price_b)
            except:
                pass
            squared_data[0] += np.asarray(myObj.simulate_control_a)**2
            squared_data[1] += np.asarray(myObj.simulate_control_b)**2
            squared_data[2] += np.asarray(myObj.s)**2
            squared_data[3] += np.asarray(myObj.q)**2
            try:
                squared_data[4] += np.asarray(myObj.s_drift)**2
            except:
                pass
            try: 
                squared_data[5] += np.asarray(myObj.simulate_price_a)**2
                squared_data[6] += np.asarray(myObj.simulate_price_b)**2
            except:
                pass
            
            #tmp_data.append( np.asarray(myObj.s).std())
    mean_data = [np.true_divide(array,successful_simulate_num) for array in mean_data]
    var_data = []
    for i in xrange(len(squared_data)):
        var_data.append(np.true_divide(squared_data[i],successful_simulate_num) - mean_data[i]**2)        
        
    print "done with the simulation"
    #return [fileName, myObj, mean_data, var_data]
    
    return[mean_data, var_data, myObj.failed_simulation]

def summary_mean_var_helper(myObj, simulate_num, options, fileName, randomOpt, dataCheckingOption=False):

    mean_data = [0,0,0,0,0,0,0,0,0]  #[simulate_control_a_mean, simulate_control_b_mean, simulate_s_mean, simulate_q_mean]
    squared_data = [0, 0, 0, 0,0,0,0,0,0] #[simulate_control_a_squared, simulate_control_b_squared, simulate_s_squared, simulate_q_squared]
    
    successful_simulate_num = 0
    q_0_origin = myObj.q_0
    q_0s = []
    tmp_data = []
    for _ in xrange(simulate_num):
        if not randomOpt:
            tmp_result = myObj.simulate_forward()
        else:
            q_0_random = 2*random.random()*q_0_origin - q_0_origin
            tmp_result = myObj.simulate_forward(q_0 = q_0_random)
            q_0s.append(q_0_random)
        if tmp_result[0]:
            successful_simulate_num += 1
            
            mean_data[0] += np.asarray(myObj.simulate_control_a)
            mean_data[1] += np.asarray(myObj.simulate_control_b)
            mean_data[2] += np.asarray(myObj.s)
            mean_data[3] += np.asarray(myObj.q)
            
            try:
                mean_data[4] += np.asarray(myObj.s_drift)
            except:
                pass
            try: 
                mean_data[5] += np.asarray(myObj.simulate_price_a)
                mean_data[6] += np.asarray(myObj.simulate_price_b)
                mean_data[7] += np.asarray(myObj.x) + np.asarray(myObj.s) * np.asarray(myObj.q)
            except:
                pass
            try:
                mean_data[8] += np.asarray(myObj.x) + (np.asarray(myObj.s) - myObj.lambda_tilde * np.asarray(myObj.q))* np.asarray(myObj.q)
            except:
                pass

            squared_data[0] += np.asarray(myObj.simulate_control_a)**2
            squared_data[1] += np.asarray(myObj.simulate_control_b)**2
            squared_data[2] += np.asarray(myObj.s)**2
            squared_data[3] += np.asarray(myObj.q)**2
            try:
                squared_data[4] += np.asarray(myObj.s_drift)**2
            except:
                pass
            try: 
                squared_data[5] += np.asarray(myObj.simulate_price_a)**2
                squared_data[6] += np.asarray(myObj.simulate_price_b)**2
                squared_data[7] += (np.asarray(myObj.x) + np.asarray(myObj.s) * np.asarray(myObj.q))**2
            except:
                pass
            try:
                squared_data[8] += (np.asarray(myObj.x) + (np.asarray(myObj.s) - myObj.lambda_tilde * np.asarray(myObj.q))* np.asarray(myObj.q))**2
            except:
                pass            
            #tmp_data.append( np.asarray(myObj.s).std())
            

            
    mean_data = [np.true_divide(array,successful_simulate_num) for array in mean_data]
    var_data = []
    for i in xrange(len(squared_data)):
        var_data.append(np.true_divide(squared_data[i],successful_simulate_num) - mean_data[i]**2)        
    print mean_data[4]
    print squared_data[4]
    data_for_checking = [[], [], []]
    if dataCheckingOption:
        for i in np.arange(0, len(myObj.value_function),int(len(myObj.value_function)/20)):
            data_for_checking[0].append(myObj.value_function[i])
            data_for_checking[1].append(myObj.a_control[i])
            data_for_checking[2].append(myObj.b_control[i])

    fixed_q_control=[options, successful_simulate_num, \
             myObj.multi_fixed_q_control(myObj.q_space[0], myObj.q_space[-1], 1)]\
             if dataCheckingOption else [options, successful_simulate_num,[[], []]]
    print "done with the simulation"
    #return [fileName, myObj, mean_data, var_data]
    return [data_for_checking, fileName,\
            fixed_q_control, mean_data, var_data, myObj.failed_simulation, q_0s, tmp_data]
def dumpData(data, fileName=None):
    if fileName is None:
        fileHandler = open(data[0], 'w')
    else:
        fileHandler = open(fileName, 'w')
    pickle.dump(data, fileHandler)

class LoadSingleData(object):

    def __init__(self, loaded_data):
        self._loaded_data = loaded_data
        self.data_for_checking = self._loaded_data[0]
        self.result_for_checking = self.data_for_checking[0]
        self.a_control_for_checking = self.data_for_checking[1]
        self.b_control_for_checking = self.data_for_checking[2]
        
        
        self._fileName = self._loaded_data[1]
        try:
            self.options = self._loaded_data[2][0]
            self.successful_simulate_num = self._loaded_data[2][1]
            self.simulated_a_control = [arr[0] for arr in self._loaded_data[2][2]]
            self.simulated_b_control = [arr[1] for arr in self._loaded_data[2][2]]
        except:
            self.options = self._loaded_data[2][0]
            self.successful_simulate_num = self._loaded_data[2][1]
            self.simulated_a_control = []
            self.simulated_b_control = []
        self.simulated_a_control_mean = self._loaded_data[3][0]
        self.simulated_b_control_mean = self._loaded_data[3][1]
        self.simulated_s_mean = self._loaded_data[3][2]
        self.simulated_q_mean = self._loaded_data[3][3]
        try:
            self.simulated_s_drift_mean = self._loaded_data[3][4]
        except:
            self.simulated_s_drift_mean = None

        self.simulated_a_control_var = self._loaded_data[4][0]
        self.simulated_b_control_var = self._loaded_data[4][1]
        self.simulated_s_var = self._loaded_data[4][2]
        self.simulated_q_var = self._loaded_data[4][3]
        try:
            self.simulated_s_drift_var = self._loaded_data[4][4]
        except:
            self.simulated_s_drift_var = None
        
        
        
        
        try:
            self.failed_simulation = self._loaded_data[5]
        except:
            self.failed_simulation = None
        

def summary_mean_var_constantPrice_helper(myObj, simulate_num, options, fileName, randomOpt, dataCheckingOption=False):

    mean_data_optimal = [0,0,0,0,0,0,0,0,0]  #[simulate_control_a_mean, simulate_control_b_mean, simulate_s_mean, simulate_q_mean]
    squared_data_optimal = [0, 0, 0, 0,0,0,0,0,0] #[simulate_control_a_squared, simulate_control_b_squared, simulate_s_squared, simulate_q_squared]
    
    mean_data_constant = [0,0,0,0,0,0,0,0,0] 
    squared_data_constant = [0, 0, 0, 0,0,0,0,0,0]
    
    
    
    successful_simulate_num = 0
    q_0_origin = myObj.q_0
    q_0s = []
    tmp_data = []
    for _ in xrange(simulate_num):
        randomness = myObj.generate_random_source()
        tmp_result = myObj.simulate_forward(useGivenRandom = True, randomSource = randomness)

        if tmp_result[0]:
            successful_simulate_num += 1
            
            mean_data_optimal[0] += np.asarray(myObj.simulate_control_a)
            mean_data_optimal[1] += np.asarray(myObj.simulate_control_b)
            mean_data_optimal[2] += np.asarray(myObj.s)
            mean_data_optimal[3] += np.asarray(myObj.q)
            
            try:
                mean_data_optimal[4] += np.asarray(myObj.s_drift)
            except:
                pass
            try: 
                mean_data_optimal[5] += np.asarray(myObj.simulate_price_a)
                mean_data_optimal[6] += np.asarray(myObj.simulate_price_b)
                mean_data_optimal[7] += np.asarray(myObj.x) + np.asarray(myObj.s) * np.asarray(myObj.q)
            except:
                pass
            try:
                mean_data_optimal[8] += np.asarray(myObj.x) + (np.asarray(myObj.s) - myObj.lambda_tilde * np.asarray(myObj.q))* np.asarray(myObj.q)
            except:
                pass

            squared_data_optimal[0] += np.asarray(myObj.simulate_control_a)**2
            squared_data_optimal[1] += np.asarray(myObj.simulate_control_b)**2
            squared_data_optimal[2] += np.asarray(myObj.s)**2
            squared_data_optimal[3] += np.asarray(myObj.q)**2
            try:
                squared_data_optimal[4] += np.asarray(myObj.s_drift)**2
            except:
                pass
            try: 
                squared_data_optimal[5] += np.asarray(myObj.simulate_price_a)**2
                squared_data_optimal[6] += np.asarray(myObj.simulate_price_b)**2
                squared_data_optimal[7] += (np.asarray(myObj.x) + np.asarray(myObj.s) * np.asarray(myObj.q))**2
            except:
                pass
            try:
                squared_data_optimal[8] += (np.asarray(myObj.x) + (np.asarray(myObj.s) - myObj.lambda_tilde * np.asarray(myObj.q))* np.asarray(myObj.q))**2
            except:
                pass            
            #tmp_data.append( np.asarray(myObj.s).std())
            tmp_result = myObj.simulate_forward_constantPrice(useGivenRandom = True, randomSource = randomness)
            mean_data_constant[0] += np.asarray(myObj.simulate_control_a)
            mean_data_constant[1] += np.asarray(myObj.simulate_control_b)
            mean_data_constant[2] += np.asarray(myObj.s)
            mean_data_constant[3] += np.asarray(myObj.q)
            
            try:
                mean_data_constant[4] += np.asarray(myObj.s_drift)
            except:
                pass
            try: 
                mean_data_constant[5] += np.asarray(myObj.simulate_price_a)
                mean_data_constant[6] += np.asarray(myObj.simulate_price_b)
                mean_data_constant[7] += np.asarray(myObj.x) + np.asarray(myObj.s) * np.asarray(myObj.q)
            except:
                pass
            try:
                mean_data_constant[8] += np.asarray(myObj.x) + (np.asarray(myObj.s) - myObj.lambda_tilde * np.asarray(myObj.q))* np.asarray(myObj.q)
            except:
                pass

            squared_data_constant[0] += np.asarray(myObj.simulate_control_a)**2
            squared_data_constant[1] += np.asarray(myObj.simulate_control_b)**2
            squared_data_constant[2] += np.asarray(myObj.s)**2
            squared_data_constant[3] += np.asarray(myObj.q)**2
            try:
                squared_data_constant[4] += np.asarray(myObj.s_drift)**2
            except:
                pass
            try: 
                squared_data_constant[5] += np.asarray(myObj.simulate_price_a)**2
                squared_data_constant[6] += np.asarray(myObj.simulate_price_b)**2
                squared_data_constant[7] += (np.asarray(myObj.x) + np.asarray(myObj.s) * np.asarray(myObj.q))**2
            except:
                pass
            try:
                squared_data_constant[8] += (np.asarray(myObj.x) + (np.asarray(myObj.s) - myObj.lambda_tilde * np.asarray(myObj.q))* np.asarray(myObj.q))**2
            except:
                pass 
            
    mean_data_optimal = [np.true_divide(array,successful_simulate_num) for array in mean_data_optimal]
    
    var_data_optimal = []
    
    mean_data_constant = [np.true_divide(array,successful_simulate_num) for array in mean_data_constant]
    var_data_constant = []
    
    for i in xrange(len(squared_data_optimal)):
        var_data_optimal.append(np.true_divide(squared_data_optimal[i],successful_simulate_num) - mean_data_optimal[i]**2)     
        
    for i in xrange(len(squared_data_constant)):
        var_data_constant.append(np.true_divide(squared_data_constant[i],successful_simulate_num) - mean_data_constant[i]**2)     
         
    print mean_data_optimal[4]
    print squared_data_optimal[4]
    data_for_checking = [[], [], []]
   
    fixed_q_control=[options, successful_simulate_num, \
             myObj.multi_fixed_q_control(myObj.q_space[0], myObj.q_space[-1], 1)]\
             if dataCheckingOption else [options, successful_simulate_num,[[], []]]
    print "done with the simulation"
    #return [fileName, myObj, mean_data_optimal, var_data_optimal]
    return [data_for_checking, fileName,\
            fixed_q_control, mean_data_optimal, var_data_optimal, mean_data_constant, var_data_constant, myObj.failed_simulation]        
    

if __name__ == '__main__':
    pass

"""
Really a bad idea. The file created is too large. Have to compress on the air.
def simulateImplicitSingle(options,simulate_num,fileName):
    myReader = ImplicitMethodReader()
    myObj = myReader.constructModelFromOptions_helper(options)
    simulate_control_a = []
    simulate_control_b = []
    simulate_s = []
    simulate_q = []
    myObj.run()
    print "done with run"
    for _ in xrange(simulate_num):
        tmp_result = myObj.simulate_forward()
        if tmp_result[0]:
            simulate_control_a.append(myObj.simulate_control_a)
            simulate_control_b.append(myObj.simulate_control_b)
            simulate_s.append(myObj.s)
            simulate_q.append(myObj.q)
    print "done with the simulation"
    return [fileName, myObj, simulate_control_a, simulate_control_b, simulate_s, simulate_q]
"""