'''
Created on Feb 8, 2015

@author: weiluo
'''

from remoteExeFiles.SaveObj_helpers import ImplicitMethodReader
import pickle
import re, random
import numpy as np

def constructFileName(options, directory):
    print directory
    if len(options)==0:
        return "defaultParams"
    return directory + re.sub( r'[:,]',"_", re.sub(r'[\'{} ]', "", str(options)))


def prepareOptions():
    myReader = ImplicitMethodReader()
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
    
    
    options = myReader.parserToArgsDict(parser)
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
def summary_mean_var_helper(myObj, simulate_num, options, fileName, randomOpt):

    mean_data = [0,0,0,0,0]  #[simulate_control_a_mean, simulate_control_b_mean, simulate_s_mean, simulate_q_mean]
    squared_data = [0, 0, 0, 0,0] #[simulate_control_a_squared, simulate_control_b_squared, simulate_s_squared, simulate_q_squared]
    
    successful_simulate_num = 0
    q_0_origin = myObj.q_0
    q_0s = []
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
            
            squared_data[0] += np.asarray(myObj.simulate_control_a)**2
            squared_data[1] += np.asarray(myObj.simulate_control_b)**2
            squared_data[2] += np.asarray(myObj.s)**2
            squared_data[3] += np.asarray(myObj.q)**2
            try:
                squared_data[4] += np.asarray(myObj.s_drift)**2
            except:
                pass
    mean_data = [array/successful_simulate_num for array in mean_data]
    var_data = []
    for i in xrange(len(squared_data)):
        var_data.append(squared_data[i]/successful_simulate_num - mean_data[i]**2)        
   
   
    data_for_checking = [[], [], []]
    
    for i in np.arange(0, len(myObj.result),int(len(myObj.result)/20)):
        data_for_checking[0].append(myObj.result[i])
        data_for_checking[1].append(myObj.a_control[i])
        data_for_checking[2].append(myObj.b_control[i])

        
    print "done with the simulation"
    #return [fileName, myObj, mean_data, var_data]
    return [data_for_checking, fileName,\
            [options, successful_simulate_num, \
             myObj.multi_fixed_q_control(myObj.q_space[0], myObj.q_space[-1], 1)], mean_data, var_data, myObj.failed_simulation, q_0s]
def dumpData(data):
    fileHandler = open(data[0], 'w')
    pickle.dump(data, fileHandler)

class LoadSingleData(object):

    def __init__(self, loaded_data):
        self._loaded_data = loaded_data
        self.data_for_checking = self._loaded_data[0]
        self.result_for_checking = self.data_for_checking[0]
        self.a_control_for_checking = self.data_for_checking[1]
        self.b_control_for_checking = self.data_for_checking[2]
        
        
        self._fileName = self._loaded_data[1]
        
        self.options = self._loaded_data[2][0]
        self.successful_simulate_num = self._loaded_data[2][1]
        self.simulated_a_control = [arr[0] for arr in self._loaded_data[2][2]]
        self.simulated_b_control = [arr[1] for arr in self._loaded_data[2][2]]
        
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