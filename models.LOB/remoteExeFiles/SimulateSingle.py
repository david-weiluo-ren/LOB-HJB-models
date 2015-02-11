'''
Created on Feb 8, 2015

@author: weiluo
'''

from remoteExeFiles.SaveObj_helpers import ImplicitMethodReader
import pickle
import re
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
    
    options = myReader.parserToArgsDict(parser)
    directory = options['dump_dir']
    options.pop('dump_dir')

    fileName = constructFileName(options, directory)
    simulate_num = options['simulate_num']
    options.pop('simulate_num')
    return [options,  simulate_num, fileName]

def summary_mean_var(options,simulate_num,fileName):
    myReader = ImplicitMethodReader()
    myObj= myReader.constructModelFromOptions_helper(options)

    mean_data = [0,0,0,0]  #[simulate_control_a_mean, simulate_control_b_mean, simulate_s_mean, simulate_q_mean]
    squared_data = [0, 0, 0, 0] #[simulate_control_a_squared, simulate_control_b_squared, simulate_s_squared, simulate_q_squared]
    
    myObj.run()
    print "done with run"
    
    successful_simulate_num = 0
    for _ in xrange(simulate_num):
        tmp_result = myObj.simulate_forward()
        if tmp_result[0]:
            successful_simulate_num += 1
            
            mean_data[0] += np.asarray(myObj.simulate_control_a)
            mean_data[1] += np.asarray(myObj.simulate_control_b)
            mean_data[2] += np.asarray(myObj.s)
            mean_data[3] += np.asarray(myObj.q)

            squared_data[0] += np.asarray(myObj.simulate_control_a)**2
            squared_data[1] += np.asarray(myObj.simulate_control_b)**2
            squared_data[2] += np.asarray(myObj.s)**2
            squared_data[3] += np.asarray(myObj.q)**2
    
    mean_data = [array/successful_simulate_num for array in mean_data]
    var_data = []
    for i in xrange(4):
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
             myObj.multi_fixed_q_control(myObj.q_space[0], myObj.q_space[-1], 1)], mean_data, var_data]
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


        self.simulated_a_control_var = self._loaded_data[4][0]
        self.simulated_b_control_var = self._loaded_data[4][1]
        self.simulated_s_var = self._loaded_data[4][2]
        self.simulated_q_var = self._loaded_data[4][3]
        
        
        
        
        self.successful_simulate_num = None
        


        
    

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