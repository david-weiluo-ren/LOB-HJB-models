'''
Created on Feb 8, 2015

@author: weiluo
'''

from remoteExeFiles.SaveObj_helpers import ImplicitMethodReader
import pickle
import re
def constructFileName(options, directory):
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
    parser.add_argument(_dump_dir, type = str, default = _simulate_num,\
                        nargs = '?', help="the directory containing the dumped objs")
    
    options = myReader.parserToArgsDict(parser)
    fileName = constructFileName(options, options['dump_dir'])
    options.pop('dump_dir')
    simulate_num = options['simulate_num']
    options.pop('simulate_num')
    return [options,  simulate_num, fileName]

def simulateImplicitSingle(options,simulate_num,fileName):
    myReader = ImplicitMethodReader()
    myObj = myReader.constructModelFromOptions_helper(options)
    simulate_control_a = []
    simulate_control_b = []
    simulate_s = []
    simulate_q = []
    myObj.run()
    for _ in xrange(simulate_num):
        tmp_result = myObj.simulate_forward()
        if tmp_result[0]:
            simulate_control_a.append(myObj.simulate_control_a)
            simulate_control_b.append(myObj.simulate_control_b)
            simulate_s.append(myObj.s)
            simulate_q.append(myObj.q)
    
    return [fileName, simulate_control_a, simulate_control_b, simulate_s, simulate_q]
    

def dumpData(data):
    fileHandler = open(data[0], 'w')
    pickle.dump(data, fileHandler)
    
    

if __name__ == '__main__':
    pass