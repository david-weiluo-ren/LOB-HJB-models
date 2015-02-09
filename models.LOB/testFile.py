'''
Created on Feb 8, 2015

@author: weiluo
'''
from remoteExeFiles.SaveObj_helpers import constructObj_wrapper
from remoteExeFiles.SimulateSingle import dumpData, simulateImplicitSingle, prepareOptions
from remoteExeFiles.SimulateComparison import simulateImplicitComparison
if __name__ == '__main__':
    #options,  simulate_num, fileName = prepareOptions()
    #dumpData(simulateImplicitSingle(options,  simulate_num, fileName))
    simulateImplicitComparison()