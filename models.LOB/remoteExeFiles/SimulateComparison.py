'''
Created on Feb 8, 2015

@author: weiluo
'''
from remoteExeFiles.SimulateSingle import prepareOptions, dumpData, summary_mean_var, LoadSingleData
from remoteExeFiles.SaveObj_helpers import ImplicitMethodReader
def simulateImplicitComparison():
    options,  simulate_num, fileName = prepareOptions()
    fileName += '_comparison'
    if 'beta' in options and options['beta']==0.0:
        dumpData(summary_mean_var(options,  simulate_num, fileName))
        return
    data = [fileName]
    nonZeroBetaOptions = options.copy()
    zeroBetaOptions = options.copy()
    data.append(summary_mean_var(nonZeroBetaOptions,  simulate_num, fileName))
    zeroBetaOptions['beta'] = 0.0
    data.append(summary_mean_var(zeroBetaOptions,  simulate_num, fileName))
    
    dumpData(data)

class LoadMultiData(object):
    def __init__(self, loaded_data):
        self._loaded_data = loaded_data
        self.non_zero_beta_obj = LoadSingleData(self._loaded_data[1])
        self.zero_beta_obj = LoadSingleData(self._loaded_data[2])
        self.reader = ImplicitMethodReader()

    def constructNonZeroBetaObj(self):
        return self.reader.constructModelFromOptions_helper(self.non_zero_beta_obj.options)
    
    def constructZeroBetaObj(self):
        return self.reader.constructModelFromOptions_helper(self.zero_beta_obj.options)

        
    
    
    
    
    
    

if __name__ == '__main__':
    pass