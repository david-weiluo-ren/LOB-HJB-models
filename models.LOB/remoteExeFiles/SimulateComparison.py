'''
Created on Feb 8, 2015

@author: weiluo
'''
from remoteExeFiles.SimulateSingle import prepareOptions, dumpData, summary_mean_var
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
if __name__ == '__main__':
    pass