'''
Created on Feb 8, 2015

@author: weiluo
'''


import SaveObj_helpers
if __name__ == '__main__':
    options = SaveObj_helpers.parserToArgsDict(SaveObj_helpers.basicParser_forAbstractLOB())

    model_type = options['type'].upper()
    options.pop('type')
    BC_type = options['BC'].upper()
    options.pop('BC')

    