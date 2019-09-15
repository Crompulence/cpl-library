import numpy as np
from shutil import move as mv 

class MOCK_Input(object):
    """

    Class to read and write into a Python MOCK script for either CFD or DEM/MD.

    """

    def __init__(self, input_file):
        """

        Any variable stored within 'USER INPUT - START' and 'USER INPUT - END' block
        are read in as attributes of the class. All other lines are ignored.

        """

        self.input_file = input_file
        with open(input_file, 'r') as f:
            input_block_start = False
            input_block_end = False
            for l in f:
                if 'USER INPUT - START' in l:
                    input_block_start = True

                if 'USER INPUT - END' in l:
                    input_block_end = True
                    break

                if input_block_start == True and '#' not in l and len(l.strip()) > 0:
                    try:
                        exec('self.' + l)
                    except:
                        print('Failed attempt to exec line: {}'.format(l))
                        pass

            # Check that the keywords were found within the file
            assert input_block_start, '{} not found in input file: {}'.format('USER INPUT - START', input_file)
            assert input_block_end, '{} not found in input file: {}'.format('USER INPUT - END', input_file)

            # Print extracted parameters
            print('Input parameters extracted from file {}'.format(input_file))
            print('\n'.join('{}: {}'.format(k, v) for k, v in self.__dict__.items()))


def MOCK_Writer(input_file, parameter, parameter_value):
    """

    Set the value for any of the input parameters within the 'USER INPUT -
    START' and 'USER INPUT - END' block.

    """

    with open(input_file, 'r') as old_file:
        input_block_start = False
        input_block_end = False
        parameter_updated = False
        with open(input_file + '.new', 'w') as new_file:
            for l in old_file:
                if 'USER INPUT - START' in l:
                    input_block_start = True

                if 'USER INPUT - END' in l:
                    input_block_end = True

                if input_block_start == True and input_block_end == False and '#' not in l and parameter in l:
                    parameter_updated = True
                    if isinstance(parameter_value, str):
                        new_file.write('{} = \'{}\'\n'.format(parameter, parameter_value))
                    else:
                        new_file.write('{} = {}\n'.format(parameter, parameter_value))
                else:
                    new_file.write(l)

    assert input_block_start, '{} not found in input file: {}'.format('USER INPUT - START', input_file)
    assert input_block_end, '{} not found in input file: {}'.format('USER INPUT - END', input_file)
    assert parameter_value, '{} not found in input file: {}'.format(parameter, input_file)

    mv(input_file + '.new', input_file)
