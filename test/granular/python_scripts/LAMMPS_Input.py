from shutil import move as mv 

class LAMMPS_Input(object):
    """

    Class to read and write into LAMMPS input script.
    
    """

    def __init__(self, input_file):
        """ 

        Initialise with all the variables in the variable block at the top of
        the input script. Note that this terminates as soon as the 'units'
        command is found, so any variables found afterward will not be found.
        Use more general 'read_script' if variables or properties afterwards
        are required.

        """
        self.input_file = input_file

        with open(self.input_file, 'r') as f:
            var_block = True
            for l in f:
                if 'units' in l:
                    var_block = False

                if var_block == True and l.startswith('variable'):
                    if l.split()[2] == 'equal':
                        exec('self.{} = {}'.format(l.split()[1], l.split()[3]))
                    elif l.split()[2] == 'string':
                        exec('self.{} = \'{}\''.format(l.split()[1], l.split()[3]))
                    else:
                        raise ValueError('Unknown keyword found in variable line. Except either equal or string.')

            if var_block == True:
                raise ValueError('Units command not found in input script.')

        # Print extracted parameters
        print('Input parameters extracted from file {}'.format(input_file))
        print('\n'.join('{}: {}'.format(k, v) for k, v in self.__dict__.items()))

    def read_variable(self, variable_name):
        """ 
        
        Read the variable block at the top of the input script. Note that this
        terminates as soon as the 'units' command is found, so any variables found
        afterward will not be found. Use more general 'read_script' if variables or
        properties afterwards are required. 
        
        """

        with open(self.input_file, 'r') as f:
            var_block = True
            for l in f:
                if 'units' in l:
                        var_block = False

                if 'variable' in l and var_block == True and l.split()[1] == variable_name:
                    try:
                        variable_value = float(l.split()[-1])
                    except:
                        variable_value = l.split()[-1]
                        print('{} returned as string instead of float'.format(variable_name))

                    break

        try:
            return variable_value
        except:
            print('{} not found in script {}'.format(variable_name, self.input_file))
            raise

    def read_script(self, input_file, string_to_search, ignore_variables=False):
        """

        Searches complete input script for a string variable and then returns the
        entire line as a string. It is too difficult to handle every possible
        combination within an input script, so use this as a basis for more
        targetted functions. The search string should be specific enough with the
        correct whitespace characters. The function will only find the first
        instance of the string. An optional ignore variables is provided to ignore
        the line if it is variable command. 

        """

        with open(input_file, 'r') as f:
            for l in f:
                if 'variable' in l and ignore_variables:
                    pass
                elif string_to_search in l:
                    string_value = l
                    break

        try:
            return string_value
        except:
            print(string_to_search + 'not found in script ' + input_file)
            raise

def LAMMPS_Writer(input_file, variable_name, variable_value):
    """

    Write input parameter into the variable block. 

    """

    variable_found = False
    with open(input_file + '.new', 'w') as nf:
        with open(input_file, 'r') as of:
            var_block = True
            for l in of:
                if 'units' in l:
                    var_block = False

                if 'variable' in l and var_block == True and l.split()[1] == variable_name:
                    variable_found = True
                    lstr = l.split()
                    nf.write('{}        {} {} {}\n'.format(lstr[0], lstr[1], lstr[2], variable_value))
                else:
                    nf.write(l)

            if variable_found == False:
                raise ValueError('Variable {} not found in file {}\n'.format(variable_name, input_file))


    mv(input_file + '.new', input_file)



