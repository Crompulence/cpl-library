class LAMMPS_Input(object):
    """

    Class to read and write into LAMMPS input script.
    
    """

    def __init__(self, case_dir='./lammps'):
        self.case = case_dir


    def read_variables(self, input_file, variable_name):
        """ 
        
        Read the variable block at the top of the input script. Note that this
        terminates as soon as the 'units' command is found, so any variables found
        afterward will not be found. Use more general 'read_script' if variables or
        properties afterwards are required. 
        
        """

        with open(input_file, 'r') as f:
            var_block = True
            for l in f:
                if 'units' in l:
                        var_block = False

                if 'variable' in l and var_block == True and l.split()[1] == variable_name:
                    try:
                        variable_value = float(l.split()[-1])
                    except:
                        variable_value = l.split()[-1]
                        print(variable_name + ' returned as string instead of float.')

                    break

        try:
            return variable_value
        except:
            print(variable_name + 'not found in script ' + input_file)
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


