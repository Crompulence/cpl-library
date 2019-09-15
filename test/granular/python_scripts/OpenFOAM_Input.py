import re
import PyFoam.Basics.DataStructures
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile as ppf
from shutil import move as mv

class OpenFOAM_Input(object):
    """

    Class to read and write into a OpenFOAM case structure. Relies on PyFoam
    which can be installed by 'pip install PyFoam'.

    """

    def __init__(self, case_dir='./openfoam'):
        self.case = case_dir

    def read_0_file(self, file_name, boundaryFace):
        """

        Read input data from any of the variables in the 0 directory for a
        given boundary face. The available types of the boundary face are
        currently zeroGradient, calculated and fixedValue. For zeroGradient
        and calculated, there is a no value. For fixedValue, either a scalar
        or vector (as a list) is returned. Only uniform values are handled.

        """

        f = ppf(file_name)
        typeFace = f['boundaryField'][boundaryFace]['type']
        if typeFace == 'zeroGradient' or typeFace == 'calculated':
            valueFace = None
        elif typeFace == 'fixedValue':
            try:
                valueFace = f['boundaryField'][boundaryFace]['value']['uniform']
                if isinstance(valueFace, PyFoam.Basics.DataStructures.Vector):
                    valueFace = [valueFace[0], valueFace[1], valueFace[2]]
            except:
                print('Attempted to read a uniform boundary face value of:')
                print(str(f['boundaryField'][boundaryFace]['value']))
                print('Non-uniform input currently not handled.')
                raise
        else:
            print('Unknown type {} found for boundary face {} in file {}'.format(
                typeFace, boundaryFace, file_name))

        return typeFace, valueFace


    def read_constant_file(self, file_name, property_name):
        """

        Read input data from any of the files in the constant directory,
        including transportProperties and environmentalProperties. For the
        property of interest, either the scalar (e.g. fluid density) or vector
        (e.g. gravity) values are returned.

        """

        f = ppf(file_name)
        prop = f[property_name][-1]
        if isinstance(prop, PyFoam.Basics.DataStructures.Vector):
            prop = [prop[0], prop[1], prop[2]]

        return prop

    def read_system_file(self, file_name, property_name):
        """

        Read controlDict is the system directory to extract any property. 

        """

        f = ppf(file_name)
        prop = f[property_name]

        return prop

def OpenFOAM_Writer_0_File(file_name, boundaryFace, valueFace, isScalar=True, axis_val=0):
    """

    Write input parameter into any of the variables in the 0 directory.
    The value of a boundary face must be provided, and where the parameter
    is a vector quantity, the element of the vector should be specified.
    This is currently limited to boundary faces with type fixedValue.

    """

    f = ppf(file_name)
    if isScalar:
        f['boundaryField'][boundaryFace]['value']['uniform'] = valueFace
    else:
        f['boundaryField'][boundaryFace]['value']['uniform'][axis_val] = valueFace

    f.writeFile()

    # Format the written file to similar format as previous
    with open(file_name + '.new', 'w') as nf:
        nf.write('{}\n'.format(r'/*--------------------------------*- C++ -*----------------------------------*\ '))
        nf.write('{}\n'.format(r'| =========                 |                                                 | '))
        nf.write('{}\n'.format(r'| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | '))
        nf.write('{}\n'.format(r'|  \\    /   O peration     | Version:  2.2.2                                 | '))
        nf.write('{}\n'.format(r'|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | '))
        nf.write('{}\n'.format(r'|    \\/     M anipulation  |                                                 | '))
        nf.write('{}\n'.format(r'\*---------------------------------------------------------------------------*/ '))

        with open(file_name, 'r') as of:
            Foamfile_block = False
            boundaryField_block = False
            for l in of:
                if 'FoamFile' in l:
                    Foamfile_block = True

                if Foamfile_block:
                    if l.startswith(' '):
                        nf.write('    ' + l.lstrip())
                    else:
                        nf.write(l)

                    if l.startswith('}'):
                        nf.write('{}\n'.format(r'// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'))
                        nf.write('\n')
                        Foamfile_block = False

                if 'dimensions' in l:
                    nf.write(l)
                    nf.write('\n')

                if 'internalField' in l:
                    nf.write(l)
                    nf.write('\n')

                if 'boundaryField' in l:
                    boundaryField_block = True

                if boundaryField_block:
                    if l.startswith('    '):
                        nf.write('    ' + '    ' + l.lstrip())
                    elif l.startswith('  '):
                        nf.write('    ' + l.lstrip())
                    elif l.startswith('}'):
                        nf.write('}\n')
                        nf.write('{}\n'.format(r'// ************************************************************************* //'))
                        boundaryField_block = False
                    else:
                        nf.write(l)

    mv(file_name + '.new', file_name)