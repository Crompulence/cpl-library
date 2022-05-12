import string
import os
import numpy as np

class InputMod:

    def __init__(self,filename):
        self.filename = filename

    def read_input(self, keyword):    

        """ 
            Read all values underneath the first appearance of
            keyword until blank line
        """
        readvals = False
        vals = []
        with open(self.filename) as f:
            for line in f.readlines():
                v = line.strip("\n")
                if readvals:
                    try:
                        vals.append(int(v))
                    except ValueError:
                        try:
                            vals.append(float(v))
                        except ValueError:
                            readvals = False
                            return vals
                elif (keyword == v): 
                    readvals = True

            
            print(('Input string ' + keyword +' not found'))

    
    def replace_input(self, keyword, keyvals):    

        """ 
            Replace N values underneath the first appearance of
            keyword with keyvals (length N)

        """

        found = False
        key_lno = 0 # Keyword linenumber

        for line in open(self.filename):
        
            key_lno += 1

            # Take into account keyword might be "turned off" by a
            # comment character before it
            if (line[0:len(keyword)]   == keyword or 
                line[1:len(keyword)+1] == keyword ): 

                # Mark the keyword as found 
                found = True

                # Ensure keyword is activated (i.e. not commented out)
                sedstr = ( "sed -i '" + str(key_lno) + "s/.*/" + keyword + 
                           "/' " + self.filename ) 
                os.system(sedstr)

                # Values start on next line
                val_lno = key_lno + 1

                if type(keyvals) is list:

                    for val in keyvals:

                        if (val != None):

                            sedstr = ( "sed -i '" + str(val_lno) + "s/.*/" + 
                                        str(val) + "/' " + self.filename ) 
                            os.system(sedstr)
                    
                        val_lno += 1

                else:

                    sedstr = ( "sed -i '" + str(val_lno) + "s/.*/" + 
                                str(keyvals) + "/' " + self.filename ) 
                    os.system(sedstr)

                # Stop looping through the file
                break
        
        if ( found == False ):

            with open(self.filename,'a') as f:
                f.write(keyword+'\n')
                for keyval in keyvals:
                    f.write(keyval+'\n')
    
            quit('Input string ' + keyword + 
                 ' not found, appended to file instead.')


#List of Dictonary classes with added routines to add and multiple inputs
class InputList(list):
    def __init__(self,*arg,**kw):
        super(InputList, self).__init__(*arg, **kw)

        self.valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)

    #Define Addition operation as elementwise addition
    def __add__(self, x):

        if (type(x) == InputDict):
            Listx = x.expand()
        elif (type(x) == InputList):
            Listx = x
        else:
            raise TypeError("Unsupported type " + str(type(x))  
                            + " for input addition" + 
                            " -- must be InputList or InputDict type")

        returnlist = InputList()
        ziplist = list(zip(self,Listx))
        for entry in ziplist:
            tempdict = {}
            for dic in entry:
                tempdict.update(dic)
            returnlist.append(tempdict)

        return returnlist

    __radd__=__add__

    #Define multiplication operation as all permutations of inputs
    def __mul__(self,x):

        if (type(x) == InputDict):
            Listx = x.expand()  
        elif (type(x) == InputList):
            Listx = x
        else:
            raise TypeError("Unsupported type " + str(type(x))  
                            + " for input multiplication" + 
                            " -- must be InputList or InputDict type")

        returnlist = InputList()
        for entry1 in self:
            for entry2 in Listx:
                newdict = {}
                newdict.update(entry1)
                newdict.update(entry2)
                returnlist.append(newdict)

        return returnlist

    __rmul__=__mul__

    #Define Addition operation as elementwise addition
    def zip_inputs(self, x):

        returnlist = self + x

        return returnlist

    #Generate all permutations of inputs and return filenames
    def outer_product_inputs(self,x):

        returnlist = self*x

        return returnlist

    def filenames(self):

        #Generate list containing filenames
        filenames = []
        for name in self:
            filename = ''
            for key, value in list(name.items()):
                if type(value) is float:
                    value = np.round(value,2)
                if type(value) is int:
                    value = value

                # Combine key and value with invalid characters removed
                kvstr = str((key,value))
                kvstr = kvstr.replace('.','p')
                filename = filename + (''.join(c for c in kvstr 
                                       if c in self.valid_chars[6:]))

            filenames.append(filename)

            #print(name,filename)


        #First remove any invalid characters
        #filenames=[(''.join(c for c in str(name.items()) 
        #                  if c in self.valid_chars[6:]))
        #                  for name in self]
         
        return filenames

#Dictonary class with added routines to add and multiple inputs
class InputDict(dict):
    def __init__(self,*arg,**kw):
        super(InputDict, self).__init__(*arg, **kw)

    #Expand InputDict with multiple values per entry into InputList
    def expand(self):

        expansion = list(self.values())[0]
        returnlist = InputList({list(self.keys())[0]:e} for e in expansion)
        
        return returnlist       

    #Wrapper to convert InputDict to InputList then add
    def __add__(self, x):

        # Convert to lists
        templist = self.expand()
        returnlist = templist + x

        return returnlist

    __radd__=__add__

    #Wrapper to convert InputDict to InputList then multiply
    def __mul__(self,x):

        # Convert to lists
        templist = self.expand()
        returnlist = templist * x

        return returnlist

    __rmul__=__mul__
