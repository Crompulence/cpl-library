import cpl
import numpy as np

lib = cpl.CPL()

int_p = 5
doub_p = 4.0
bool_p = False
int_ptr = np.array([[9,2,3],[4,5,6],[7,8,9], [1,3,4]], order='F',dtype=np.int32)
print int_ptr.dtype
doub_ptr = np.array([[2.,5.,4.], [2.,3.,4.], [1.,2.,3.],[4.,5.,6.],[7.,8.,9.]], order='F')
print doub_ptr.dtype
#doub_ptr = np.asanyarray([1.,2.,3.,4.,5.,6.,7.,8.,9.], order='F')
lib.test_python(int_p, doub_p, bool_p, int_ptr, doub_ptr)
