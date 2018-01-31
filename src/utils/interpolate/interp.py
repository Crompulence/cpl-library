c_float_p = ctypes.POINTER(ctypes.c_float)
data = numpy.array([[0.1, 0.1], [0.2, 0.2], [0.3, 0.3]])
data = data.astype(numpy.float64)
data_p = data.ctypes.data_as(c_float_p)
