import numpy as np
from cplpy.cpl import CPL


__all__ = ["Packer", "Unpacker"]

class CPLPacker(object):
    def __init__(self, name, region, callback, data_len):
        self.region = region
        self.callback = callback
        self.name = name
        self.data_len = data_len
        self.buffer = None


    def __call__(self, socket):
        self.callback(socket)

    def createBuffs(self, socket):
        ncxl = socket.proc_region.ncxl
        ncyl = socket.proc_region.ncyl
        nczl = socket.proc_region.nczl
	#if socket._realm == CPL.MD_REALM:
        self.buffer = -1*np.ones((self.data_len, ncxl, ncyl, nczl), order='F', dtype=np.float64)
        #else:
        #    self.buffer = 2*np.ones((self.data_len, ncxl, ncyl, nczl), order='F', dtype=np.float64)

    def prepareBuffs(self, socket):
        pass

class Unpacker(CPLPacker):
    def __init__(self, name, region, callback, data_len):
        CPLPacker.__init__(self, name, region, callback, data_len)

    def prepareBuffs(self, socket):
        socket.recv_buff = self.buffer
        socket.send_buff = np.zeros(3, order='F', dtype=np.float64)

class Packer(CPLPacker):
    def __init__(self, name, region, callback, data_len):
        CPLPacker.__init__(self, name, region, callback, data_len)

    def prepareBuffs(self, socket):
        socket.send_buff = self.buffer
        socket.recv_buff = np.zeros(3, order='F', dtype=np.float64)


