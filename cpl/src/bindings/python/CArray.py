import ctypes
from abc import ABCMeta, abstractmethod

class CArray:
	__metaclass__ = ABCMeta

	def __init__(self):
		self.dims = []
		self.type = None

	
	@abstractmethod
	def pack(obj, data_type, dims): pass

	@abstractmethod
	def unpack(obj): pass


