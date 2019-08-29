import numpy as np 

class DragForce(object):
	""" 

	Base class for various drag force models implemented in coupled simulations
	with attributes drag coefficients (B), drag force (F) on an individual
	particle and pressure gradient (dP) across a sample.

	"""

	def __init__(self, B=None):
		self.B = B 

	def calculate_drag_force(self, Uf, Vp):
		"""

		Drag force is equal to the drag coefficient multiplied by the relative
		velocity between the particle velocity (Vp) and fluid velocity (Uf).

		"""

		if self.B == None: 
			raise('The drag coefficient, B, must be set to calculate the drag force.')

		self.F = self.B*(Uf - Vp)
		return self.F

	def calculate_pressure_gradient(self, epsf, Uf, Vp, dp):
		"""

		The pressure gradient is calculated as a function of the porosity or void
		fraction (epsf), the mean particle volume which is a function of the mean
		particle diameter (dp), and the calculated drag force (F).

		"""

		self.F = self.B*(Uf - Vp)
		self.dP = ((1. - epsf)/epsf)*(self.F/(np.pi*dp**3/6))
		return self.dP

class Stokes(DragForce):
	""" 

	Derived class for Stokes drag Law.

	"""

	def __init__(self, muf=1.e-3, dp=0.001):
		super(Stokes, self).__init__()
		self.B = 3*np.pi*muf*dp

class DiFelice(DragForce):
	"""

	Derived class for Di Felice drag force model

	"""

	def __init__(self, muf=1.e-3, rhof=1000., Uf=0.01, dp=0.001, Vp=0., epsf=0.36):
		super(DiFelice, self).__init__()
		if abs(Vp - Uf) == 0.:
			self.B = 0.
		else:
			Re = rhof*dp*abs(Uf - Vp)/muf
			chi = -(3.7 - 0.65*np.exp(-0.5*((1.5 - np.log10(Re))**2)))
			Cd = (0.63 + 4.8/np.sqrt(Re))**2
			self.B = 0.5*Cd*rhof*(np.pi/4)*(dp**2)*(epsf**chi)*abs(Uf - Vp)