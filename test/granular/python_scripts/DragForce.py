import numpy as np 

class DragForce(object):
	""" 

	Base class for various drag force models implemented in coupled simulations
	with attributes drag coefficients (B), drag force (F) on an individual
	particle and pressure gradient (dP) across a sample.

	"""

	def __init__(self, B=None):
		self.B = B
		self.MIN_REL_VELOCITY = 1.e-6 

	def calculate_drag_force(self, Uf, Vp):
		"""

		Drag force is equal to the drag coefficient multiplied by the relative
		velocity between the particle velocity (Vp) and fluid velocity (Uf).

		"""

		if self.B == None: 
			raise('The drag coefficient, B, must be set to calculate the drag force.')

		F = self.B*(Uf - Vp)
		return F

	def calculate_pressure_gradient(self, Uf, Vp):
		"""

		The pressure gradient is calculated as a function of the porosity or void
		fraction (epsf), the mean particle volume which is a function of the mean
		particle diameter (dp), and the calculated drag force (F).

		"""

		F = self.B*(Uf - Vp)
		volume = (np.pi/6.)*(self.dp**3)
		dP = (1. - self.epsf)*(F/volume)
		return dP

class Stokes(DragForce):
	""" 

	Derived class for Stokes drag Law.

	"""

	def __init__(self, muf, epsf, dp):
		super(Stokes, self).__init__()
		self.muf = muf
		self.dp = dp
		self.epsf = epsf
		self.B = 3*np.pi*muf*dp*epsf

class DiFelice(DragForce):
	"""

	Derived class for Di Felice drag force model

	"""

	def __init__(self, muf, rhof, epsf, dp, Uf, Vp):
		super(DiFelice, self).__init__()
		self.muf = muf
		self.rhof = rhof
		self.dp = dp
		self.epsf = epsf

		if abs(Vp - Uf) < self.MIN_REL_VELOCITY:
			self.B = 0.
		else:
			Re = rhof*dp*epsf*abs(Uf - Vp)/muf
			chi = 3.7 - 0.65*np.exp(-0.5*((1.5 - np.log10(Re))**2))
			Cd = (0.63 + 4.8/np.sqrt(Re))**2
			self.B = 0.5*Cd*rhof*(np.pi/4)*(dp**2)*(epsf**(2-chi))*abs(Uf - Vp)

	def update_drag_coefficient(self, Uf, Vp):
		self.__init__(self.muf, self.rhof, self.epsf, self.dp, Uf, Vp)

class Ergun(DragForce):
	""" 

	Derived class for Ergun drag force model

	"""

	def __init__(self, muf, rhof, epsf, dp, Uf, Vp):
		super(Ergun, self).__init__()
		self.muf = muf
		self.rhof = rhof
		self.dp = dp
		self.epsf = epsf

		if abs(Vp - Uf) < self.MIN_REL_VELOCITY:
			self.B = 0.
		else:
			self.B = (150.*np.pi*muf*dp/6)*((1-epsf)/epsf) + (1.758*np.pi*rhof*(dp**2)/6)*abs(Vp - Uf)

	def update_drag_coefficient(self, Uf, Vp):
		self.__init__(self.muf, self.rhof, self.epsf, self.dp, Uf, Vp)