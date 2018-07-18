#! /usr/bin/env python
import numpy as np

class CouetteAnalytical:
    
    def __init__(self, Re, U, Lmin, Lmax, nmodes=1000,npoints=20):
        self.nmodes = nmodes    
        self.npoints = npoints
        self.Re = Re*(Lmax - Lmin)
        self.U = U
        self.L = 1. #Lmax - Lmin
        self.Lmax = Lmax
        self.Lmin = Lmin
        
    def get_vprofile(self,t):    
        """
            Get velocity profile at time t 
            for start-up Couette flow where 
            the top wall starts moving with 
            velocity U at time t=0.
        """

        y = np.linspace(0., self.L, num=self.npoints)
        yspace = np.linspace(self.Lmin, self.Lmax, num=self.npoints)
        vprofile = np.zeros(self.npoints)
        pi = np.pi
        k = 1./self.Re

        #Using matrix operations should be faster than loop
        n = np.arange(1,self.nmodes+1)
        l = (n*pi/self.L)**2.
        yn = np.sin(np.pi*np.outer(y,n))/self.L
        tn = -(-1.)**n * (2/(n*pi)) * self.U * (1. - np.exp(-l*k*t))
        vprofile = np.sum(np.multiply(yn,tn.T),1)
#        for n in range(1,self.nmodes):
#            l = (n*pi/self.L)**2.
#            vprofile = vprofile + (- (-1.)**n              * 
#                                      (2/(n*pi))            * 
#                                      self.U                *
#                                      (1. - np.exp(-l*k*t)) * 
#                                      np.sin(n*pi*y/self.L))

#        print(_vprofile - vprofile)

        vprofile[-1] = self.U # set top value to BC

        # Normalise
        #yspace = yspace + self.lwall
        #yspace = yspace / yspace[-1] 
        #vprofile = vprofile# / vprofile[-1] 

        return yspace, vprofile
    
    def get_Pprofile(self,t,viscosity):

        yspace = np.linspace(self.Lmin, self.Lmax,num=self.npoints)
        tau = np.zeros(self.npoints)    
        pi = np.pi
        k = 1./self.Re

        # Add zeroth mode
        tau = tau + 0.5*self.U*2./self.L

        # Higher frequency modes 
        for n in range(1,self.nmodes):
            l = (n*pi/self.L)**2.
            tau = tau + (
                          (2./self.L)*((-1.)**(n*self.U))*(np.exp(-l*k*t)) *
                          np.cos(n*pi*yspace/self.L)
                        )

        # Normalise yspace
        yspace = yspace + self.lwall
        yspace = yspace / yspace[-1] 
        tau = viscosity * tau
    
        return yspace,tau 
