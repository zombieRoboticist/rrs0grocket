import convert as c;
import fluids as f;
import phys as p;
import numpy as np;
import propeller as pr;
from math import pi;
import nozzle as nz;
from scipy.optimize import fsolve,root

class Rocket:

	drag = np.zeros((1,1));
	prop = pr.Propeller();
	diamater=0;
	noz = nz.Nozzle();
	# Ain=0;
	# Aout=0;


	def __init__(self):


		return None;

	def setDiameter(self, d):
		"""sets diameter of the rocket

		inputs: self, d (diameter m) 
		output: none"""

		self.diameter=d;
		return None;

	def setCd (self, cd):
		""" """
		self.drag = np.array(cd);
		return None;

	def getDrag(self, v,z):
		""" finds the drag on the rocket at a given speed

		Inputs:	self, z (altitude m), v (velocity m/s)
		Output: drag force (N)	"""

		i=0;
		m=f.mach(v,z);
		# print(m);
		for j in range (self.drag.shape[1]-1):
			if m >= self.drag[0][j] and m < self.drag[0][j+1]:
				i = j;

			if m>self.drag[0][-1] or m<self.drag[0][0]:
				raise ValueError("mach number out of range of cd values")
				return None

	    #finds dc using linear interpolation
		r= (self.drag[0][i] -m)/(self.drag[0][i] - self.drag[0][i+1]);

		cd= self.drag[1][i+1]*r +self.drag[1][i]*(1-r);
		

		return (self.diameter**2)*pi* cd*f.qComp(v,z)/4;



	def setProp(self):
		""" """

	def setNozzle (self, nozzle):
		""" """ 
		self.noz=nozzle
		# self.Ain = nozzle.ain
		# self.Aout=nozzle.aex

	# def setAin(self,A):
	# 	""" """
	# 	self.Ain=A;
	# 	return None;

	# def setAex(self,A):
	# 	""" """
	# 	self.Aout=A;
	# 	return None;

	def getMinVout(self,v,z, thrust):
		""" """
		mdot=self.noz.ain*p.correctedDensity(v,z)*v;

		return fsolve(( lambda x:thrust+(p.correctedPresure(v,z)-p.correctedPresure(x,z))*self.noz.aex/mdot+v-x),v*.9)[0];

	# thrust/ax/rhoout +ain*vin^2/ax
	def getPower(self,vin,vout,z):
		""" """
		mdot=self.noz.ain*p.correctedDensity(vin,z)*vin;
		cp = 1007 #j/kg/k
		return mdot*(vout**2-vin**2)/2+mdot*cp*(p.correctedTemperature(vout,z)-p.correctedTemperature(vin,z))
