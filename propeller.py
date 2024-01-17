import convert as c;
import fluids as f;
import phys as p;
import numpy as np;
from math import pi;

class Propeller:

	a=0; #velocity increment coeficient at propeller (V[1+a], p'+delta p)
	b= 0; #velocity increment coeficient at (V[1+b],p)
	r=0; #radius of propeller
	n=0;#num of blades
	width = 0; #width of the blade at r=.75R
	TcList=np.zeros((1,1));
	QcList=np.zeros((1,1));
	VND= lambda self, M, z, w:  M*f.vSound(z)/w/self.r/2;  #computes V/nD, inputs: self mach number, altitude, angular velocity
	S= lambda self: 2*pi*self.r/self.n/self.width;# somputse solidity factor, inputs: self


	def __init__ (self):
		""" initializer """

	def setr(self,r):
		""" """
		self.r=r;
		return None;

	def setA(self,a):
		""" sets the values of a and b

		input: self, a (velocity increment coeficient at propeller (V[1+a], p'+delta p))
		output: none"""
		self.a=a
		self.b=2*b
		return None

	def setB(self,b):
		""" sets the values of a and b

		input: self, (b velocity increment coeficient at (V[1+b],p))
		output: none"""
		self.b=b
		self.a=b/2
		return None

	def ADisk(self):
		""" returns the area of the propeller as a disk

		input: self
		output: area or disk"""
		return pi*self.r**2

	def efficiencyIdeal(self):
		""" returns the efficiency of the propeller

		input: self
		output: ideal efficiency of the propeller"""
		return 1/(1+self.a)

	def thrustIdeal(self,v,z):
		""" computes the energy produced by the thrust of the propeller

		input: self, v (velocity of fluid m/s), z (altitude of propeller m)
		output: thrust (mass flow rate*velocity added to fluid)"""
		return .5*f.mdotIdeal(v,z,self.ADisk())*(self.b*v)**2

	def setTc( self,c ):
		"""sets the thrust coefficients

		inputs: self, c (a 2d array organized as [ [velocity list], [coefficient list] ])
		output: null """
		self.Tc = np.array(c);
		return None

	def thrustFromTc (self, M, z, w):
		"""computes the propeller thrust from the thrust coefficient

		Input: self, M (mach number), z (altitude (m)), w (angular velocity (rad/s))
		output: thrust (N)"""

		v= M*f.vSound(z); 
		vnd = self.VND(M,z,w);

		i=0;
		for j in range (self.Tc.shape(1)-1):
			if vnd >= self.Tc[0][j] and vnd < self.Tc[0][j+1]:
				i = j;

		if vnd>self.Tc[0][-1] or vnd<self.Tc[0][0]:
			raise ValueError("velocity/angular velocity/diameter out of range of Tc Values")
			return None

		#finds Tc using linear interpolation
		r= (self.Tc[0][i]-vnd)/(self.Tc[0][i] - self.Tc[0][i+1]);

		Ct= self.Tc[1][i+1]*r +self.Tc[1][i]*(1-r);

		return Ct*((self.r*2)**2)*f.qComp(v,z)*2;

	def setQc( self,c ):
		"""sets the torque coefficients

		inputs: self, c (a 2d array organized as [ [velocity list], [coefficient list] ])
		output: null """
		self.Tc = np.array(c);
		return None

	def torqueFromQc (self, M, z, w):
		"""computes the propeller torque from the torque coefficient

		Input: self, M (mach number), z (altitude (m)), w (angular velocity (rad/s))
		output: torque (N*m)"""

		v= M*f.vSound(z); 
		
		i=0;
		for j in range (self.Qc.shape(1)-1):
			if v >= self.Qc[0][j] and v < self.Qc[0][j+1]:
				i = j;

		if v>self.Qc[0][-1] or v<self.Qc[0][0]:
			raise ValueError("velocity/angular velocity/diameter out of range of Qc Values")
			return None

		#finds Qc using linear interpolation
		r= (self.Qc[0][i]-v)/(self.Qc[0][i] - self.Qc[0][i+1]);

		Cq= self.Qc[1][i+1]*r +self.Qc[1][i]*(1-r);

		return Cq*((self.r*2)**3)*f.qComp(v,z)*2;

	def efficiencyFromCtCp( self, M, z, w ):
		"""computes the propeller torque from the torque coefficient

		Input: self, M (mach number), z (altitude (m)), w (angular velocity (rad/s))
		output: efficiency"""

		v= M*f.vSound(z); 
		vnd = self.VND(M,z,w);

		i=0;
		for j in range (self.Tc.shape(1)-1):
			if vnd >= self.Tc[0][j] and vnd < self.Tc[0][j+1]:
				i = j;

		if vnd>self.Tc[0][-1] or vnd<self.Tc[0][0]:
			raise ValueError("velocity/angular velocity/diameter out of range of Tc Values")
			return None

		#finds Tc using linear interpolation
		r= (self.Tc[0][i]-vnd)/(self.Tc[0][i] - self.Tc[0][i+1]);

		Ct= self.Tc[1][i+1]*r +self.Tc[1][i]*(1-r);

		i=0;
		for j in range (self.Qc.shape(1)-1):
			if v >= self.Qc[0][j] and v < self.Qc[0][j+1]:
				i = j;

		if v>self.Qc[0][-1] or v<self.Qc[0][0]:
			raise ValueError("velocity/angular velocity/diameter out of range of Qc Values")
			return None

		#finds Qc using linear interpolation
		r= (self.Qc[0][i]-v)/(self.Qc[0][i] - self.Qc[0][i+1]);

		Cq= self.Qc[1][i+1]*r +self.Qc[1][i]*(1-r);

		return vnd*Ct/Cq/2/pi;
