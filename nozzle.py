import convert as c;
import fluids as f;
import phys as p;
import numpy as np;
import propeller as pr;
from math import pi;
# import rocket;
from scipy.optimize import fsolve,root


class Nozzle:

	cp= 1.007 #kJ/kg/K for -136 to 152 F
	gama=1.4
	ain=0;
	ath=0;
	aex=0;
	mdot=0;
	vinf=0;
	L=0;
	R=287;

	def __init__(self):
		""" """
		return None;
	def setAin(self, ain):
		self.ain=ain;
		return None

	def setAth (self, Ath):
		""" """ 
		self.ath=Ath;
		return None;

	def setAex(self,Aex):
		""" """
		self.aex=Aex;
		return None;

	def setA (self, Ain, Ath, Aex):
		""" """
		self.setAex(Aex);
		self.setAth(Ath);
		self.setAin(Ain);
		return None
	def setL(self,l):
		self.L=l;
		return None

	def getmdot(self,vinf,z):

		rho= p.correctedDensity(vinf,z)
		return rho*self.ain*vinf

	def setmdot(self,m):
		self.mdot=m;
		return None;

	def isThrottled(self, vinf, z):

		arat=f.Aratio(vinf,z)
		if(self.ain/self.ath)>=arat:
			return True
		return False

	def getVout(self, vinf,z, fd):
		""" depricated... not sure if valid"""
		if self.isThrottled(vinf,z):
			print("throtled")

		else:
			print('not')
			t1=p.correctedTemperature(vinf,z)
			w= -fd*vinf
			return fsolve( lambda x: w/self.getmdot(vinf,z)-(vinf**2)/2-self.cp*t1+(x**2)/2+self.cp*p.correctedTemperature(x,z),vinf*.9 )[0]

	def vOutnvst(self, pin, pout):
		re2 = self.aex/pi;
		return re2*((pin)-(pout))/16/.0000171/self.L/.5

	def vouttank(self,pt,tt,mdot):
		t=lambda x: tt*(x/pt)**((self.gama-1)/self.gama)
		m2=lambda v,t: (v**2)/(self.gama*self.R*t)
		p2 = lambda v,t: pt/((1+(self.gama-1)*m2(v,t)/2)**(self.gama/(self.gama-1)))
		v1= mdot/self.ain/pt*self.R*tt
		# return root((lambda x: (x**2)-((mdot*287)**2)*(((tt**(7))/(tt-(x**2)/2/self.cp)**(5)))/(self.aex/pt)**2), (mdot/1.28/self.aex*10)**2)
		# return root((lambda x: [pt-x[0]*(tt/x[1])**(self.gama/(self.gama-1)), tt-x[1]-(x[2]**2)/2/self.cp, mdot-x[0]*self.aex*x[2]/287/x[1]]), [pt*10,tt*10,mdot/self.aex/1.28*10])
		# return root((lambda x:mdot-x/287/t(x)*self.aex*(2*self.cp*(tt-t(x)))**.5), pt*.9)
		# return root((lambda x: mdot/self.aex-pt*(2/287/tt*self.gama/(self.gama-1)*(((x/pt)**(2/self.gama))-((x/pt)**((self.gama+1)/self.gama))))**.5),pt*.3)
		# return fsolve((lambda x: [mdot-p2(x[0],x[1])*x[0]*self.aex/self.R/x[1], tt-x[1]*(1+(self.gama-1)*m2(x[0],x[1])/2)]), [mdot/1.28/self.aex,tt])
		# return root((lambda x: [mdot-p2(x[0],x[1])*x[0]*self.aex/self.R/x[1], p2(x[0],x[1])/pt + self.ain/self.aex+ tt*(x[0]**2)/x[1]/(v1**2)-x[0]/v1 ]), [mdot/1.28/self.aex,tt])

	# def rhoOutnvst(self, mdot,z,rtube):
	# 	re2 = self.aex/pi;
	# 	return mdot/(3*(f.vSound(z)**2)/((rtube**4)-(re2**2))**.5)/pi/re2**2
	def stateCon(self, vo, po,mo,rhoo,to, cf=, n=10000;):
		v= np.array([vo]);
		p=np.array([po]);
		m=np.array([mo]);
		rho=np.array([rhoo]);
		T=np.array([to]);
		dx=self.L/n;
		rin= (self.ain/pi)**.5;
		rth=(self.ath/pi)**.5
		dr=(rth-rin)/n
		da=pi*dr**2
		for i in range(n):
			v.append(v[-1]+v[-1]*(-1/(1-m[-1]**2)*da/(self.ain+i*da)+(1.4*(m[-1]**2)*(1+.4/2*m[-1]**2)/(1-m[-1]**2))*4*cf*dx/(rth+i*dr)/2));
			p.append(p[-1]+p[-1]*(1.4*(m[-1]**2)/(1-m[-1]**2)*da/(self.ain+i*da)-1.4(m[-1]**2)*(1+.4*m[-1]**2)/(1-m[-1]**2)/2*4*cf*dx/(rth+i*dr)/2));
			rho.append(rho[-1]+rho[-1]*((m[-1]**2)/(1-m[-1]**2)*da/(self.ain+i*da)-1.4(m[-1]**2)/2/(1-m[-1]**2)*4*cf*dx/(rth+i*dr)/2));
			T.append(T[-1]+T[-1]*(.4*(m[-1]**2)/(1-m[-1]**2)*da/(self.ain+i*da)-1.4*.4*(m[-1]**4)/2/(1-m[-1]**2)*4*cf*dx/(rth+i*dr)/2));
			m.append(m[-1]+m[-1]*(-2*(1+.2*m[-1]**2)/(1-m[-1]**2)*da/(self.ain+i*da)+1.4*(m[-1]**2)*(1+.2*m[-1]**2)/(1-m**2)*4*cf*dx/(rth+i*dr)/2)**.5);
		return v,p,m,rho,T;


