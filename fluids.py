import numpy as np;
import convert as c;
from math import sqrt;
import phys as p;

gama = 1.4 ; #ratio of specific heats for air

def ReX (rho, v, mu, x):
	""" calculates the reynalds number at a point

	inputs: rho (fluid density kg/m/m/m), v (free streem fluid velocity m/s), mu (fluid viscosity kg/m/s), x (length of member at point m) 
	output: reynalds number at that point"""
	return rho*v*x/mu;

def ReList(rho, v, mu, x):
	""" calculates the reynalds number at a list of points

	inputs: rho (fluid density kg/m/m/m), v (free streem fluid velocity m/s), mu (fluid viscosity kg/m/s), x (list of length of member at point m) 
	output: list reynalds number at the input points"""
	re = np.zeros_like(x);
	for i in range(len(x)):
		re[i]= ReX(rho, v, mu, x[i]);
	return re;

def TempR(alt):
	""" estimates temperature at different altitudes up to 262448

	Input: alt (altitude in ft)
	output: absolute temperature in R"""
	T = 518.6;
	if (alt <= 36152.):  #          // Troposphere
		T = 518.6 - 3.56 * alt/1000. ;
	if (alt >= 36152. and alt <= 82345.):  #  // Stratosphere
		T = 389.98 ;
	if (alt >= 82345. and alt <= 155348.):           
		T = 389.98 + 1.645 * (alt-82345.)/1000. ;
	if (alt >= 155348. and alt <= 175346.):           
		T = 508.788 ;
	if (alt >= 175346. and alt <= 262448.):           
		T = 508.788 - 2.46888 * (alt-175346.)/1000. ;
	if alt>262448.:
		T=0;
	return T; #rankine

def Temp(z):
	"""  estimates temperature at different altitudes up to 80,000 m

	Input: z (altitude in m)
	output: absolute temperature in K"""

	alt = z/c.ft
	T = TempR(alt)
	return T*c.rankine #k

def vSound(z):
	"""calculates the velocity of sound at a given altitude
	source: NASA https://www.grc.nasa.gov/WWW/k-12/rocket/sound.html
	Input: z (altitude in m)
	output: Speed of sound in m/s"""
	rgas = 1718. ;        #        /* ft2/sec2 R */
	alt = z/c.ft
	T = TempR(alt) #rankine

	a0 = sqrt(gama*rgas*T) ;  # feet /sec
	return a0 *c.ft ;   # m/s

def mach(v, z):
	"""calculates the mach number of a object traveling though earth's atmosphere

	input: v (velocity of object in m/s) and z (altitude of the object in m)
	output: mach number"""
	return v/vSound(z);

def press (z):
	""" estimates atmospheric presures at different altitudes for ideal gas

	input: z (altitude in m)
	output: pressure in Pa"""
	R = 286. #J/kg/K
	return R*Temp(z)*p.density(z) #Pascals

#for isentropic gasses only 

def tRatio (v,z):
	"""calculates the compressibility ratio for temperature

	input: v (the velocity of the fluid m/s) and z (the altitude of the fluid m)
	output: the temperature ratio T/Tt"""
	M = mach(v,z)
	return 1/(1+.5*(gama-1)*M*M)

def pRatio(v,z):
	"""calculates the compressibility ratio for presure

	input: v (the velocity of the fluid m/s) and z (the altitude of the fluid m)
	output: the pressure ratio P/Pt"""
	return tRatio(v,z)**(gama/(gama-1))

def rhoRatio(v,z):
	""" calculates the compressibility ratio for density

	input: v (the velocity of the fluid m/s) and z (the altitude of the fluid m)
	output: the density ratio rho/rhot"""
	return tRatio(v,z)**(-1/(gama-1));

def q(v,z):
	""" computes the dynamic pressure of incompressible air

	input: v (the velocity of the fluid m/s) and z (the altitude of the fluid m)
	output: the dynamic pressure"""
	return p.density(z)*v*v/2

def qComp(v,z):
	""" computes the dynamic pressure of compressible air (?) using density comperesibility ratio
*check with woodcock*
	input: v (the velocity of the fluid m/s) and z (the altitude of the fluid m)
	output: the dynamic pressure"""
	return .5*p.density(z)/rhoRatio(v,z)*v*v;

def qRatio(v,z):
	""" computes the dynamic pressure ratio of air based on mach number (not sure if is commpressible)

	input: v (the velocity of the fluid m/s) and z (the altitude of the fluid m)
	output: the dynamic pressure divided by the gas pressure (q/p)"""
	M=mach(v,z);
	return gama*M*M/2

def mdotIdeal(v,z,A):
	""" calculates ideal gas mass flow rate of atmosphere

	input: v (velocity m/s), z (altitude m), and A (area m^2)
	output: ideal gas mass flow rate kg/s"""
	R = 1716 *c.ft*c.ft/c.rankine
	return (A*(press(z)/pRatio(v,z))/sqrt(Temp(z)/tRatio(v,z)))*sqrt(gama/R)*mach(v,z)*((tRatio(v,z))**(-1*(gama+1)/(2*(gama-1))));