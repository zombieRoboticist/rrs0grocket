import convert as c;
import fluids as f;
import phys as p;
import numpy as np;
import matplotlib.pyplot as plt;
from math import pi;

altitude = np.linspace(0,79000, 10000);
velocity= np.linspace(0,350,35);
temperature = np.zeros_like(altitude);
#temp2=np.zeros_like(altitude);
presure= np.zeros_like(altitude);
vsound = np.zeros_like(altitude);
rho= np.zeros_like(altitude);
gravity = np.zeros_like(altitude);
tRatio = np.zeros([len(velocity),len(altitude)]); 
pRatio = np.zeros([len(velocity),len(altitude)]);
mach = np.zeros([len(velocity),len(altitude)]);
mdotIdeal = np.zeros([len(velocity),len(altitude)]);

for i in range(len(altitude)):
	temperature[i]= f.Temp(altitude[i]);
	# temp2[i]=f.TempR(altitude[i]/c.ft)*c.rankine
	presure[i]=f.press(altitude[i]);
	vsound[i]=f.vSound(altitude[i]);
	rho[i]=p.density(altitude[i]);
	gravity[i]=p.g(altitude[i]);
	for j in range(len(velocity)):
		#print(j)
		tRatio[j][i]= temperature[i]/f.tRatio(velocity[j],altitude[i])
		pRatio[j][i] = presure[i]/f.pRatio(velocity[j],altitude[i])
		mach [j][i]=f.mach(velocity[j],altitude[i])
		mdotIdeal[j][i]=f.mdotIdeal(velocity[j],altitude[i], (9*pi*c.inch*c.inch))

def graphTemp():
	plt.figure(1, figsize=(10, 6));
	plt.plot(altitude, temperature)   
	# plt.plot(altitude,temp2)    
	plt.title( "aproximate temperature at altitudes" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("Temp (K)", fontsize=12)
	plt.show();

def graphPress():
	plt.figure(1, figsize=(10, 6));
	plt.plot(altitude, presure)       
	plt.title( "aproximate Presure at altitudes" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("Presure (Pa)", fontsize=12)
	plt.show();

def graphVsound():
	plt.figure(1, figsize=(10, 6));
	plt.plot(altitude, vsound)       
	plt.title( "aproximate velocity of sound at altitudes" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("V sound (m/s)", fontsize=12)
	plt.show();

def graphRho():
	plt.figure(1, figsize=(10, 6));
	plt.plot(altitude, rho)       
	plt.title( "aproximate densities at altitudes" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("Density (kg/m^3)", fontsize=12)
	plt.show();

def graphG():
	plt.figure(1, figsize=(10, 6));
	plt.plot(altitude, gravity)       
	plt.title( "aproximate graviational acceleration at altitudes" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("g (m/s^2)", fontsize=12)
	plt.show();

def graphTt():
	plt.figure(1, figsize=(10, 6));
	for i in range(len(velocity)):
		name = "V = "+ str(velocity[i])+" m/s"
		plt.plot(altitude, tRatio[i], label = name);       
	plt.title( "aproximate temperatures at altitudes and speeds" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("Temp (K)", fontsize=12)
	plt.legend();
	plt.show();

def graphPt():
	plt.figure(1, figsize=(10, 6));
	for i in range(len(velocity)):
		name = "V = "+ str(velocity[i])+" m/s"
		plt.plot(altitude, pRatio[i], label = name);       
	plt.title( "aproximate presures at altitudes and speeds" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("Presure (Pa)", fontsize=12)
	plt.legend();
	plt.show();

def graphMach():
	plt.figure(1, figsize=(10, 6));
	for i in range(len(velocity)):
		name = "V = "+ str(velocity[i])+" m/s"
		plt.plot(altitude, mach[i], label = name);       
	plt.title( "aproximate Mach numbers at altitudes and speeds" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("Mach", fontsize=12)
	plt.legend();
	plt.show();

def graphMdotIdeal():
	plt.figure(1, figsize=(10, 6));
	for i in range(len(velocity)):
		name = "V = "+ str(velocity[i])+" m/s"
		plt.plot(altitude, mdotIdeal[i], label = name);       
	plt.title( "aproximate mass flows at altitudes and speeds" , fontsize=14)
	plt.xlabel(" altitude (m)"        , fontsize=12)
	plt.ylabel("Mass Flow (kg/s)", fontsize=12)
	plt.legend();
	plt.show();

graphG();
graphRho();
graphTemp();
graphPress();
graphVsound();
graphTt();
graphPt();
graphMach();
graphMdotIdeal();
