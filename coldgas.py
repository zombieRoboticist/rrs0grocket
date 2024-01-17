from math import log,e,pi
import nozzle as nz;
from scipy.optimize import fsolve,root
import  fluids as f
import numpy as np;
from sim import fdvt;
import rocket;
import matplotlib.pyplot as plt;
import convert as c
#table of density vs temp and pressure of air from https://www.engineeringtoolbox.com/air-temperature-pressure-density-d_771.html
density = [['p',30,40,50,60,70,80,90,100,120,140,150,200,250,300,400,500,600],[0,0.081, 0.080, 0.078, 0.076, 0.075, 0.074, 0.072, 0.071, 0.069, 0.066, 0.065, 0.060, 0.056, 0.052, 0.046, 0.041, 0.038], [5, 0.109, 0.107, 0.105, 0.102, 0.101, 0.099, 0.097, 0.095, 0.092, 0.089, 0.087, 0.081, 0.075, 0.070, 0.062, 0.056, 0.05],[10, 0.136, 0.134, 0.131, 0.128, 0.126, 0.124, 0.121, 0.119, 0.115, 0.111, 0.109, 0.101, 0.094, 0.088, 0.078, 0.070, 0.063],[20, 0.192, 0.188, 0.185, 0.180, 0.177, 0.174, 0.171, 0.168, 0.162, 0.156, 0.154, 0.142, 0.132, 0.123, 0.109, 0.098, 0.089],[30, 0.247, 0.242, 0.238, 0.232, 0.228, 0.224, 0.220, 0.216, 0.208, 0.201, 0.198, 0.183, 0.170, 0.159, 0.141, 0.126, 0.114],[40, 0.302, 0.295, 0.291, 0.284, 0.279, 0.274, 0.269, 0.264, 0.255, 0.246, 0.242, 0.225, 0.208, 0.195, 0.172, 0.154, 0.140],[50, 0.357, 0.350, 0.344, 0.336, 0.330, 0.324, 0.318, 0.312, 0.302, 0.291, 0.287, 0.265, 0.246, 0.23, 0.203, 0.182, 0.165],[60, 0.412, 0.404, 0.397, 0.388, 0.381, 0.374, 0.367, 0.361, 0.348, 0.337, 0.331, 0.306, 0.284, 0.266, 0.235, 0.210, 0.190],[70, 0.467, 0.458, 0.451, 0.440, 0.432, 0.424, 0.416, 0.409, 0.395, 0.382, 0.375, 0.347, 0.322, 0.301, 0.266, 0.238, 0.216],[80, 0.522, 0.512, 0.504, 0.492, 0.483, 0.474, 0.465, 0.457, 0.441, 0.427, 0.420, 0.388, 0.361, 0.337, 0.298, 0.267, 0.241],[90, 0.578, 0.566, 0.557, 0.544, 0.534, 0.524, 0.515, 0.505, 0.488, 0.472, 0.464, 0.429, 0.399, 0.372, 0.329, 0.295, 0.267],[100, 0.633, 0.620, 0.610, 0.596, 0.585, 0.574, 0.564, 0.554, 0.535, 0.517, 0.508, 0.470, 0.437, 0.408, 0.360, 0.323, 0.292],[120, 0.743, 0.728, 0.717, 0.700, 0.687, 0.674, 0.662, 0.650, 0.628, 0.607, 0.597, 0.552, 0.513, 0.479, 0.423, 0.379, 0.343],[140, 0.853, 0.836, 0.823, 0.804, 0.789, 0.774, 0.760, 0.747, 0.721, 0.697, 0.686, 0.634, 0.589, 0.550, 0.486, 0.436, 0.394],[150, 0.909, 0.890, 0.876, 0.856, 0.840, 0.824, 0.809, 0.795, 0.768, 0.742, 0.730, 0.675, 0.627, 0.586, 0.518, 0.464, 0.420],[200, 1.185, 1.161, 1.142, 1.116, 1.095, 1.075, 1.055, 1.036, 1.001, 0.967, 0.951, 0.879, 0.817, 0.764, 0.675, 0.604, 0.547],[250, 1.46, 1.431, 1.408, 1.376, 1.350, 1.325, 1.301, 1.278, 1.234, 1.193, 1.173, 1.084, 1.008, 0.941, 0.832, 0.745, 0.675],[300, 1.736, 1.702, 1.674, 1.636, 1.605, 1.575, 1.547, 1.519, 1.467, 1.418, 1.395, 1.289, 1.198, 1.119, 0.989, 0.886, 0.802],[400, 2.29, 2.24, 2.21, 2.16, 2.12, 2.08, 2.04, 2.00, 1.933, 1.868, 1.838, 1.698, 1.579, 1.475, 1.303, 1.167, 1.057],[500, 2.84, 2.78, 2.74, 2.68, 2.63, 2.58, 2.53, 2.48, 2.4, 2.32, 2.28, 2.11, 1.959, 1.83, 1.618, 1.449, 1.312],[700, 3.94, 3.86, 3.8, 3.72, 3.65, 3.58, 3.51, 3.45, 3.33, 3.22, 3.17, 2.93, 2.72, 2.54, 2.25, 2.01, 1.822],[800, 4.49, 4.4, 4.33, 4.24, 4.16, 4.08, 4.00, 3.93, 3.8, 3.67, 3.61, 3.34, 3.1, 2.9, 2.56, 2.29, 2.08],[900, 5.05, 4.95, 4.87, 4.76, 4.67, 4.58, 4.50, 4.42, 4.26, 4.12, 4.05, 3.75, 3.48, 3.25, 2.87, 2.58, 2.33],[1000, 5.6, 5.49, 5.4, 5.28, 5.18, 5.08, 4.99, 4.9, 4.73, 4.57, 4.5, 4.16, 3.86, 3.61, 3.19, 2.86, 2.59]]

rho = lambda c,b, p: c*p+b;
pressfromrho = lambda c,b,rho: (rho-b)/c;
ctof=lambda c: 9*c/5+32;
ftoc=lambda f:(f-32)*5/9;
patopsi= lambda pa: pa/6894.757;
psitopat= lambda psi:psi*6894.757;
in2tom2= lambda i: i*0.00064516
m2toin2 = lambda i:i/0.00064516

def getB(temp):
	index=0;
	for i in range(1,len(density[0])):
		if density[0][i] ==temp:
			return density[1][i]
		elif density[0][i]>temp:
			index=i;
			break
	return density[1][index-1]+ ((temp-density[0][index-1])/(density[0][index]-density[0][index-1]))*(density[1][index]-density[1][index-1]);

def getC(temp, press):
	indexTemp=0;
	for i in range(1,len(density[0])):
		if density[0][i]>temp:
			indexTemp=i;
			break
	indexPress=0;
	for i in range(1,len(density)):
		if density[i][0]>press:
			indexPress=i;
			break
	tratio=(temp-density[0][indexTemp-1])/(density[0][indexTemp]-density[0][indexTemp-1])
	dphi=  density[indexPress][indexTemp-1]+ tratio*(density[indexPress][indexTemp]-density[indexPress][indexTemp-1])
	dplo=density[indexPress-1][indexTemp-1]+ tratio*(density[indexPress-1][indexTemp]-density[indexPress-1][indexTemp-1])
	return (dphi-dplo)/(density[indexPress][0]-density[indexPress-1][0])

def getDensityUS(temp, press):
	"""temp in F, press in psi, lb/ft3 out """
	# print(temp)
	# print(press)
	if press>1000 or press<0:
		raise ValueError("pressure out of range; must be between 0 psi and 1000 psi");
	elif temp>600 or temp<30:
		raise ValueError("temperature out of range; temp must be between 30 F and 600 F");
	else:
		indexTemp=0;
		for i in range(1,len(density[0])):
			if density[0][i]>temp:
				indexTemp=i;
				break
		indexPress=0;
		for i in range(1,len(density)):
			if density[i][0]>press:
				indexPress=i;
				break
		tratio=(temp-density[0][indexTemp-1])/(density[0][indexTemp]-density[0][indexTemp-1])
		dphi=  density[indexPress][indexTemp-1]+ tratio*(density[indexPress][indexTemp]-density[indexPress][indexTemp-1])
		# print(dphi)
		dplo=density[indexPress-1][indexTemp-1]+ tratio*(density[indexPress-1][indexTemp]-density[indexPress-1][indexTemp-1])
		# print(dplo)
		# print(press-density[indexPress-1][0])
		return dplo+ (press-density[indexPress-1][0])/(density[indexPress][0]-density[indexPress-1][0])*(dphi-dplo);

def getDensitySI(temp,press):
	""" temp in c, press in Pa, kg/m3 out"""
	return getDensityUS(ctof(temp),patopsi(press))*16.01846337396

def getVa(temp,pressin, Ain, Ath, vth):
	"""use si: c, Pa, m2, m2, m/s"""
	# c=getC(ctof(temp),patopsi(pressin))
	a=f.vSound(0)
	rhoin=getDensitySI(temp,pressin)
	# const= (((rhoin**2)*(Ain**2))/((vth**2)*(Ath**2)))*((vth**2)+ log(rhoin,e)*a*a)
	# coeff=((rhoin**2)*(Ain**2))/((vth**2)*(Ath**2))*a*a
	# print(coeff)
	# print(const)
	return fsolve( lambda x: (x**2)/2+log((Ain*x*rhoin/Ath/vth),e)*a*a -(vth**2)/2-log(rhoin,e)*a*a,.05 )[0];

def rhob(rhoa, Aa, va,Ab,Vb):
	return rhoa*va*Aa/Vb/Ab;

def mdot(rho,A,v):
	return rho*A*v;

def getVexss(temp, vth, rhoth):

	rhoex=getDensitySI(temp,0)
	# print(rhoex)
	# return (2*((vth**2)/2+ (1/getC(temp,0))*(-log(rhoex,e)+log(rhoth,e))))**.5
	return (2*((vth**2)/2+ (f.vSound(0)**2)*(-log(rhoex,e)+log(rhoth,e))))**.5  #use: best one
	# return fsolve(lambda x: log(x,e)-(x**2)/(2*f.vSound(0)**2)+(vth**2)/(2*f.vSound(0)**2)-log(vth,e)+log(ath,e)-log(aex,e),300)[0]
	# return mdot(rhoth,ath,vth)/rhoex/aex

def eqns (T,ac,vc,aa,presa,temp):
	a = f.vSound(0)
	vv = getDensitySI(temp,0)*vc*ac/getDensitySI(temp,presa)/aa
	# return [T-vc*vc*ac*getDensitySI(temp,0), log(vv,e)-(vv**2)/2/(a**2) +(vc**2)/(2*a**2)-log(vc,e)-log(ac,e)+log(aa,e), (vc**2)-(vv**2)+2*(a**2)*(log(getDensitySI(temp,0),e)-log(getDensitySI(temp,presa),e)) ]
	return [T-vc*vc*ac*getDensitySI(temp,0), log(vv,e)-(vv**2)/2/(a**2) +(vc**2)/(2*a**2)-log(vc,e)-log(ac,e)+log(aa,e) ]

def test():
	temp=60
	pressa=psitopat(1)
	Aa=in2tom2(9*pi)
	Ab=in2tom2(.125*pi)
	Aex=in2tom2(4*pi)
	vb=f.vSound(0)
	va=getVa(temp,pressa,Aa,Ab,vb)
	rha=getDensitySI(temp,pressa)
	rhb=rhob(rha,Aa,va,Ab,vb)
	vexss=getVexss(temp,vb,rhb)
	Aexss=mdot(rha,Aa,va)/vexss/getDensitySI(temp,0)
	print(rha)
	print(rhb)
	print(va)
	print(mdot(rha,Aa,va))
	print(vexss)
	print(m2toin2(Aexss)) 
	print("idealized thrust (supersonic exit velocity): "+ str(mdot(rha,Aa,va)*vexss))
	# ccmair=mdot(rha,Aa,va)/getDensitySI(temp,psitopat(100))*20
	# print("requires "+ str(ccmair)+" cubic meters of air for 20 seconds of opperation")
	# print("length of required 6 in diameter air canister: "+ str((ccmair/in2tom2((3**2)*pi))))
	rock = rocket.Rocket();

	cdlist=np.loadtxt('jpltestcd.csv', skiprows=0, delimiter=',', unpack=True);

	rock.setCd(cdlist);
	rock.setDiameter(8*c.inch);

	fvt =fdvt(.3,.01,.5,rock)
	# print(fvt)

	# thrust = np.linspace(54,0,5);
	thrust = fvt[:,2];
	out=[]
	tmdot=[];
	for i in range(len(thrust)):
		func = lambda x: eqns(thrust[i],x[1],x[0],Aa,pressa,temp)
		root=0
		found=False;
		for j in range(1,1000):
			for k in range(1,1000):
				try: 
					x0=(j,in2tom2(k))
					root = fsolve(func,x0)
				except:
					1==1
				else:
					out.append(root)
					found=True
					break
			if found:
				tmdot.append(mdot(getDensitySI(temp,0),root[1], root[0]))
				break
	out = np.array(out)
	plt.plot(fvt[:,0],fvt[:,2], label='fd v t');
	plt.show()
	plt.plot(fvt[:,2], out[:,0], label='vc vs fd', color='blue')
	plt.tick_params(axis ='y', labelcolor = 'blue') 
	plt.legend()
	plt.twinx()
	plt.plot(fvt[:,2], tmdot[:], color='red', label = 'mdot vs fd')
	plt.tick_params(axis ='y', labelcolor = 'red') 
	plt.legend(loc='lower right')
	plt.show()
	mass = sum(tmdot)*fvt[1][0]
	vol = mass/getDensitySI(temp,psitopat(100))
	print(str(mass)+" kg of air = "+ str(vol)+" m^3 air used in " +str(fvt[-1][0])+ " s")


test();
