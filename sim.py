import convert as c;
import fluids as f;
import phys as p;
import numpy as np;
import propeller as pr;
from math import pi;
import nozzle as nz;
import rocket;
import matplotlib.pyplot as plt;

def fdVM( mi,mf, n,rock,y=0):
	""" calculates the drag force as a function of the mach number at a given altitude

	inputs: mi-initial mach number, mf-final mach number, n- numer of steps, rock - the rocket to test, y- altitude in m
	output: the curve of Fd vs V"""
	v = np.linspace(mi*f.vSound(y), mf*f.vSound(y),n)
	dat =  np.zeros((2, len(v)));
	dat[0] = v;
	for i in range(len(v)):
		dat[1][i]=rock.getDrag(v[i],y)
	return dat

def vvst(mi,mf,dt):
	""" calculates a curve of velocity vs time"""
	n = abs(int((mf*f.vSound(0)-mi*f.vSound(0))/(9.81*dt)))
	temp =np.linspace(mi*f.vSound(0), mf*f.vSound(0),n)
	out= np.zeros((2, len(temp)));
	out[1]=temp
	out[0]=np.linspace(0,abs((mf*f.vSound(0)-mi*f.vSound(0))/(9.81)),n)
	return out

def fdvt (mi,mf,dt,rock):
	""" calculates fd vs time"""
	vvt=vvst(mi,mf,dt);
	# print(vvt)
	out = [];
	for i in range(len(vvt[0])):
		out.append([vvt[0][i],vvt[1][i],rock.getDrag(vvt[1][i],0)])
	return np.array(out)

# build a rocket
def test ():
	testRocket = rocket.Rocket();
		
		# import drag data from file

	cdlist=np.loadtxt('jpltestcd.csv', skiprows=0, delimiter=',', unpack=True);

	testRocket.setCd(cdlist);
	testRocket.setDiameter(8*c.inch);
	noz = nz.Nozzle();
	noz.setA(12*c.inch**2,1.96*pi*c.inch**2,1.96*pi*c.inch**2)
	noz.setL(8*c.inch)
	testRocket.setNozzle(noz)
	# print(cdlist)
	z=0;
	curve = fdVM(0.01,.3,1000,testRocket, y=z);

	# testRocket.setAex(9*pi*c.inch**2);
	# testRocket.setAin(12*c.inch**2);
	# plt.plot(cdlist[0],cdlist[1]);
	# plt.show();
	out = [];
	for i in range(len(curve[0])):
		vout = testRocket.getMinVout(curve[0][i],z, curve[1][i] )
		out.append([f.mach(curve[0][i],z),curve[1][i], vout, vout-curve[0][i], testRocket.getPower(curve[0][i],vout,curve[1][i]) ])
	out = np.array(out)

	plo=False

	if plo:
		plt.plot(f.mach(curve[0],z),curve[1], color = 'blue', label = 'Fd (N)');
		# plt.plot(out[:,0], out[:,3], color='green', label='min delta V (m/s)')
		plt.tick_params(axis ='y', labelcolor = 'blue') 
		plt.legend();
		plt.twinx()
		plt.plot(out[:,0], out[:,2], color='red', label='min Vout (m/s)')
		plt.tick_params(axis ='y', labelcolor = 'red') 
		plt.legend(loc='lower right')
		plt.show();
		plt.plot(out[:,0], out[:,3], color='green', label='min delta V (m/s)')
		plt.legend();
		plt.show();
		plt.plot(out[:,1], out[:,3], color='green', label='min delta V (m/s)')
		plt.tick_params(axis ='y', labelcolor = 'green') 
		plt.legend(loc = 'upper center');
		plt.twinx()
		plt.plot(out[:,1], out[:,2], color='red', label='min Vout (m/s)')
		plt.tick_params(axis ='y', labelcolor = 'red') 
		plt.legend(loc='lower center')
		plt.show();
		plt.plot(out[:,1],out[:,4], color='c', label='power (w)')
		plt.show()

	
	print(testRocket.getDrag(.3*f.vSound(z),z))

	v= f.vSound(z)*.3;
	pres=p.correctedPresure(v,z);
	thrust=testRocket.getDrag(v,z);
	rvout =testRocket.getMinVout(v,z,thrust)
	print(p.correctedDensity(v,z))
	print(rvout);
	print(rvout*p.correctedDensity(rvout,z)*1.96*pi*c.inch**2)
	nvout=noz.getVout(v,z,thrust)
	print(nvout)
	print(nvout*p.correctedDensity(nvout,z)*1.96*pi*c.inch**2)
	print(testRocket.getPower(v,testRocket.getMinVout(v,z,thrust),z))
	print()
	print(v)
	mdot = v*p.correctedDensity(v,z)*12*c.inch**2
	vo=testRocket.noz.vouttank(700,273,mdot)
	print(vo)
	print(mdot/vo/testRocket.noz.aex)
	# print(f.thrust(v,z,12*c.inch**2,9*pi*c.inch**2,testRocket.getMinVout(v,z,thrust),pres))
	# print(28.979407568711018/25)
	# print(57.958815137422036/50)
	print(f.vSound(z))
test();





def test2():
	noz = Nozzle();
	noz.setA(12*c.inch**2,2*pi*c.inch**2,2*pi*c.inch**2)
	M = .3
	z=0
	vinf = M*f.vSound(z)
	print(vinf)
	# print(60*vinf)
	# print(60*vinf-vinf*vinf-cp*p.correctedTemperature(vinf,z))
	print(noz.getVout(vinf,z,60))
	# print(noz.getmdot(vinf,z))
# test2();
# print('hello world')