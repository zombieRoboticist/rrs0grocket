%runID = Simulink.sdi.createRun('cd',"file",'/jpltestcd.csv');
diam=[8,6,1/8]; %[od,tube id,throat id] inches
z0=0;%m
v0=175;%m/s
gama=1.4;
g=9.81;%m/s/s
m=15;%kg
tfin=1.5;

sim("internalprop.slx")
plot(t,v)
plot(t,z)
plot(t,fd)
plot(t,atmo)