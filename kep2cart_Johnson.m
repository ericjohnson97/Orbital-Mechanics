function [rVec,vVec] = kep2cart_Johnson(a,e,mu,i,OMEGA,omega,theta)
r = (a*( 1- e^2 ))/(1+e*cos(theta));
X = r*(cos(OMEGA)*cos(omega+theta)-sin(OMEGA)*sin(omega+theta)*cos(i));
Y = r*(sin(OMEGA)*cos(omega+theta)+cos(OMEGA)*sin(omega+theta)*cos(i));
Z=r*(sin(i)*sin(omega+theta));

rVec = [ X, Y, Z];

h = sqrt(mu*a*(1-e^2));
p = a*(1-e^2);

Vx = ( (X*h*e)/(r*p) )*sin(theta) - (h/r)*( cos(OMEGA)*sin(omega+theta)+sin(OMEGA)*cos(omega+theta)*cos(i) );
Vy = ( (Y*h*e)/(r*p) )*sin(theta) - (h/r)*( sin(OMEGA)*sin(omega+theta)-cos(OMEGA)*cos(omega+theta)*cos(i) );
Vz = ( (Z*h*e)/(r*p) )*sin(theta) + (h/r)*sin(i)*cos(omega+theta);

vVec = [ Vx, Vy, Vz];
end