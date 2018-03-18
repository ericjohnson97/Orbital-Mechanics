function ydot = funcJ2J3(t,y)
Re = 6378.1363;
mu = 398600.4415;
J2 = 0.0010826267;
J3 = -0.0000025327;

%y = [x, y, z, vx, vy, vz]
%ydot = [vx, vy, vz, ax, ay, az]
ydot(1) = y(4);
ydot(2) = y(5);
ydot(3) = y(6);
r(1) = y(1);
r(2) = y(2);
r(3) = y(3);
ap = aJ2J3(mu, Re, J2, J3, r);

ydot(4) = -(mu/(norm(r)^3))*y(1) + ap(1);
ydot(5) = -(mu/(norm(r)^3))*y(2) + ap(2);
ydot(6) = -(mu/(norm(r)^3))*y(3) + ap(3);
ydot = ydot';