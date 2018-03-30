function ydot = funcDrag(t,y)

%y = [x, y, z, vx, vy, vz]
%ydot = [vx, vy, vz, ax, ay, az]
mu = 398600.4415;
rvec = [y(1), y(2), y(3)];
vvec = [y(4), y(5), y(6)];
ad = drag(vvec, rvec);
r = sqrt( y(1).^2 + y(2).^2 + y(3).^2 );



ydot(1) = y(4);
ydot(2) = y(5);
ydot(3) = y(6); 

ydot(4) = -(mu/(r^3))*y(1) + ad(1);
ydot(5) = -(mu/(r^3))*y(2) + ad(2);
ydot(6) = -(mu/(r^3))*y(3) + ad(3);

ydot = ydot';




 