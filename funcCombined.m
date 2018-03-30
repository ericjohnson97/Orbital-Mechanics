function ydot = funcCombined(t,y)

%y = [x, y, z, vx, vy, vz]
%ydot = [vx, vy, vz, ax, ay, az]
muEarth = 398600.4415;
REarth = 6378.1363;
J2 = 0.0010826267;
J3 = -0.0000025327;

rvec = [y(1), y(2), y(3)];
vvec = [y(4), y(5), y(6)];
ad = drag(vvec, rvec);
ap = aJ2J3(muEarth, REarth, J2, J3, rvec);
r = sqrt( y(1).^2 + y(2).^2 + y(3).^2 );
[asun, rSun] = accelSun(rvec, t);
gamma =1;
Cr = 1.5;
am = .01;

asrp = accelSRP(gamma, Cr, am, rSun);

ydot(1) = y(4);
ydot(2) = y(5);
ydot(3) = y(6);

ydot(4) = -(muEarth/(r^3))*y(1) + ad(1) + ap(1) + asun(1) + asrp(1);
ydot(5) = -(muEarth/(r^3))*y(2) + ad(2) + ap(2) + asun(2) + asrp(2);
ydot(6) = -(muEarth/(r^3))*y(3) + ad(3) + ap(3) + asun(3) + asrp(3);

ydot = ydot';

end