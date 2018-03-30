%Eric Johnson 
%Homework 7 ASE 366L

clc; clear all; close all;

rEarth = 6378.1363;
muEarth = 398600.4415;
J2 = 0.0010826267;
J3 = -0.0000025327;
J4 = -0.0000016196;
CD = 2.0;
omegaEarth = 7.2921158553E-5;


altT = [0; 25; 30; 40; 50; 60; 70; 80; 90; 100; 110; 120; 130; 140; 150; 180; 200; 250; 300; 350; 400; 450; 500; 600; 700; 800; 900; 1000];
rho0T = [1.225; 3.899E-2; 1.774E-2; 3.972E-3; 1.057E-3; 3.206E-4; 8.77E-5; 1.905E-5; 3.396E-6; 5.297E-7; 9.661E-8; 2.438E-8; 8.484E-9; 3.845E-9; 2.070E-9; 5.464E-10; 2.789E-10; 7.248E-11; 2.418E-11; 9.518E-12; 3.725E-12; 1.585E-12; 6.967E-13; 1.454E-13; 3.614E-14; 1.170E-14; 5.245E-15; 3.019E-15];
scaleHT = [7.249; 6.349; 6.682; 7.554; 8.382; 7.714; 6.549; 5.799; 5.382; 5.877; 7.263; 9.473; 12.636; 16.149; 22.523; 29.740; 37.105; 45.546; 53.628; 53.298; 58.515; 60.828; 63.822; 71.835; 88.667; 124.64; 181.05; 268];

EAM = [altT, rho0T, scaleHT];
alt=0:10:1000;
for i=1:length(alt)
    h = alt(i);
    for ii=1:length(EAM)
       if h < EAM(ii,1)
           h0 = EAM(ii-1,1);
           H = EAM(ii-1,3);
           rho0 = EAM(ii-1,2);
           break;
       end
    end
    rho(i) = rho0*exp(-(h-h0)/H);
end

semilogx(rho, alt)
title('Altitude vs Density');
xlabel('Density [kg/m^3]');
ylabel('Altitude [km]');

clear all;

rEarth = 6378.1363;
muEarth = 398600.4415;
J2 = 0.0010826267;
J3 = -0.0000025327;
J4 = -0.0000016196;
CD = 2.0;
omegaEarth = 7.2921158553E-5;

%orbit characteristics

a=7300;
e=.01;
i=45*(pi/180);
OMEGA = 0;
omega = 0;
theta = 0;

[rVec,vVec] = kep2cart_Johnson(a,e,muEarth,i,OMEGA,omega,theta);

y = [rVec , vVec];
secPDay = 86400;
t = .5:.5*secPDay:100*secPDay;
options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-30);
[T,Y] = ode45('funcDrag',t, y, options);

Y

 E0 = (norm(Y(1,4:6))^2/2) - (muEarth./norm(Y(1,1:3)));
 for ii=1:length(Y)
    E(ii) = (norm(Y(ii,4:6)).^2/2) - (muEarth./norm(Y(ii,1:3)));
 end
dE = E0 - E
figure(2)
title('Difference in Specific Energy');
xlabel('Time [Days]');
ylabel('Specifc Energy [Km^2 / s^2]');
hold on 
plot(t/secPDay, dE)

for ii=1:length(Y)
    r(ii) = norm(Y(ii,1:3));
    [a(ii),e(ii),i,OMEGA,omega,theta] = cart2kep_Johnson(Y(ii,1:3),Y(ii,4:6),muEarth);
    
    rp(ii) = a(ii)*(1-e(ii));
    ra(ii) = a(ii)*(1+e(ii));
end
figure(3)
subplot(3,1,1)
title('Orbit Radius');
xlabel('Time [Days]');
ylabel('Radius [km]');
hold on 
plot(t/secPDay, r)

subplot(3,1,2)
title('Radius of Pariapse');
xlabel('Time [Days]');
ylabel('Radius [Km]');
hold on 
plot(t/secPDay, rp)

subplot(3,1,3)
title('Radius of Apoapsis');
xlabel('Time [Days]');
ylabel('Radius [Km]');
hold on 
plot(t/secPDay, ra)

figure(4)
subplot(2,1,1)
title('Semi Major Axis Over Time');
xlabel('Time [Days]');
ylabel('Semi Major axis [Km]');
hold on 
plot(t/secPDay, a)
 
subplot(2,1,2)
title('Eccentricity Over Time');
xlabel('Time [Days]');
ylabel('Semi Major axis [Km]');
hold on 
plot(t/secPDay, e)





