%Eric Johnson 
%Homework 8 ASE 366L

clc; clear all; close all;

%Problem 4

r = [2781.000, 5318.000, -5629.000];
rSun = [2379260.000, 148079334.000, -1009936.000];

am = .01;
Cr = 1.5;
gamma = 1;
asrp = accelSRP(gamma, Cr, am, rSun);

fprintf('Problem 4.) \n');
fprintf('The acceleration from the solar radiation pressure is: \n%g i %g j %g k [km/s^2]\n', asrp);

%Problem 5

RSun = 432288;
REarth = 6378.1363;
muEartrh = 398600.4415;
muSun = 1.3271244E11;
omegaEarth = 7.2921158553E-5;
J2 = 0.0010826267;
J3 = -0.0000025327;

a = 6800; 
e = 0.005; 
i = 71;
OMEGA = 300;
omega = 78;
theta = 0;
[rVec,vVec] = kep2cart_Johnson(a,e,muEartrh,i,OMEGA,omega,theta);

y = [rVec , vVec];
t = 0:10:86400;
options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-30);
[T,Y] = ode45('funcCombined',t, y, options);
Y(length(Y),1:3);


figure(3)
subplot(3,1,1)
title('Orbit Radius i Component');
xlabel('Time [Days]');
ylabel('Radius [km]');
hold on 
plot(t, Y(:,1) )

subplot(3,1,2)
title('Orbit Radius j Component');
xlabel('Time [Days]');
ylabel('Radius [km]');
hold on 
plot(t, Y(:,2) )

subplot(3,1,3)
title('Orbit Radius k Component');
xlabel('Time [Days]');
ylabel('Radius [km]');
hold on 
plot(t, Y(:,3) )
[T,Ytrue] = ode45('func',t, y, options);
poseErr = norm(Y(length(Y),1:3) - Ytrue(length(Ytrue),1:3));
fprintf('Problem 5.) \n');
fprintf('The magnitude of the position error is %f [km] \n', poseErr);










