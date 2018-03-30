%Eric Johnson
%Homework 6 ASE 366L
clc; clear all; close all;

Re =  6378.1363;
mue = 398600.4415;
J2 = 0.0010826267;
J3 = -0.0000025327;
J4 = -0.0000016196;
%Problem 2
%orbit charateristics
a = 7303.14;
e = .003;
i = 98;
OMEGAsec2 = -(3/2)*((sqrt(mue)*(Re^2)*J2)/((a^(7/2))*(1-e^2)))*cosd(i);
n = sqrt(mue/(a^3));
p = a*(1 - e^2);
OMEGAsec4 = (15*J4*n*cosd(i))/(32*p^4)*(8+12*e^12 - (14 +21*e^2)*sind(i)^2); 
fprintf('Problem 2 \n');
fprintf('A.)OMEGA dot due to J2 is %e [rad/s]\n', OMEGAsec2);
fprintf('B.)OMEGA dot due to J4 is %e [rad/s]\n', OMEGAsec4);
percentDiff = OMEGAsec4/OMEGAsec2*100;
fprintf('C.)the percent difference between OMEGA dot J2 and J3 is %e \n', percentDiff);
fprintf('The perturbation due to J4 relative to J2 is essentially negligable\n');
clear all;
%constants
Re = 6378.1363;
mu = 398600.4415;
J2 = 0.0010826267;
J3 = -0.0000025327;
J4 = -0.0000016196;

%problem 5

fprintf('\nProblem 5: \n');
%r = [0.009, 4595.737, 4595.731];
r = [6092.032, 2487.062, 2487.062];

ap = aJ2J3(mu, Re, J2, J3, r);

fprintf('ap = %g i %g j %g k [km/s] \n', ap );

%Problem 6 
clear all;
%constants
Re = 6378.1363;
mu = 398600.4415;
J2 = 0.0010826267;
J3 = -0.0000025327;
J4 = -0.0000016196;

fprintf('\nProblem 6: \n');

%part a
a = 7000;
e = 0.01;
i = (pi/4);
OMEGA = 0;
omega = -(pi/2);
theta = (2*pi)/3;

[rVec,vVec] = kep2cart_Johnson(a,e,mu,i,OMEGA,omega,theta);

y = [rVec , vVec];
t = 0:1:86400;
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-20);
[T,Y] = ode45('funcJ2J3',t, y, options);
rf = [Y(length(Y),1), Y(length(Y),2), Y(length(Y),3)];
fprintf('A.) r = %g i %g j %g k [km] \n', rf);

%part b

a = 6800.0; 
e = 0.03;
i = 1.309;
OMEGA = 3.49;
omega = 5.24;
theta = 0;

[rVec,vVec] = kep2cart_Johnson(a,e,mu,i,OMEGA,omega,theta);

y = [rVec , vVec];
t = 0:60:86400;
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-20);
[T,Y] = ode45('funcJ2J3',t, y, options);
rf = [Y(length(Y),1), Y(length(Y),2), Y(length(Y),3)];
fprintf('B.) r = %g i %g j %g k [km] \n', rf);

figure(1)
subplot(3,1,1)
title('Position Vector Over Time');
xlabel('time [s]');
ylabel('radius [km]');
hold on
r = (Y(:,1).^2 + Y(:,2).^2 + Y(:,3).^2 ).^(1/2);
plot(t, r)

subplot(3,1,2)
title('Velocity Vector Over Time');
xlabel('time [s]');
ylabel('velocity [km/s]');
hold on
v = ( Y(:,4).^2 + Y(:,5).^2 + Y(:,6).^2 ).^(1/2);
plot(t, v)


ap = aJ2J3(mu, Re, J2, J3, r);

rt = [Y(:,1), Y(:,2), Y(:,3)];
for i=1:length(rt)
    
    
    ap = aJ2J3(mu, Re, J2, J3, rt(i,:));
    a(1) = -(mu/(norm(rt(i,:))^3))*rt(i,1) + ap(1);
    a(2) = -(mu/(norm(rt(i,:))^3))*rt(i,2) + ap(2);
    a(3) = -(mu/(norm(rt(i,:))^3))*rt(i,3) + ap(3);
    at(i) = norm(a); 
end


subplot(3,1,3)
title('Acceleration Vector Over Time');
xlabel('time [s]');
ylabel('acceleration [km/s^2]');
hold on
plot(t, at)






