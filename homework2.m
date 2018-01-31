%Eric Johnson 366L
%Homework 2
clc; clear all;
close all;

mu = 398600.4415;

a = 26000;
e = .72;
i = 75;
OMEGA = 90;
omega = -90;
tp = 0;

%problem 1
[rVec,vVec] = kep2cart_Johnson(a,e,mu,i,OMEGA,omega,0);

fprintf('Problem 1.) \n');
fprintf('r = %g i %g j %g k \n',rVec);
fprintf('v = %g i %g j %g k \n',vVec);

%problem 2

theta = 0;
%eccentric anomally - angle from center to spacecraft
E = acos( (1/e) - (norm(rVec)/(e*a) ) );
Me = E - e*sin(E);

T = (2*pi) * sqrt( (a^3)/(mu)  );
tp = 0 - ( (T*Me) / (2*pi) );

r(1) = norm(rVec);
v(1) = norm(vVec);
theta(1) = theta;
t = 0:20:3*86400;
for j=1:length(t)
    Me = ( (2*pi) / T ) * ( t(j) - tp );
    for i=1:100
        fEo = E - e*sin(E) - Me;
        fEoPrime = 1 - e*cos(E);
        E1 = E - (fEo/fEoPrime);
        if ( E1 - E < 10^10)
            E = E1;
            break;
        end
        E = E1;
    end
    theta(j) = acos( (cos(E) - e) / ( 1 - e*cos(E) ));
    r(j) = ( a*(1 - e^2) ) / ( 1 + e*cos(theta(j)) );
    v(j) = sqrt( (mu/a) + ((2*mu)/r(j)) );
end
fprintf('\nProblem 2.) \n');
fprintf('the true anomally is %f [rad] at tf\n',theta(length(t)));
fprintf('the eccentric anomally is %f [rad] at tf \n',E);

%problem 3

fprintf('\nProblem 3.) \n');
[rVectf,vVectf] = kep2cart_Johnson(a,e,mu,i,OMEGA,omega,theta(length(t)));
fprintf('the position vector is %g i %g j %g k [km] at tf\n',rVectf);
fprintf('the velocity vector is %g i %g j %g k [km] at tf \n',vVectf);

%problem 4 

ht = 20; %step size for t
y = [rVec , vVec];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-20);
[T,Y] = ode45('func',t,y, options);

rt = sqrt(Y(:,1).^2 + Y(:,2).^2 + Y(:,3).^2);
vt = sqrt(Y(:,4).^2 + Y(:,5).^2 + Y(:,6).^2);
ax = -(mu./(rt.^3)).*Y(:,1);
ay = -(mu./(rt.^3)).*Y(:,2);
az = -(mu./(rt.^3)).*Y(:,3);
at = sqrt(ax.^2 + ay.^2 + az.^2);

fprintf('\nProblem 4.) \n');
fprintf('relative tolerance = 1e-13 and absolute tolerance = 1e-20 \n'); 

Eo = ((vt(1)^2)/2) - (mu/rt(1));
Et = ((vt.^2)/2) - (mu./rt);

figure(1)
title('Change in Specific Energy');
xlabel(' Time [s]');
ylabel('Delta Specific Energy [km^2/s^2]');
hold on
plot(t, Et - Eo);
fprintf('the error inspecific energy increases with each iteration due to propagted error. \n');
fprintf('the relative error increases when the velocity increses at pariapse since the change in positions are larger \n');

%problem 5

fprintf('\nProblem 5.) \n');
posErr = (rt(length(t)) - r(length(t))) / r(length(t));
velErr = (vt(length(t)) - v(length(t))) / v(length(t));
fprintf('the position percent error is %f at tf\n', posErr); 
fprintf('the velocity percent error is %f at tf\n', velErr);

%problem 6
for i=1:length(t)
    [rVectf,vVectf] = kep2cart_Johnson(a,e,mu,i,OMEGA,omega,theta(length(t)));
    riErr(i) = (Y(i,1) - rVec(1)) / rVec(1);
    rjErr(i) = (Y(i,2) - rVec(2)) / rVec(2);
    rkErr(i) = (Y(i,3) - rVec(3)) / rVec(3);
    viErr(i) = (Y(i,4) - vVec(1)) / vVec(1);
    vjErr(i) = (Y(i,5) - vVec(2)) / vVec(2);
    vkErr(i) = (Y(i,6) - vVec(3)) / vVec(3);
end
figure(2)
subplot(3,1,1)
title('Change in position error in i direction');
xlabel(' Time [s]');
ylabel('Percent error');
hold on
plot(t,riErr);
subplot(3,1,2)
title('Change in position error in j direction');
xlabel(' Time [s]');
ylabel('Percent error');
hold on
plot(t,rjErr);
subplot(3,1,3)
title('Change in position error in k direction');
xlabel(' Time [s]');
ylabel('Percent error');
hold on
plot(t,rkErr);

%velocity 
figure(3)
subplot(3,1,1)
title('Change in velocity error in i direction');
xlabel(' Time [s]');
ylabel('Percent error');
hold on
plot(t,viErr);
subplot(3,1,2)
title('Change in velocity error in j direction');
xlabel(' Time [s]');
ylabel('Percent error');
hold on
plot(t,vjErr);
subplot(3,1,3)
title('Change in velocity error in k direction');
xlabel(' Time [s]');
ylabel('Percent error');
hold on
plot(t,vkErr);

%problem 7















