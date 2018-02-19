%366L homework 3
%Eric Johnson
clc; clear all;


%Problem 2
% 
% UTC = [14, 58, 9];
% dUT1 = 0.184798;
% UT1 = UTC;
% UT1(3) = UT1(3) + dUT1
% TUT1 = UT1(1)*60*60 + UT1(2)*60 + UT1(3);
% 

%Problem 3
clear all;
fprintf('Problem 3 \n');
% October 24 2016 06:30:30 TAI

numOfDays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
Month = 10;
DOM = 30;
Year = 2016;
DOY = sum(numOfDays(1:10)+DOM);
Hour = 6;
Min = 30;
Sec = 30;

fprintf('Date Y/DOY/H/M/S: %4.0d:%3.0d:%2.0d:%2.0d:%2.0d \n', Year, DOY, Hour, Min, Sec);

JD = 367*Year - round((7*(Year+round((Month+9)/12)))/4) + round((275*Month)/9) + DOM + 1721013.5 + (1/24)*(Hour +(1/60)*(Min+(Sec/60)));

fprintf('Julian Date: %f \n',JD);
MJD = JD - 2400000.5;
fprintf('Modified Julian Date %f \n', MJD);
theta = Hour*6 + Min*(15/60) + Sec*(15/3600);
fprintf('Date Y/M/D/R: %4.0d:%2.0d:%2.0d:%2.0d \n', Year, Month, DOM, theta);

%problem 4 
clear all;
fprintf('\nProblem 4 \n');
dUT1 = 0.184798;

%observation 1
UTC1 = 2458165.4375;
%observation 2
TT2 = 2458165.4375;
UTC2 = TT2 - 32.184 - dUT1;

dObserv = UTC1 - UTC2;

rt1 = [-1140.697, 6126.252, -926.199]; 
vt1 = [-5.706065, -1.934816, -5.7701134];

rt2 = [-1350.501, 6048.007, -1138.605];
vt2 = [-5.632593, -2.293872, -5.709216];


t = 0:.2:dObserv;
y = [rt1 , vt1];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-20);
[T,Y] = ode45('func',t,y, options);

rt1atT2 = [Y(length(t),1), Y(length(t),2), Y(length(t),3)];
vt1atT2 = [Y(length(t),4), Y(length(t),5), Y(length(t),6)];
fprintf('The Position and velocity of operator A spacecraft at t 2458165.4375 TT is \n');
fprintf('r = %g i %g j %g k [km] \n',rt1atT2);
fprintf('v = %g i %g j %g k [km/s] \n',vt1atT2);

dr = norm((rt2 - rt1atT2) ./ rt2)*100; 
dv = norm((vt2 - vt1atT2) ./ vt2)*100;
fprintf('No because there is a %f percent error in position \n', dr);
fprintf('and a %f percent error in velocity \n', dv);

%problem 5













