%homework 4 
clc; clear all;
close all;

%compare results to 
%rSun0 = [ 26493081.968, -132755488.191, -57556564.835] %test data
%JDut1 = 2451545.0; %test data


JDut1 = 2457793.5;
Tut1 = (JDut1 - 2451545)/36525;

lambdaMSun = 280.46 + 36000.771*Tut1;
MSun = 357.5291092 + 35999.0534*Tut1;

%assume Tut1 is essentially = Ttbd
Ttdb = Tut1;
lambdaE = lambdaMSun + 1.914666471*sind(MSun) + 0.019994643*sind( 2*MSun);

rSunMag = 1.000140612 - 0.016708617*cosd(MSun) - 0.000139589*cosd(2*MSun);
epsilon = 23.439291 - 0.0130042*Ttdb;

rSun = [ rSunMag*cosd(lambdaE), rSunMag*cosd(epsilon)*sind(lambdaE), rSunMag*sind(epsilon)*sind(lambdaE)] ;


deltaUT1 = 67.6439;
TT = Tut1;
zeta= (2306.2181*TT + 0.30188*TT^2 + 0.017998*TT^3)*(1/3600);
phi = (2004.3109*TT - 0.42665*TT^2 - 0.041833*TT^3)*(1/3600);
z = (2306.2181*TT + 1.09468*TT^2 + 0.018203*TT^3)*(1/3600);

R3zeta = [ cosd(zeta), -sind(zeta), 0; sind(zeta), cosd(zeta), 0; 0, 0, 1];
R2phi = [ cosd(-phi), 0, sind(-phi); 0, 1, 0; -sind(-phi), 0, cosd(-phi)];
R3z = [ cosd(z), -sind(z), 0; sind(z), cosd(z), 0; 0, 0, 1];

P = R3zeta*R2phi*R3z;
rSunGCRF = P*rSun';
fprintf('Problem 4.)\n');
fprintf('a.) rSun = %g i %g j %g k [AU] \n', rSunGCRF);

%part b

for i=1:365
    
    JDut1 = 2457793.5 + i;
    Tut1 = (JDut1 - 2451545)/36525;

    lambdaMSun = 280.46 + 36000.771*Tut1;
    MSun = 357.5291092 + 35999.0534*Tut1;

    %assume Tut1 is essentially = Ttbd
    Ttdb = Tut1;
    lambdaE = lambdaMSun + 1.914666471*sind(MSun) + 0.019994643*sind( 2*MSun);

    rSunMag = 1.000140612 - 0.016708617*cosd(MSun) - 0.000139589*cosd(2*MSun);
    epsilon = 23.439291 - 0.0130042*Ttdb;

    rSun = [ rSunMag*cosd(lambdaE), rSunMag*cosd(epsilon)*sind(lambdaE), rSunMag*sind(epsilon)*sind(lambdaE)] ;


    deltaUT1 = 67.6439;
    TT = Tut1;
    zeta= (2306.2181*TT + 0.30188*TT^2 + 0.017998*TT^3)*(1/3600);
    phi = (2004.3109*TT - 0.42665*TT^2 - 0.041833*TT^3)*(1/3600);
    z = (2306.2181*TT + 1.09468*TT^2 + 0.018203*TT^3)*(1/3600);

    R3zeta = [ cosd(zeta), -sind(zeta), 0; sind(zeta), cosd(zeta), 0; 0, 0, 1];
    R2phi = [ cosd(-phi), 0, sind(-phi); 0, 1, 0; -sind(-phi), 0, cosd(-phi)];
    R3z = [ cosd(z), -sind(z), 0; sind(z), cosd(z), 0; 0, 0, 1];

    P = R3zeta*R2phi*R3z;
    rSunGCRF = P*rSun';  
    data(i,:) = rSunGCRF';
    
    
end
figure(1)
subplot(3,1,1);
title('Sun i Position Over Time');
xlabel('Julian Date');
ylabel('i Position [AU]');
hold on
t=[1:365] + 2457793.5; 
plot(t, data(:,1));

subplot(3,1,2);
title('Sun j Position Over Time');
xlabel('Julian Date');
ylabel('j Position [AU]');
hold on
t=[1:365] + 2457793.5; 
plot(t, data(:,2));

subplot(3,1,3);
title('Sun k Position Over Time');
xlabel('Julian Date');
ylabel('k Position [AU]');
hold on
t=[1:365] + 2457793.5; 
plot(t, data(:,3));

%problem 5
clc; clear all;


%spacecraft position
a = 500000.000; 
e = 0.1;
i = 50;
OMEGA = 350;
omega = 90;
theta = 0;
muEarth = 398600.4415;
muSun = 1.3271244*10^11;

[rVec,vVec] = kep2cart_Johnson(a,e,muEarth,i,OMEGA,omega,theta)

%sun position
JDut1 = 2451545.0;

y = [rVec , vVec];
t = 1:86400:365*86400;
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-20);
[T,Y] = ode45('func3',t,y, options);


for ii=1:length(Y)
    
    
    [a,e,i,OMEGA,omega,theta] = cart2kep_Johnson(Y(ii,1:3),Y(ii,4:6),muEarth);
    aT(ii) = a;
    eT(ii) = e;
    iT(ii) = i;
    OMEGAT(ii) = OMEGA;
    omegaT(ii) = omega;
    thetaT(ii) = theta;
    
    
end

figure(2)
subplot(6,1,1);
title('Semi Major Axis');
xlabel('Julian Date');
ylabel('a [km]');
hold on
plot((t/86400 + 2451545.0), aT)

subplot(6,1,2);
title('Eccentricity');
xlabel('Julian Date');
ylabel('e ');
hold on
plot((t/86400 + 2451545.0), eT)

subplot(6,1,3);
title('Inclination');
xlabel('Julian Date');
ylabel('i');
hold on
plot((t/86400 + 2451545.0), eT)







