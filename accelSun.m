function [asun, rsSun] = accelSun(r, t)
    muSun = 1.3271244E11;

    JDut1 = 2451545.0 + (t/86400);
    Tut1 = (JDut1 - 2451545)/36525;
    
    lambdaMSun = 280.46 + 36000.771*Tut1;
    MSun = 357.5291092 + 35999.0534*Tut1;
    
    %assume Tut1 is essentially = Ttbd
    Ttdb = Tut1;
    lambdaE = lambdaMSun + 1.914666471*sind(MSun) + 0.019994643*sind( 2*MSun);
    
    rSunMag = 1.000140612 - 0.016708617*cosd(MSun) - 0.000139589*cosd(2*MSun);
    epsilon = 23.439291 - 0.0130042*Ttdb;
    
    rSun = [ rSunMag*cosd(lambdaE), rSunMag*cosd(epsilon)*sind(lambdaE), rSunMag*sind(epsilon)*sind(lambdaE)] ;
    
    TT = Tut1;
    zeta= (2306.2181*TT + 0.30188*TT^2 + 0.017998*TT^3)*(1/3600);
    phi = (2004.3109*TT - 0.42665*TT^2 - 0.041833*TT^3)*(1/3600);
    z = (2306.2181*TT + 1.09468*TT^2 + 0.018203*TT^3)*(1/3600);
    
    R3zeta = [ cosd(zeta), -sind(zeta), 0; sind(zeta), cosd(zeta), 0; 0, 0, 1];
    R2phi = [ cosd(-phi), 0, sind(-phi); 0, 1, 0; -sind(-phi), 0, cosd(-phi)];
    R3z = [ cosd(z), -sind(z), 0; sind(z), cosd(z), 0; 0, 0, 1];
    
    P = R3zeta*R2phi*R3z;
    rSunGCRF = P*rSun'*149600000;
    rsSun(1) = -r(1) + rSunGCRF(1);
    rsSun(2) = -r(2) + rSunGCRF(2);
    rsSun(3) = -r(3) + rSunGCRF(3);

    rsSunMag = sqrt( (rsSun(1)).^2 + (rsSun(2)).^2 + (rsSun(3)).^2 );
    
    asun(1) = (muSun/(rsSunMag^3))*r(1);
    asun(2) = (muSun/(rsSunMag^3))*r(2);
    asun(3) = (muSun/(rsSunMag^3))*r(3);
end