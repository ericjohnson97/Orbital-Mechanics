function asrp = accelSRP(gamma, Cr, am, rSun)
    rhoSRP = 4.57E-6;
    asrp = -rhoSRP*gamma*Cr*am*(rSun/norm(rSun));
end