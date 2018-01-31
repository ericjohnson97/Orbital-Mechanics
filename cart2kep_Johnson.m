function [a,e,i,OMEGA,omega,theta] = cart2kep_Johnson(rVec,vVec,mu)
k = [0,0,1];
hVec = cross(rVec,vVec);
nVec = cross(k,hVec);
eVec = (1/mu)*(cross(vVec,hVec) - ((mu*rVec)/norm(rVec)) );
e = norm(eVec);
r = norm(rVec);
v = norm(vVec);
h = norm(hVec);
n = norm(nVec);
energy = ((v^2)/(2)) - ((mu)/(r));
if( e <= 1)
   a = -mu/(2*energy);
end
if( e>1 )
   a = inf;
end

theta = acos( dot(eVec,rVec)/(e*r) );
if (dot(rVec, vVec) < 0)
   theta = 2*pi - theta;
end
i = acos( hVec(3)/h);
OMEGA = acos(nVec(1)/n);
omega = acos(dot(nVec,eVec)/(n*e));


end