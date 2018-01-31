function [theta, r, v] = calcThetaRV(t,Me,E,eMag,a,r,v,theta,u)
%PART C

T = (2*pi) * sqrt( (a^3)/(u)  );
tp = 0 - ( (T*Me) / (2*pi) );
fprintf('C.) \n');
fprintf('the time of pariapse is %f [s]\n',tp);

r(1) = norm(r);
v(1) = norm(v);
theta(1) = theta;
for j=1:length(t)
    Me = ( (2*pi) / T ) * ( t(j) - tp );
    for i=1:100
        fEo = E - eMag*sin(E) - Me;
        fEoPrime = 1 - eMag*cos(E);
        E1 = E - (fEo/fEoPrime);
        if ( E1 - E < 10^10)
            E = E1;
            break;
        end
        E = E1;
    end
    theta(j) = acos( (cos(E) - eMag) / ( 1 - eMag*cos(E) ));
    r(j) = ( a*(1 - eMag^2) ) / ( 1 + eMag*cos(theta(j)) );
    v(j) = sqrt( (u/a) + ((2*u)/r(j)) );
end
fprintf('D.) \n');
fprintf('the true anomally at t=15,000 is %f [s]\n',theta(30));
end
