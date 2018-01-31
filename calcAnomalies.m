function [E, Me, theta] = calcAnomalies(eMag,e,r,a)
theta = acos( ( dot(e,r) ) / ( eMag * norm(r) ) );
%eccentric anomally - angle from center to spacecraft
E = acos( (1/eMag) - (norm(r)/(eMag*a) ) );
Me = E - eMag*sin(E);
fprintf('B.) \n');
fprintf('the true anomally is %f [rad]\n',theta);
fprintf('the eccentric anomally is %f [rad]\n',E);
fprintf('the mean anomally is %f [rad]\n',Me);
end
