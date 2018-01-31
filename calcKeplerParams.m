function [h1,e,eMag,a] = calcKeplerParams(r,v,u) 

h1 =  cross( r,v );
e = ( cross( v,h1 ) / u ) - ( r / norm(r) );
eMag = norm(e);
a = ( norm(h1)^2 ) / ( u*( 1 - eMag^2 ) );

end
