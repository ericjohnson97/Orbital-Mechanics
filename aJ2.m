function a = aJ2(mu, J2, Re, r)

    a(1) = ((3*mu*J2*Re^2)/(2*(norm(r)^5)))*r(1)*(5*(r(3)/norm(r))^2 - 1);
    a(2) = ((3*mu*J2*Re^2)/(2*(norm(r)^5)))*r(2)*(5*(r(3)/norm(r))^2 - 1);
    a(3) = ((3*mu*J2*Re^2)/(2*(norm(r)^5)))*r(3)*(5*(r(3)/norm(r))^2 - 3);

end