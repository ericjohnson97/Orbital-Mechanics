function rho = density(h)
    rho = 0;
    if(h <= 40000)
        T = 518.69 - (3.5662E-3 * h);
        rho = (6.6277E-15)*T^4.256;
    end
    if(h > 40000)
       p = 2678.4*exp((-4.8063E-5) *h);
       rho = (1.4939E-6)*p;
    end
end