function ap = aJ2J3(mu, Re, J2, J3, r)
    a2 = aJ2(mu, J2, Re, r);
    a3 = aJ3(mu, J3, Re, r);

    ap = a2 + a3;

end
