%Orbital Mechanics Homework 1
clc; clear all;

%Problem 2b:
r = [.4717, 1.2960, 1.2418];
v = [-.1601, -.43987, .59881];

T = v/norm(v);
W = cross(v,r)/(norm(cross(v,r)));
N = cross(T,W);
fprintf("The trnasformation matrix from NTW to ijk is \n");
Tntw2ijk = [N; T; W] 


%Problem 3
clear all;

r1 = [-.58309, -1.60202, 1.33197];
v1 = [-.1471, -.404165, -.489338];

r2 = [-.59154, -1.62524, 1.30331];
v2 = [-.14374, -.39492, -.49688];

r21ijk = r2 - r1;

R = r1/norm(r1);
W = cross(r1,v1)/(norm(cross(r1,v1)));
S = cross(W,R);

Tijk2RST = [R;S;W];

r21rsw = Tijk2RST*r21ijk';

fprintf("Sat 1 is %g R %g S %g W [DU] away from sat 2 in Sat 1's RSW frame \n", r21rsw);




