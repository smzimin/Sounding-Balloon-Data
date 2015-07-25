function [u,v,u1,v1, izl] = partKoefFunction(U, z0, kLeft, kRight, t, l)
qLeft = sqrt(l/(2*kLeft));
qRight = sqrt(l/(2*kRight));
wl = qLeft*z0;
wr = qRight*z0;

a11 = sin(wr)*exp(-wr);
a12 = cos(wr)*exp(-wr);
a13 = -cos(wl)*sinh(wl);
a14 = -sin(wl)*cosh(wl);

a21 = cos(wr)*exp(-wr);
a22 = -sin(wr)*exp(-wr);
a23 = -sin(wl)*cosh(wl);
a24 = cos(wl)*sinh(wl);

a31 = qRight*kRight*(cos(wr)-sin(wr))*exp(-wr);
a32 = -qRight*kRight*(cos(wr)+sin(wr))*exp(-wr);
a33 = -qLeft*kLeft*(cos(wl)*cosh(wl) - sin(wl)*sinh(wl));
a34 = -qLeft*kLeft*(cos(wl)*cosh(wl) + sin(wl)*sinh(wl));

a41 = -qRight*kRight*(cos(wr)+sin(wr))*exp(-wr);
a42 = -qRight*kRight*(cos(wr)-sin(wr))*exp(-wr);
a43 = -qLeft*kLeft*(cos(wl)*cosh(wl) + sin(wl)*sinh(wl));
a44 = -qLeft*kLeft*(-cos(wl)*cosh(wl) + sin(wl)*sinh(wl));

r1 = cos(wl)*cosh(wl);
r2 = sin(wl)*sinh(wl);
r3 = qLeft*kLeft*(cos(wl)*sinh(wl) - sin(wl)*cosh(wl));
r4 = qLeft*kLeft*(cos(wl)*sinh(wl) + sin(wl)*cosh(wl));

A = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];
R = [r1; r2; r3; r4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C =  A^(-1) * (-U) *R;

u = zeros(1, length(t));
v = zeros(1, length(t));

u1 = zeros(1, length(t));
v1 = zeros(1, length(t));

p = 0; izl = 0;

for i = 1:length(t)
    
    if t(1,i) >= z0 && p == 0
        p = 1;
        izl = i;
    end
    
    if t(1,i) <= z0
        
        x = qLeft*t(1,i);
        
        u(1,i) = C(3,1)*cos(x)*sinh(x) + C(4,1)*sin(x)*cosh(x) - U*cos(x)*cosh(x);
        v(1,i) = C(3,1)*sin(x)*cosh(x) - C(4,1)*cos(x)*sinh(x) - U*sin(x)*sinh(x);
        
        u1(1,i) = qLeft * (C(3,1)*(cos(x)*cosh(x) - sin(x)*sinh(x)) + C(4,1)*(cos(x)*cosh(x) + sin(x)*sinh(x)) - U*(cos(x)*sinh(x) - sin(x)*cosh(x)));
        v1(1,i) = qLeft * (C(3,1)*(cos(x)*cosh(x) + sin(x)*sinh(x)) + C(4,1)*(-cos(x)*cosh(x) + sin(x)*sinh(x)) - U*(cos(x)*sinh(x) + sin(x)*cosh(x)));
        
    else
        x = qRight*t(1,i);
        
        u(1,i) = (C(1,1)*sin(x) + C(2,1)*cos(x))*exp(-x);
        v(1,i) = (C(1,1)*cos(x) - C(2,1)*sin(x))*exp(-x);
        
        u1(1,i) = qRight *( (C(1,1)*(cos(x) - sin(x)) - C(2,1)*(cos(x) + sin(x)))*exp(-x));
        v1(1,i) = qRight *( (-C(1,1)*(cos(x) + sin(x)) - C(2,1)*(cos(x) - sin(x)))*exp(-x));
    end
end