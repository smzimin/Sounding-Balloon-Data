function [u, v, u1, v1, izl] = konst_linear_k(U, z0, k, z, l)
q = sqrt(l/(2*k));
        
s = sin(q*z0);
c = cos(q*z0);
sh = sinh(q*z0);
ch = cosh(q*z0);

ss = s*sh;
sc = s*ch;
cc = c*ch;
cs = c*sh;

Z = -2*sqrt(-1i*l*k);

rH = real(besselh(0,1,Z));
imH = imag(besselh(0,1,Z));

rdH = real((-1i*l)*besselh(1,1,Z)/sqrt(-1i*l*k));
imdH = imag((-1i*l)*besselh(1,1,Z)/sqrt(-1i*l*k));

A = [cs sc -rH imH;    
    sc -cs -imH -rH;
    q*(cc-ss) q*(cc+ss) -rdH imdH;
    q*(cc+ss) q*(ss-cc) -imdH -rdH];

B = [U*cc; U*ss; q*U*(cs-sc); q*U*(cs+sc)];

C = A^(-1) * B;

p = 0; izl = 0;

for i = 1:length(z)
    
    if z(i) >= z0 && p == 0
        p = 1;
        izl = i;
    end
    
    if z(i) < z0
        
        x = q*z(i);
        
        u(i) = C(1)*cos(x)*sinh(x) + C(2)*sin(x)*cosh(x) - U*cos(x)*cosh(x);
        v(i) = C(1)*sin(x)*cosh(x) - C(2)*cos(x)*sinh(x) - U*sin(x)*sinh(x);
        
         u1(1,i) = q * (C(1)*(cos(x)*cosh(x) - sin(x)*sinh(x)) + C(2)*(cos(x)*cosh(x) + sin(x)*sinh(x)) - U*(cos(x)*sinh(x) - sin(x)*cosh(x)));
         v1(1,i) = q * (C(1)*(cos(x)*cosh(x) + sin(x)*sinh(x)) + C(2)*(-cos(x)*cosh(x) + sin(x)*sinh(x)) - U*(cos(x)*sinh(x) + sin(x)*cosh(x)));
        
    else
        x = -2*sqrt(-1i*l*(k-z0+z(i)));
        

        u(i) = C(3) * real(besselh(0,1,x)) - C(4)* imag(besselh(0,1,x));
        v(i) = C(3) * imag(besselh(0,1,x)) + C(4)* real(besselh(0,1,x));
        
        u1(1,i) = real((C(3) + 1i*C(4))* (-1i*l)*besselh(1,1,x)/sqrt(-1i*l*(k-z0+z(i))) ) + U;
        v1(1,i) = imag((C(3) + 1i*C(4))* (-1i*l)*besselh(1,1,x)/sqrt(-1i*l*(k-z0+z(i))) );
    end
end


end

