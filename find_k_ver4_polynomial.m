function k = find_k_ver4_polynomial( u, v, h, l)

opt = odeset('AbsTol', 1e-12, 'RelTol', 1e-12, 'MaxStep', (h(end)-h(1))/1000);

du1 = polyder(u);
du2 = polyder(du1);
du3 = polyder(du2);
dv1 = polyder(v);
dv2 = polyder(dv1);
dv3 = polyder(dv2);

ua = polyval(u,h(1));
d1ua = polyval(du1,h(1));
d2ua = polyval(du2,h(1));

va = polyval(v,h(1));
d1va = polyval(dv1,h(1));
d2va = polyval(dv2,h(1));


ub = polyval(u,h(end));
d1ub = polyval(du1,h(end));
d2ub = polyval(du2,h(end));

vb = polyval(v,h(end));
d1vb = polyval(dv1,h(end));
d2vb = polyval(dv2,h(end));


[~, k1] = ode15s(@(t, y) odefun(t, y), h, [0 1], opt);
k1 = k1';

[~, k2] = ode15s(@(t, y) odefun(t, y), h, [1 0], opt);
k2 = k2';

a11 = (d1ua^2+d1va^2)*k1(2,1) + (d2ua*d1ua+d2va*d1va)*k1(1,1);
a12 = (d1ua^2+d1va^2)*k2(2,1) + (d2ua*d1ua+d2va*d1va)*k2(1,1);
a21 = (d1ub^2+d1vb^2)*k1(2,end) + (d2ub*d1ub+d2vb*d1vb)*k1(1,end);
a22 = (d1ub^2+d1vb^2)*k2(2,end) + (d2ub*d1ub+d2vb*d1vb)*k2(1,end);
b1 = l*(ua*d1va-va*d1ua);
b2 = l*(ub*d1vb-vb*d1ub);


A = [a11 a12; a21 a22];
b = [b1; b2];

C = A \ b;

k = C(1) * k1 + C(2) * k2;

nev0 = (d1ua^2+d1va^2)*k(2,1) + (d2ua*d1ua+d2va*d1va)*k(1,1) - l*(ua*d1va-va*d1ua);
nev2000 = (d1ub^2+d1vb^2)*k(2,end) + (d2ub*d1ub+d2vb*d1vb)*k(1,end) - l*(ub*d1vb-vb*d1ub);

% figure(2);
subplot(122)
plot(k(1,:), h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on'); hold on;
xlabel('k(z)','FontName','Arial Cyr','FontSize',20);
ylabel('Высота, метры','FontName','Arial Cyr','FontSize',20);
title('Восстановленный k(z)','FontName','Arial Cyr','FontSize',20);


function dydx = odefun(x,y)
    
d1uz = polyval(du1,x);
d2uz = polyval(du2,x);
d3uz = polyval(du3,x);

d1vz = polyval(dv1,x);
d2vz = polyval(dv2,x);
d3vz = polyval(dv3,x);

dydx = [y(2);
    -(2*(d2uz*d1uz+d2vz*d1vz)*y(2) + (d3uz*d1uz+d3vz*d1vz)*y(1) ) / (d1uz^2+d1vz^2)];

end

end