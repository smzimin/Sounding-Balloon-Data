function k = find_k_ver2( u, v, h, l, smallk, bigk, step, print_errors )

opt = odeset('AbsTol', 1e-7, 'RelTol', 1e-7, 'MaxStep', (h(end)-h(1))/1000);

du1 = fnder(u, 1);
du2 = fnder(u, 2);
du3 = fnder(u, 3);

dv1 = fnder(v, 1);
dv2 = fnder(v, 2);
dv3 = fnder(v, 3);


ua = ppval(u, h(1));
d1ua = ppval(du1, h(1));
d2ua = ppval(du2, h(1));

va = ppval(v, h(1));
d1va = ppval(dv1, h(1));
d2va = ppval(dv2, h(1));


ub = ppval(u, h(end));
d1ub = ppval(du1, h(end));
d2ub = ppval(du2, h(end));

vb = ppval(v, h(end));
d1vb = ppval(dv1, h(end));
d2vb = ppval(dv2, h(end));


eps = 10^-3;

close all
figure(1);

k0 = smallk:step:bigk;
errors = size(k0);
length(k0)

k_opt = -999;
for i = 1:length(k0)
    i
    dk0 = - ( (d2ua*d1ua+d2va*d1va)*k0(i) - l*(ua*d1va-va*d1ua) ) / (d1ua^2+d1va^2);
    [~, k] = ode45(@(t, y) odefun(t, y), h, [k0(i) dk0], opt);
    k = k';
    
    plot(k(1,h<1000), h(h<1000), 'LineSmoothing', 'on'); hold on;
    errors(i) = (d1ub^2+d1vb^2)*k(2,end) + (d2ub*d1ub+d2vb*d1vb)*k(1,end) - l*(ub*d1vb-vb*d1ub);
    if abs(errors(i)) < eps
        k_opt = k0(i)
        %errors(i)
        plot(k(1,h<1000), h(h<1000), 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on'); hold on;
    end
end

if k_opt ~= -999
dk0 = - ( (d2ua*d1ua+d2va*d1va)*k_opt - l*(ua*d1va-va*d1ua) ) / (d1ua^2+d1va^2);
[~, k] = ode45(@(t, y) odefun(t, y), h, [k_opt dk0], opt);
k = k';
end

xlabel('k(z)','FontName','Arial Cyr','FontSize',20);
ylabel('Высота, метры','FontName','Arial Cyr','FontSize',20);
title('Нужный k(z) - красный','FontName','Arial Cyr','FontSize',20);

if print_errors == 1
    figure(2);
    plot(k0, errors, 'LineWidth', 2.5, 'LineSmoothing', 'on');
    xlabel('k(0)','FontName','Arial Cyr','FontSize',20);
    title('Невязка на правом конце','FontName','Arial Cyr','FontSize',20);
end

function dydx = odefun(x,y)
    
d1uz = ppval(du1, x);
d2uz = ppval(du2, x);
d3uz = ppval(du3, x);

d1vz = ppval(dv1, x);
d2vz = ppval(dv2, x);
d3vz = ppval(dv3, x);

dydx = [y(2);
    -(2*(d2uz*d1uz+d2vz*d1vz)*y(2) + (d3uz*d1uz+d3vz*d1vz)*y(1)  - 1e-7 / y(1)^8 ) / (d1uz^2+d1vz^2)];

end

end