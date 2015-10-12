function find_k_splines_plot_Ik(alpha, beta, u, v, h, l, a, b, eps )
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


k0 = a:eps:b;
Ik = zeros(size(k0));

for i = 1:length(k0)
    i
    dk0 = ( l*(ua*d1va-va*d1ua) - (d2ua*d1ua+d2va*d1va)*k0(i) ) /...
        (d1ua^2+d1va^2);
    [~, k_exact] = ode45(@(t, y) odefun(t, y), h, [k0(i) dk0]);
    k_exact = k_exact';
    
    I1 = (ppval(du2,h).*k_exact(1,:)' + ppval(du1,h).*k_exact(2,:)' ...
    + l*ppval(v,h)).^2 + (ppval(dv2,h).*k_exact(1,:)' ...
    + ppval(dv1,h).*k_exact(2,:)' - l*ppval(u,h)).^2 ;

    Ik(i) = trapz(h,I1);
end

figure();
plot(k0, Ik, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on'); hold on;

xlabel('k(0)','FontName','Arial Cyr','FontSize',20);
ylabel('I[k]','FontName','Arial Cyr','FontSize',20);
title('Значение функционала I[k]','FontName','Arial Cyr','FontSize',20);


function dydx = odefun(x,y)
    
d1uz = ppval(du1, x);
d2uz = ppval(du2, x);
d3uz = ppval(du3, x);

d1vz = ppval(dv1, x);
d2vz = ppval(dv2, x);
d3vz = ppval(dv3, x);

dydx = [y(2);
    -(2*(d2uz*d1uz+d2vz*d1vz)*y(2) + (d3uz*d1uz+d3vz*d1vz)*y(1) - ...
    alpha*beta / (y(1)^(beta+1))) / (d1uz^2+d1vz^2)];

end

end