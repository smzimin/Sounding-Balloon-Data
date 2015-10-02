function k = find_k_ver2_polynomial(alpha, beta, u, v, h, l, smallk, bigk, step, print_errors )

opt = odeset('AbsTol', 1e-7, 'RelTol', 1e-7, 'MaxStep', (h(end)-h(1))/1000);


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

u_temp = polyval(u,h);
v_temp = polyval(v,h);
du1_temp = polyval(du1,h);
dv1_temp = polyval(dv1,h);
du2_temp = polyval(du2,h);
dv2_temp = polyval(dv2,h);


eps = 10^-3;

close all
figure(1);

k0 = smallk:step:bigk;
errors = zeros(size(k0));
max_res = zeros(size(k0));
length(k0)

k_opt = -999;
for i = 1:length(k0)
    i
    dk0 = - ( (d2ua*d1ua+d2va*d1va)*k0(i) - l*(ua*d1va-va*d1ua) ) / (d1ua^2+d1va^2);
    [~, k] = ode45(@(t, y) odefun(t, y), h, [k0(i) dk0], opt);
    k = k';
    
    if length(k(1,:)) == length(h)
    
      errors(i) = (d1ub^2+d1vb^2)*k(2,end) + (d2ub*d1ub+d2vb*d1vb)*k(1,end) - l*(ub*d1vb-vb*d1ub);
      
%       q = (du1_temp.^2+dv1_temp.^2).*k(2,:)';
%       q = q + (du2_temp.*du1_temp+dv2_temp.*dv1_temp).*k(1,:)';
%       q = q - l*(u_temp.*dv1_temp-v_temp.*du1_temp);
%       max_res(i) = max(abs(q));
      
       plot(k(1,:), h, 'LineWidth', 2, 'LineSmoothing', 'on'); hold on;
    
    if abs(errors(i)) < eps
        k_opt = k0(i)
        %errors(i)
     %   plot(k(1,:), h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on'); hold on;
    end
    end
end

if k_opt ~= -999
dk0 = - ( (d2ua*d1ua+d2va*d1va)*k_opt - l*(ua*d1va-va*d1ua) ) / (d1ua^2+d1va^2);
[~, k] = ode45(@(t, y) odefun(t, y), h, [k_opt dk0], opt);
k = k';
end

xlabel('k(z)','FontName','Arial Cyr','FontSize',20);
ylabel('Высота, километры','FontName','Arial Cyr','FontSize',20);
title('Восстановленный k(z)','FontName','Arial Cyr','FontSize',20);

if print_errors == 1
    figure(2);
    plot(k0, errors, 'LineWidth', 2.5, 'LineSmoothing', 'on');
    xlabel('k(0)','FontName','Arial Cyr','FontSize',20);
    title('Невязка на правом конце','FontName','Arial Cyr','FontSize',20);
    
%     figure(3);
%     plot(k0, max_res, 'LineWidth', 2.5, 'LineSmoothing', 'on');
%     xlabel('k(0)','FontName','Arial Cyr','FontSize',20);
%     title('Максимальная невязка','FontName','Arial Cyr','FontSize',20);
end

function dydx = odefun(x,y)
    
d1uz = polyval(du1,x);
d2uz = polyval(du2,x);
d3uz = polyval(du3,x);

d1vz = polyval(dv1,x);
d2vz = polyval(dv2,x);
d3vz = polyval(dv3,x);

dydx = [y(2);
    -(2*(d2uz*d1uz+d2vz*d1vz)*y(2) + (d3uz*d1uz+d3vz*d1vz)*y(1) - alpha / y(1)^beta ) / (d1uz^2+d1vz^2)];

end

end