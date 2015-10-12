function [k, wrong] = find_k_ver5(alpha, beta, u, v, h, l, a, b, eps, print_errors )
wrong = 0;
k = [];

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

% aset = -[-10.^(10:-1:-10) 10.^(-10:10)];
% fa = f(aset(1));
% for i = 2:length(aset)
%     fa(i) = f(aset(i));
%     if fa(i-1)*fa(i)<0
%         b = aset(i-1);
%         a = aset(i);
%         break;
%     end
% end
% 
% f_a = fa(i-1);
% f_b = fa(i);

f_a = f(a);
f_b = f(b);

while true
    if (f_a * f_b > 0 && f_a > 0) || b >= 10^100 % || f_a > f_b
        wrong = 1;
        [f_a  f_b]
        break
    elseif f_a * f_b > 0 && f_a < 0
        b = 10*b
        f_b = f(b);
    else
        c = (a + b) / 2;
        break
    end
end



if wrong == 0
    errors = [];
    numOfIter = 0;
    
    f_c = 10;

    while (b - a) / 2 >= eps
        f_c = f(c);
        if f_c*f_b <= 0
            a = c;
        else
            b = c;
        end
        c = (a + b) / 2;
        numOfIter = numOfIter + 1;
        errors = [errors f_c];
       % [numOfIter a b f_c]
    end
    K0_found = c;

    dK0 = ( l*(ua*d1va-va*d1ua) - (d2ua*d1ua+d2va*d1va)*K0_found  ) / (d1ua^2+d1va^2);
    [~, k] = ode45(@(t, y) odefun(t, y), h, [K0_found dK0], opt);

    k = k';

    %figure();
    subplot(122)
    plot(k(1,:), h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on'); hold on;

    xlabel('k(z)','FontName','Arial Cyr','FontSize',20);
    ylabel('������, ���������','FontName','Arial Cyr','FontSize',20);
    title('k(z)','FontName','Arial Cyr','FontSize',20);
    f_c = f_c

    if print_errors == 1
        figure(2);
        plot(errors, 'LineWidth', 2.5, 'LineSmoothing', 'on');
        xlabel('k(0)','FontName','Arial Cyr','FontSize',20);
        title('������� �� ������ �����','FontName','Arial Cyr','FontSize',20);
    end
end

function residual = f(k0)
    
dk0 = ( l*(ua*d1va-va*d1ua) - (d2ua*d1ua+d2va*d1va)*k0  ) / (d1ua^2+d1va^2);
[~, k_temp] = ode45(@(t, y) odefun(t, y), h, [k0 dk0]);
k_temp = k_temp';

residual = (d1ub^2+d1vb^2)*k_temp(2,end) + ...
    (d2ub*d1ub+d2vb*d1vb)*k_temp(1,end) - l*(ub*d1vb-vb*d1ub);
end


function dydx = odefun(x,y)
    
d1uz = polyval(du1,x);
d2uz = polyval(du2,x);
d3uz = polyval(du3,x);

d1vz = polyval(dv1,x);
d2vz = polyval(dv2,x);
d3vz = polyval(dv3,x);

dydx = [y(2);
    -(2*(d2uz*d1uz+d2vz*d1vz)*y(2) + (d3uz*d1uz+d3vz*d1vz)*y(1) - ...
    alpha*beta / (y(1)^(beta+1))) / (d1uz^2+d1vz^2)];

end

end