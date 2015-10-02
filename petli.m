clear
close all

%% генерируем петлю
r = 10;
a = linspace(0, 1.2*pi, 20);
u_exact = r*sin(a);
v_exact = r*a.*(2*pi - a).*cos(a);
p = linspace(0, 1, length(a));
l = 10;



figure(1);
%subplot(121)


% %% приближаем ее многочленами
if length(p) <= 9
    n = 3;
elseif length(p) <= 21
    n = round(length(p)/3);
else
    n = 7;
end
u = polyfit(p, u_exact, n);
v = polyfit(p, v_exact, n);
h = linspace(p(1), p(end), 1000);
h = h';
q = h;

plot(polyval(u,h), polyval(v,h), 'LineWidth', 2.5, 'LineSmoothing', 'on');  hold on;
plot(u_exact, v_exact, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(u_exact(1), v_exact(1), 'g:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');


%% задаем параметры и восстанавливаем k(z)
eps = 10^-10;
smallk = 1e-10; % 3.5559e-06
bigk = 1000000;
alpha = 0*20*10^-3;
beta = 2;


%[k_exact, wrong] = find_k_ver5(alpha, beta, u, v, h, l, smallk, bigk, eps, 0 );
k_exact = find_k_ver4_polynomial( u, v, h, l);
%k_exact = find_k_ver2_polynomial(alpha, beta, u, v, h, l, smallk, bigk, eps, 1 );
% 
% % n = 2;
% % k = polyfit(h, k_exact(1,:)', n);
% % dk = polyfit(h, k_exact(2,:)', n);
% 
% % h = h(100:900);
% % veter0 = [polyval(u,h(1)) polyval(polyder(u),h(1)) polyval(v,h(1)) polyval(polyder(v),h(1))];
% % [u_vosst, v_vosst] = vosstanovitelnaya_function( k, dk, h, l, veter0);
% % 
% % 
% % figure(3);
% % plot(polyval(u, h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');  hold on;
% % plot(u_vosst, h, 'LineWidth', 1.2, 'LineSmoothing', 'on');
% % 
% % figure(4);
% % plot(polyval(v, h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');  hold on;
% % plot(v_vosst, h, 'LineWidth', 1.2, 'LineSmoothing', 'on');
% % 
% % 
% % figure(5);
% % plot(polyval(polyder(u),q), polyval(polyder(v),q), 'LineWidth', 2.5, 'LineSmoothing', 'on');
% %figure; plot(atan(polyval(v,h)./polyval(u,h)), h)
% 
% du1 = polyder(u);
% du2 = polyder(du1);
% du3 = polyder(du2);
% dv1 = polyder(v);
% dv2 = polyder(dv1);
% dv3 = polyder(dv2);
% 
%  
%  dk2 = polyval(du1, h).^2 + polyval(dv1, h).^2;
%  dk1 = polyval(du1, h).*polyval(du2, h) + polyval(dv1, h).*polyval(dv2, h);
%  dk0 = polyval(du1, h).*polyval(du3, h) + polyval(dv1, h).*polyval(dv3, h);
%  
%  figure;
%  plot(dk2, h,'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%  plot(dk2(1), h(1), 'ro', 'LineWidth', 2.5, 'LineSmoothing', 'on');
% %  figure;
% %  plot(dk1, h,'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
% %   plot(dk1(1), h(1), 'ro', 'LineWidth', 2.5, 'LineSmoothing', 'on');
% %  figure;
% %  plot(dk0, h,'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
% %   plot(dk0(1), h(1), 'ro', 'LineWidth', 2.5, 'LineSmoothing', 'on');