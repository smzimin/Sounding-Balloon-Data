clear;
close all;

set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arial Cyr');

highestLevel = 2000;
l = 0.00001;
h = linspace(0, highestLevel, 1000);

const = 1;
z0 = 300;
K = 4;
K_high = 4;

if const == 1
    [u_exact, v_exact, u1, v1, izl] = partKoefFunction(10, z0, K, K_high, h, l);
else
    [u_exact, v_exact, u1, v1, izl] = konst_linear_k(10, z0, K, h, l);
end


% figure(777);
% plot(u_exact, h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on');
% u = spline(h, u_exact);
% v = spline(h, v_exact);
n = 7;
u = polyfit(h, u_exact, n);
v = polyfit(h, v_exact, n);
plot(polyval(u,h), polyval(v,h), 'r:', 'LineWidth', 3.5);   hold on;
plot(u_exact, v_exact, 'LineWidth', 1.5, 'LineSmoothing', 'on');


% plot(ppval(u,1:2000), 1:2000, 'ro'); hold on;


h = linspace(0, highestLevel, 1000);

step = 10^-1;
smallk = K;
bigk = K;

% k_exact = find_k_ver2( u, v, h, l, smallk, bigk, step, 0 );
 k_exact = find_k_ver4_polynomial( u, v, h, l);
%k_exact = find_k_ver4( u, v, h, l);
%k_exact = find_k_ver5(alpha, beta, u, v, h, l, smallk, bigk, eps, 0 );
%k_exact = find_k_ver2_polynomial(0, 1, u, v, h, l, smallk, bigk, step, 0 );

du1 = polyder(u);
       du2 = polyder(du1);
       du3 = polyder(du2);
       
       dv1 = polyder(v);
       dv2 = polyder(dv1);
       dv3 = polyder(dv2);

k_const = l*sum(polyval(u,h).*polyval(dv2,h)-polyval(v,h).*polyval(du2,h)) / ...
            sum(polyval(dv2,h).^2+polyval(du2,h).^2)
        
       plot(k_const*ones(size(h)), h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on');
        
        I1 = (polyval(du2,h).*k_exact(1,:) + polyval(du1,h).*k_exact(2,:) + l*polyval(v,h)).^2 + ...
            (polyval(dv2,h).*k_exact(1,:) + polyval(dv1,h).*k_exact(2,:) - l*polyval(u,h)).^2;
        I2_const = (polyval(du2,h)*k_const + l*polyval(v,h)).^2 + ...
    (polyval(dv2,h)*k_const - l*polyval(u,h)).^2 ;

       integral1 = trapz(h,I1)
       integral2_const = trapz(h,I2_const)
        
        fr = integral1 / integral2_const

figure(2);
if const == 1
    plot(K*ones(size(h(h <= z0))), h(h <= z0),'LineSmoothing', 'on');
    plot(K_high*ones(size(h(h >= z0))), h(h >= z0),'LineSmoothing', 'on');
else
    plot(K*ones(size(h(h <= z0))), h(h <= z0),'LineSmoothing', 'on');
    plot(h(h >= z0) - (z0 - K), h(h >= z0),'LineSmoothing', 'on');
end


n = 4;
h = h';
k = polyfit(h, k_exact(1,:)', n);
dk = polyfit(h, k_exact(2,:)', n);
% plot(k_exact(1,:)', h, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
% plot(polyval(k,h(100:end)), h(100:end), 'LineWidth', 2.5, 'LineSmoothing', 'on'); 
% k = spline(h, k_exact(1,:));
% dk = spline(h, k_exact(2,:));

%veter0 = [ppval(u,h(1)) ppval(fnder(u, 1),h(1)) ppval(v,h(1)) ppval(fnder(v, 1),h(1))];
veter0 = [polyval(u,h(1)) polyval(polyder(u),h(1)) polyval(v,h(1)) polyval(polyder(v),h(1))];
[u_vosst, v_vosst] = vosstanovitelnaya_function( k, dk, h, l, veter0, k_exact);


figure(3);
plot(u_exact, v_exact, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');  hold on;
plot(u_vosst, v_vosst, 'LineWidth', 1.2, 'LineSmoothing', 'on');