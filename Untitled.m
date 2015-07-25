clear;
close all;

set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arial Cyr');

highestLevel = 2*10^3;
l = -2*7.2921*(10^-5)*sin(0.33*pi);

numOfLevels = highestLevel/50;
step = highestLevel / numOfLevels;
h = step*(0:numOfLevels);

const = 1;
z0 = 600;
K = 1000;
K_high = 300;

if const == 1
    [u_exact, v_exact, u1, v1, izl] = partKoefFunction(10, z0, K, K_high, h, l);
else
    [u_exact, v_exact, u1, v1, izl] = konst_linear_k(10, z0, K, h, l);
end

%figure(777);
%plot(u_exact, h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on');
u = spline(h, u_exact);
v = spline(h, v_exact);

% plot(ppval(u,1:2000), 1:2000, 'ro'); hold on;

              
numOfLevels = highestLevel/10;
step = highestLevel / numOfLevels;
h = step*(0:numOfLevels);

step = 10^-1;
smallk = K;
bigk = K;

k = find_k_ver2( u, v, h, l, smallk, bigk, step, 0 );
%k = find_k_ver4( u, v, h, l);


I = (ppval(fnder(u,2),h).*k(1,:) + ppval(fnder(u,1),h).*k(2,:) + l*ppval(v,h)).^2 + ...
    (ppval(fnder(v,2),h).*k(1,:) + ppval(fnder(v,1),h).*k(2,:) - l*ppval(u,h)).^2 ;


integral = trapz(h,I)

figure(1);
if const == 1
    plot(K*ones(size(h(h <= z0))), h(h <= z0),'LineSmoothing', 'on');
    plot(K_high*ones(size(h(h >= z0))), h(h >= z0),'LineSmoothing', 'on');
else
    plot(K*ones(size(h(h <= z0))), h(h <= z0),'LineSmoothing', 'on');
    plot(h(h >= z0) - (z0 - K), h(h >= z0),'LineSmoothing', 'on');
end