function [ k_exact, fr, shir ] = make_table_of_ks_splines( regularization, ...
    shir_vis, fin_par, data, plot_Ik )
k_exact = []; fr = 0; shir = 0;

%% params
alpha = regularization(1);
beta = regularization(2);

left = shir_vis(1);
right = shir_vis(2);
highestLevel = shir_vis(3);
heightDivisor = shir_vis(4);

smallk = fin_par(1);
bigk = fin_par(2);
eps = fin_par(3);

mon = data(1);
day = data(2);
num = data(3);

t = 35:16:4800;

%% read data
if mon < 10
    fid=fopen(['data\BUFRG00F.20140' num2str(mon*100+day) '.00']);
else
    fid=fopen(['data\BUFRG00F.2014' num2str(mon*100+day) '.00']);
end

dat1 = fread(fid, [4800 750], 'float');
fclose(fid);
ncol = find (dat1(1,:)==222222);
dat1 = dat1(:,1:ncol);
    
%% loop over the data
counter = 0;
for i = 1:length(dat1(1,:))
   if dat1(5,i) >= left*10^6   &&  dat1(5,i) <= right*10^6;
       
       %% extract good data
       [correct, a, speed, p, h] = get_the_data(dat1, highestLevel, ...
                                                    heightDivisor, t, i);
        if correct == 0
            continue
        end
       
       %% get observation you need (kinda stupid)
       counter = counter + 1;
       if counter < num
           continue;
       end
       
       %% interpolate
       u = spline(p, speed .* sin(a));
       v = spline(p, speed .* cos(a));
       
       %% get derivatives
       % [du1, du2, du3, dv1, dv2, dv3] = get_derivatives(u,v);
       
       %% some plots
       %plot_godograph(u, v, h, speed, a)
       plot_uv(u, v, h, speed, a, p)
       %plot_derivatives(u, v, h, speed, a, p)
       %plot(ppval(du1,h).^2 + ppval(dv1,h).^2, h, 'r', 'LineWidth', 2.5, 'LineSmoothing', 'on');
       
       %% main calculations
       shir = 0.00001*dat1(5,i);
       l = 2*7.2921*(10^-5)*sin(pi/180*shir);
       
%        if plot_Ik == 0
%            if alpha == 0
%                k_exact = find_k_ver4( u, v, h, l);
%            else
%                k_exact = find_k_ver5_splines(alpha, beta, u, v, h, ...
%                    l*10^8, smallk, bigk, eps, 0 )*10^-8;
%            end
%        else
%            find_k_splines_plot_Ik(alpha, beta, u, v, h, l, smallk,...
%                bigk, eps )
%        end
%        %% find k_const and compare
%        [k_const, I1, I2] = k_const_and_integrals(u, v, h, k_exact, l);
%        fr = I1 / I2

       %% find uv using k(z)
       %[u_vosst, v_vosst] = uv_twice(u, v, l, k_exact, h, speed, a, p);

       break
   end
end

end

function [correct, a, speed, p, h] = get_the_data(dat1, highestLevel, ...
    heightDivisor, t, i)
correct = 1;
a = []; p = []; h = [];
speed = [];

p = dat1(t,i);
a = pi/180*dat1(t+3,i);
speed = 0.1*dat1(t+4,i);
sel = p > -1000 & dat1(t+4,i) > 20 & p < highestLevel*1.5;
if sum (sel) < 8 || min (p(sel)) > 250 || max(p(sel)) < highestLevel-250
   correct = 0;
end

if correct == 1
    a = a(sel);
    speed = speed(sel);
    p = p(sel);       
    sel = p (2:end) > p(1:end-1);
    sel = [true; sel];     
    if sum (sel) < 8
       correct = 0;
    end
    
    if correct == 1
        a = a(sel);
        speed = speed(sel);
        p = p(sel) / heightDivisor;

        a = a(p < (highestLevel/heightDivisor));
        speed = speed(p < (highestLevel/heightDivisor));
        p = p(p < (highestLevel/heightDivisor));

        h = linspace(p(1), p(end), 1000);
        h = h';
    end
end
end

function plot_godograph(u, v, h, speed, a)
figure();
subplot(121)
plot(ppval(u,h), ppval(v,h), 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
plot(speed .* sin(a), speed .* cos(a), 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(ppval(u,h(1)), ppval(v,h(1)), 'g:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(ppval(u,h(end)), ppval(v,h(end)), 'y:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
xlabel('u','FontName','Arial Cyr','FontSize',20);
ylabel('v','FontName','Arial Cyr','FontSize',20);
title('Годограф','FontName','Arial Cyr','FontSize',20);
end

function plot_derivatives(u, v, h, speed, a, p)
[du1, du2, du3, dv1, dv2, dv3] = get_derivatives(u,v);
figure(6);
plot(ppval(v,h), h, 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
plot(ppval(dv1,h), h, 'r', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(ppval(dv2,h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(ppval(dv3,h), h, 'y', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(speed .* cos(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
end

function plot_uv(u, v, h, speed, a, p)
figure(6);
plot(ppval(u,h), h, 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
plot(speed .* sin(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
figure(7);
plot(ppval(v,h), h, 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
plot(speed .* cos(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
end

function [du1, du2, du3, dv1, dv2, dv3] = get_derivatives(u,v)
du1 = fnder(u, 1);
du2 = fnder(u, 2);
du3 = fnder(u, 3);

dv1 = fnder(v, 1);
dv2 = fnder(v, 2);
dv3 = fnder(v, 3);
end

function [k_const, I1, I2] = k_const_and_integrals(u, v, h, k_exact, l)
[du1, du2, ~, dv1, dv2, ~] = get_derivatives(u,v);
k_const = l*sum(ppval(u,h).*ppval(dv2,h)-ppval(v,h).*ppval(du2,h)) / ...
            sum(ppval(dv2,h).^2+ppval(du2,h).^2);
    
%figure(1);
plot(k_const*ones(size(h)), h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on');
        
I1 = (ppval(du2,h).*k_exact(1,:)' + ppval(du1,h).*k_exact(2,:)' ...
    + l*ppval(v,h)).^2 + (ppval(dv2,h).*k_exact(1,:)' ...
    + ppval(dv1,h).*k_exact(2,:)' - l*ppval(u,h)).^2 ;

I2_const = (ppval(du2,h)*k_const + l*ppval(v,h)).^2 + ...
(ppval(dv2,h)*k_const - l*ppval(u,h)).^2;
        
I1 = trapz(h,I1);
I2 = trapz(h,I2_const);
end

function [u_vosst, v_vosst] = uv_twice(u, v, l, k_exact, h, speed, a, p)
k = spline(h, k_exact(1,:));
dk = spline(h, k_exact(2,:));

veter0 = [ppval(u,h(1)) ppval(fnder(u,1),h(1)) ppval(v,h(1)) ppval(fnder(v,1),h(1))];
[u_vosst, v_vosst] = vosstanovitelnaya_function( k, dk, h, l, veter0);

figure(3);
plot(speed .* sin(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
plot(ppval(u, h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(u_vosst, h, 'LineWidth', 2.5, 'LineSmoothing', 'on');
figure(4);
plot(speed .* cos(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
plot(ppval(v, h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
plot(v_vosst, h, 'LineWidth', 2.5, 'LineSmoothing', 'on');
end

function [u_vosst, v_vosst] = vosstanovitelnaya_function( k, dk, ...
    h, l, veter0)

opt = odeset('AbsTol', 1e-10, 'RelTol', 1e-10, 'MaxStep', (h(end)-h(1))/1000);

%[~, y] = ode15s(@(t, y) odefun(t, y), h, veter0, opt);
[~, y] = ode45(@(t, y) odefun(t, y), h, veter0, opt);

u_vosst = y(:,1)';
v_vosst = y(:,3)';

function dydx = odefun(x,y)
    
kz = ppval(k, x); 
dkz = ppval(dk, x);

dydx = [y(2);
    -(y(2)*dkz + l*y(3)  ) / kz;
    y(4);
    (-y(4)*dkz + l*y(1)  ) / kz];

end

end