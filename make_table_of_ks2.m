function [ k_exact, fr, shir ] = make_table_of_ks2( alpha, beta, left, right, smallk, bigk, eps, heightDivisor, num )
highestLevel = 2000;
numOfLevels = highestLevel/20;
t = 35:16:4800;

table_of_ks = [];
for mon = 3%3:12
for day = 9%1:31
   % [mon day]
    if mon < 10
        fid=fopen(['data\BUFRG00F.20140' num2str(mon*100+day) '.00']);
    else
        fid=fopen(['data\BUFRG00F.2014' num2str(mon*100+day) '.00']);
    end
    if fid <= 0
        continue;
    end
    dat1=fread(fid, [4800 750], 'float');
    fclose(fid);
    ncol = find (dat1(1,:)==222222);
    dat1 = dat1(:,1:ncol);
    
    
    counter = 0;
table1 = -ones(length(dat1(1,:)), numOfLevels);

for i = 1:length(dat1(1,:))
   if dat1(5,i) >= left*10^6   &&  dat1(5,i) <= right*10^6;
       
       p = dat1(t,i);
       a = pi/180*dat1(t+3,i);
       speed = 0.1*dat1(t+4,i);
       sel = p > -1000 & dat1(t+4,i) > 20 & p < highestLevel*1.5;
       if sum (sel) < 8 || min (p(sel)) > 250 || max(p(sel)) < highestLevel-250
           continue;
       end
       a = a(sel);
       speed = speed(sel);
       p = p(sel);       
       sel = p (2:end) > p(1:end-1);
       sel = [true; sel];     
       if sum (sel) < 8
           continue;
       end  
       a = a(sel);
       speed = speed(sel);
       p = p(sel) / heightDivisor;
       
       a = a(p < (highestLevel/heightDivisor));
       speed = speed(p < (highestLevel/heightDivisor));
       p = p(p < (highestLevel/heightDivisor));
       
       counter = counter + 1;
       if counter < num
           continue;
       end
       
       rng('default')
       %speed = speed + 0.25*rand(size(p));
%        speed = speed(1:end-3);
%        p = p(1:end-3);
%        a = a(1:end-3);
       
        if length(p) <= 9
            n = 3;
        elseif length(p) <= 21
            n = round(length(p)/3);
        else
            n = 7;
        end
       u = polyfit(p, speed .* sin(a), n);
       v = polyfit(p, speed .* cos(a), n);
       
       
       highestLevel = 2000 / heightDivisor;
       h = linspace(p(1), p(end), 1000);
       h = h';
       
       
       
       
        figure();
        subplot(121)
        plot(polyval(u,h), polyval(v,h), 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
        plot(speed .* sin(a), speed .* cos(a), 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
        plot(polyval(u,h(1)), polyval(v,h(1)), 'g:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
        plot(polyval(u,h(end)), polyval(v,h(end)), 'y:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
        xlabel('u','FontName','Arial Cyr','FontSize',20);
        ylabel('v','FontName','Arial Cyr','FontSize',20);
        title('Годограф','FontName','Arial Cyr','FontSize',20);

       du1 = polyder(u);
       du2 = polyder(du1);
       du3 = polyder(du2);
       
       dv1 = polyder(v);
       dv2 = polyder(dv1);
       dv3 = polyder(dv2);
       
%         figure(6);
%         plot(polyval(du1,h).^2 + polyval(dv1,h).^2, h, 'r', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%         min(polyval(du1,h).^2 + polyval(dv1,h).^2)
%         mean(polyval(du1,h).^2 + polyval(dv1,h).^2)
       
%        plot(polyval(v,h), h, 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%        plot(polyval(dv1,h), h, 'r', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%        plot(polyval(dv2,h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%        plot(polyval(dv3,h), h, 'y', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%        plot(speed .* cos(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');
       
        l = 2*7.2921*(10^-5)*sin(pi/180*0.00001*dat1(5,i));
        shir = 0.00001*dat1(5,i);
%         l = sin(pi/180*0.00001*dat1(5,i));
        %l = l * 10^6;
        %[k_exact, wrong] = find_k_ver5(alpha, beta, u, v, h, l, smallk, bigk, eps, 0 );
        k_exact = find_k_ver4_polynomial( u, v, h, l);
        %  k_exact = find_k_ver2_polynomial(alpha, beta, u, v, h, l, smallk, bigk, eps, 1 );
        
%         k_const = l*sum(polyval(u,p).*polyval(dv2,p)-polyval(v,p).*polyval(du2,p)) / ...
%             sum(polyval(dv2,p).^2+polyval(du2,p).^2)
        
        k_const = l*sum(polyval(u,h).*polyval(dv2,h)-polyval(v,h).*polyval(du2,h)) / ...
            sum(polyval(dv2,h).^2+polyval(du2,h).^2)
        
       plot(k_const*ones(size(h)), h, 'r', 'LineWidth', 1.5, 'LineSmoothing', 'on');
        
        I1 = (polyval(du2,h).*k_exact(1,:)' + polyval(du1,h).*k_exact(2,:)' + l*polyval(v,h)).^2 + ...
            (polyval(dv2,h).*k_exact(1,:)' + polyval(dv1,h).*k_exact(2,:)' - l*polyval(u,h)).^2 ;
        
        I2_const = (polyval(du2,h)*k_const + l*polyval(v,h)).^2 + ...
    (polyval(dv2,h)*k_const - l*polyval(u,h)).^2 ;
        
      

       integral1 = trapz(h,I1);
       integral2_const = trapz(h,I2_const);
        
        fr = integral1 / integral2_const

%   k_t = -1:0.1:1;
%         error = zeros(size(k_t));
%         for q = 1:length(k_t)
%             
%         I2_const = (polyval(du2,h)*k_t(q) + l*polyval(v,h)).^2 + ...
%     (polyval(dv2,h)*k_t(q) - l*polyval(u,h)).^2 ;
%         error(q) = trapz(h,I2_const);
% 
%         end
%         figure(777)
%         error
%         plot(k_t, error);
        
%         q = [ones(size(k_exact(1,:)')) h];
%         bk = regress(k_exact(1,:)', [ones(size(k_exact(1,:)')) h]);
%         bdk = regress(k_exact(2,:)', [ones(size(k_exact(1,:)')) h]);
%         
%         bk = regress([k_exact(1,1); k_exact(1,end)], [1 h(1); 1 h(end)]);
%         bdk = regress([k_exact(2,1); k_exact(2,end)], [1 h(1); 1 h(end)]);
%         
%         I2 = (du2(h).*(q * bk) + du1(h).*(q * bdk) + l*v(h)).^2 + ...
%     (dv2(h).*(q * bk) + dv1(h).*(q * bdk) - l*u(h)).^2 ;
% 
%         integral2 = trapz(h,I2)
%         
        
%         k = spline(h, k_exact(1,:));
%         dk = spline(h, k_exact(2,:));
%         
%         veter0 = [polyval(u,h(1)) polyval(polyder(u),h(1)) polyval(v,h(1)) polyval(polyder(v),h(1))];
%         [u_vosst, v_vosst] = vosstanovitelnaya_function( k, dk, h, l, veter0);
%         
%         
%         figure(3);
%         plot(speed .* sin(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%         plot(polyval(u, h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%         plot(u_vosst, h, 'LineWidth', 2.5, 'LineSmoothing', 'on');
%         figure(4);
%         plot(speed .* cos(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%         plot(polyval(v, h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%         plot(v_vosst, h, 'LineWidth', 2.5, 'LineSmoothing', 'on');
        
        
        
      %  table1 (i,:) = k(1,:);

       break
   end
end

table1 = table1(sum(table1~=-1,2)>15, :);

%save([num2str(mon*100+day) '.mat'], 'table1');


table_of_ks = [table_of_ks; table1];
end
end


end

