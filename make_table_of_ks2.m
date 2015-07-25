function [ table_of_ks, bu, bv ] = make_table_of_ks2( highestLevel, numOfLevels, t, h, left, right, smallk, bigk, eps, num )
table_of_ks = [];
for mon = 3%3:12
for day = 9%1:31
    [mon day]
    if mon < 10
        fid=fopen(['BUFRG00F.20140' num2str(mon*100+day) '.00']);
    else
        fid=fopen(['BUFRG00F.2014' num2str(mon*100+day) '.00']);
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
       p = p(sel);
       
       counter = counter + 1;
       if counter < num
           continue;
       end
       
       
       
       bu = regress(speed .* sin(a), [ones(size(p)) p p.^2 p.^3 p.^4])
       u = @(z) [ones(size(z)) z z.^2 z.^3 z.^4] * bu;
       du1 = @(z) [ones(size(z)) z z.^2 z.^3] * (bu(2:5).* [1; 2; 3; 4]);
       du2 = @(z) [ones(size(z)) z z.^2] * (bu(3:5).* [2; 6; 12]);
       du3 = @(z) [ones(size(z)) z] * (bu(4:5).* [6; 24]);
       
       bv = regress(speed .* cos(a), [ones(size(p)) p p.^2 p.^3 p.^4])
       v = @(z) [ones(size(z)) z z.^2 z.^3 z.^4] * bv;
       dv1 = @(z) [ones(size(z)) z z.^2 z.^3] * (bv(2:5).* [1; 2; 3; 4]);
       dv2 = @(z) [ones(size(z)) z z.^2] * (bv(3:5).* [2; 6; 12]);
       dv3 = @(z) [ones(size(z)) z] * (bv(4:5).* [6; 24]);
       
       
       highestLevel = 1000;
       numOfLevels = highestLevel;
       step = highestLevel / numOfLevels;
       h = step*(1:numOfLevels); 
       
       h = h';
       figure(1);
%        plot(v(h), h, 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%        plot(dv1(h), h, 'r', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%        plot(dv2(h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%        plot(dv3(h), h, 'y', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%        plot(speed .* cos(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on');


%         plot(du1(h).^2, h, 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%         plot(dv1(h).^2, h, 'r', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
    %  plot(du1(h).^2+dv1(h).^2,h');
       
        l = 2*7.2921*10^-5*sin(pi/180*0.00001*dat1(5,i));
        %k_exact = find_k_ver4_polynomial( bu, bv, h, l);
        k_exact = find_k_ver2_polynomial( bu, bv, h, l, smallk, bigk, eps, 1 );
        
        
        I = (du2(h).*k_exact(1,:)' + du1(h).*k_exact(2,:)' + l*v(h)).^2 + ...
    (dv2(h).*k_exact(1,:)' + dv1(h).*k_exact(2,:)' - l*u(h)).^2 ;

        integral = trapz(h,I)
        
        q = [ones(size(k_exact(1,:)')) h];
        bk = regress(k_exact(1,:)', [ones(size(k_exact(1,:)')) h]);
        bdk = regress(k_exact(2,:)', [ones(size(k_exact(1,:)')) h]);
        
        bk = regress([k_exact(1,1); k_exact(1,end)], [1 h(1); 1 h(end)]);
        bdk = regress([k_exact(2,1); k_exact(2,end)], [1 h(1); 1 h(end)]);
        
        I2 = (du2(h).*(q * bk) + du1(h).*(q * bdk) + l*v(h)).^2 + ...
    (dv2(h).*(q * bk) + dv1(h).*(q * bdk) - l*u(h)).^2 ;

        integral2 = trapz(h,I2)
        
%         k = spline(h, k_exact(1,:));
%         dk = spline(h, k_exact(2,:));
%         
%         [u_vosst, v_vosst] = vosstanovitelnaya_function( k, dk, h, l, [u(h(1)) du1(h(1)) v(h(1)) dv1(h(1))]);
%         
%         
%         figure(3);
%         plot(speed .* sin(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%         plot(u(h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
%         plot(u_vosst, h, 'LineWidth', 2.5, 'LineSmoothing', 'on');
%         figure(4);
%         plot(speed .* cos(a), p, 'r:o', 'LineWidth', 2.5, 'LineSmoothing', 'on'); hold on;
%         plot(v(h), h, 'g', 'LineWidth', 2.5, 'LineSmoothing', 'on');
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

