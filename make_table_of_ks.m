function [ table_of_ks, u, v, ng, nw ] = make_table_of_ks( highestLevel, numOfLevels, t, h, left, right, smallk, bigk, eps, num )
table_of_ks = [];
num_of_good_k = 0;
num_of_wrong_k = 0;
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
    
    
table1 = -ones(length(dat1(1,:)), numOfLevels);

counter = 0;

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
       
       ub = regress(speed .* sin(a), [ones(size(p)) p.^2 p.^3 p.^4]);
       
              
       wrong = 0;
    
        u = spline(p, speed .* sin(a));
        v = spline(p, speed .* cos(a));
%        
%        du1 = fnder(u, 1);
% 
%         dv1 = fnder(v, 1);
%        
%        semilogx(ppval(du1,h).^2+ppval(dv1,h).^2,h);


       highestLevel = 2000;
       numOfLevels = highestLevel;
       step = highestLevel / numOfLevels;
       h = step*(1:numOfLevels); 
       
       h = h';
              
       l = 2*7.2921*10^-5*sin(pi/180*0.00001*dat1(5,i));
%        k = find_k( u, v, h, l );
       k_exact = find_k_ver2( u, v, h, l, smallk, bigk, eps, 1 );
    %  k = find_k_ver4( u, v, h, l);
     %  [k, wrong] = find_k_ver3( u, v, h, l, 0, bigk, eps, 0 );
       
%        if wrong == 0
%            num_of_good_k = num_of_good_k + 1
%            table1 (i,:) = k(1,:);
%        else
%            num_of_wrong_k = num_of_wrong_k + 1
%        end
       
       
       break
   end
end

table1 = table1(sum(table1~=-1,2)>15, :);

save([num2str(mon*100+day) '.mat'], 'table1');


table_of_ks = [table_of_ks; table1];
end
end




ng = num_of_good_k;
nw = num_of_wrong_k;

end