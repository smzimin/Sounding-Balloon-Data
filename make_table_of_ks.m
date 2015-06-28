function [ table_of_ks ] = make_table_of_ks( highestLevel, numOfLevels, t, h, left, right )
table_of_ks = [];
for mon = 3:12
for day = 1:31
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
    
       x = cputime;
       u = spline(p, speed .* cos(a));
       du1 = fnder(u, 1);
       du2 = fnder(u, 2);
       du3 = fnder(u, 3);

       v = spline(p, speed .* sin(a));
       dv1 = fnder(v, 1);
       dv2 = fnder(v, 2);
       dv3 = fnder(v, 3);
       
       l = 2*7.2921*10^-5*sin(0.00001*dat1(5,i));
       solinit = bvpinit(h, [35 0]);
       sol = bvp4c(@odefun, @twobc, solinit);
       k = deval(sol, h);
       k(:,1)
       table1 (i,:) = k(1,:);
   end
end

table1 = table1(sum(table1~=-1,2)>15, :);

save([num2str(mon*100+day) '.mat'], 'table1');


table_of_ks = [table_of_ks; table1];
end
end



function res = twobc(ya,yb)

ua = ppval(u, h(1));
va = ppval(v, h(1));
d1ua = ppval(du1, h(1));
d1va = ppval(dv1, h(1));
d2ua = ppval(du2, h(1));
d2va = ppval(dv2, h(1));

ub = ppval(u, h(end));
vb = ppval(v, h(end));
d1ub = ppval(du1, h(end));
d1vb = ppval(dv1, h(end));
d2ub = ppval(du2, h(end));
d2vb = ppval(dv2, h(end));

res = [(d1ua^2+d1va^2)*ya(2) + (d2ua*d1ua+d2va*d1va)*ya(1) + l*(ua*d1va-va*d1ua) ;
    (d1ub^2+d1vb^2)*yb(2) + (d2ub*d1ub+d2vb*d1vb)*yb(1) + l*(ub*d1vb-vb*d1ub);];

end

function dydx = odefun(x,y)
    
d1uz = ppval(du1, x);
d1vz = ppval(dv1, x);
d2uz = ppval(du2, x);
d2vz = ppval(dv2, x);
d3uz = ppval(du3, x);
d3vz = ppval(dv3, x);

dydx = [y(2);
   % ppval(pp1, x)*(y(1)-ppval(pp2, x))-ppval(pp2, x)-0.1];
    -(2*(d2uz*d1uz+d2vz*d1vz)*y(2) + (d3uz*d1uz+d3vz*d1vz)*y(1) ) / (d1uz^2+d1vz^2)];

end


end

