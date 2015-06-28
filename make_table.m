function [ table ] = make_table( highestLevel, numOfLevels, t, h, left, right )
table = [];
for mon = 3:12
    mon
for day = 1:31
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
       sel = p > -1000 & dat1(t+4,i) > 20 & p < highestLevel*1.5;
       if sum (sel) < 8 || min (p(sel)) > 250 || max(p(sel)) < highestLevel-250
           continue;
       end
       a = a(sel);
       p = p(sel);       
       sel = p (2:end) > p(1:end-1);
       sel = [true; sel];     
       if sum (sel) < 8
           continue;
       end  
       a = a(sel);
       p = p(sel);
       
       %test = p;
   
       s = interp1 (p,sin(a),h,'linear','extrap');
       c = interp1 (p,cos(a),h,'linear','extrap');
       a1 = 180/pi*atan2 (s,c);
       a1 (a1 < 0) = a1 (a1 < 0) + 360;
       table1 (i,:) = a1;
   end
end

table1 = table1(sum(table1~=-1,2)>15, :);

table = [table; table1];
end
end


end

