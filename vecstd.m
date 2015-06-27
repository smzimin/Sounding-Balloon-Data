function temp = vecstd(temp)
a = sind(temp);
        b = cosd(temp);
        
        a = std(a);
        b = std(b);
        
        temp = atan2(a,b)*180/pi;