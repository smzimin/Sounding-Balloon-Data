function temp = vecmean(temp)
a = sind(temp);
        b = cosd(temp);
        
        a = mean(a);
        b = mean(b);
        
        temp = atan2(a,b)*180/pi;
