function [ numOfTimes, matrT, matrE, matrD, interval ] = calculate_means_and_conf_interval( angle, numOfLevels )

numOfTimes = length(angle(:,1));
matrT = zeros(numOfLevels,numOfLevels,numOfTimes);
matrE = zeros(numOfLevels,numOfLevels);
matrD = zeros(numOfLevels,numOfLevels);
interval = zeros(2,numOfLevels);

for i = 1 %numOfLevels
    for k = 1:numOfLevels
        
        for p = 1:numOfTimes    
            matrT(i,k,p) = angle(p,k) - angle(p,i);
            if matrT(i,k,p) > 180
                matrT(i,k,p) = matrT(i,k,p) - 360;
            elseif matrT(i,k,p) < -180
                matrT(i,k,p) = matrT(i,k,p) + 360;
            end
            if angle(p,k) == -1 || angle(p,i) == -1
                matrT(i,k,p) = 400;
            end;
        end                
               
        temp(1,:) = matrT(i,k,:);
        
         r1 = find(abs(temp) < 400);  
         temp = temp(r1);

        
        matrE(i,k) = vecmean(temp);
        

        matrD(i,k) = vecstd(temp);
        
        if i == 1 && k ~= 1
            boost = bootstrp(1000,@vecmean,temp);
            boost = sort(boost);
            interval(1,k) = boost(26);
            interval(2,k) = boost(975);
        end
        
        clear temp;
    end   
end

end

