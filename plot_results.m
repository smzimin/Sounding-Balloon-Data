function plot_results( shirota, fgr, num, matrE, matrT, height, numOfLevels, interval )
if shirota == 0
    figure(1);
    plot(matrE(1,:), height, 'LineWidth', 1.5); hold on;
    for k = 1:numOfLevels  
        plot([interval(1,k) interval(2,k)],[height(k) height(k)],'r*','LineWidth',1);
        plot([interval(1,k) interval(2,k)],[height(k) height(k)],'r','LineWidth',1);  
    end
    xlabel('Угол поворота, градусы','FontName','Arial Cyr','FontSize',20);
    ylabel('Высота, метры','FontName','Arial Cyr','FontSize',20);
    
else
    figure(fgr);
    subplot (1,2,num)
    plot(matrE(1,:),height,'LineWidth',1.5); hold on;
    for k = 1:numOfLevels  
        plot([interval(1,k) interval(2,k)],[height(k) height(k)],'r*','LineWidth',1);
        plot([interval(1,k) interval(2,k)],[height(k) height(k)],'r','LineWidth',1);  
    end
    xlabel('Угол поворота, градусы','FontName','Arial Cyr','FontSize',20);
    ylabel('Высота, метры','FontName','Arial Cyr','FontSize',20);
    % l = legend('Выборочное среднее','Доверительный интервал');
    % set(l,'FontSize',15);

    figure(fgr + 3);
    subplot(1,2,num) 
    % hold off
    % clear r1
    % clear temp
    temp(1,:) = matrT(1,numOfLevels,:);
    r1 = find(abs(temp) < n);
    temp = temp(r1);
    size(temp)
    % h = hist(temp,x);
    % hist(temp,x); hold on;
    % xlabel('Угол поворота, градусы','FontName','Arial Cyr','FontSize',25);

    m = vecmean(temp);
    % plot([m m],[0 max(h) + 4],'r*','LineWidth',2)
    % plot([m m],[0 max(h) + 4],'r','LineWidth',2)

    plot(interval(2,:) - interval(1,:),height);
    %clear temp;
end
end

