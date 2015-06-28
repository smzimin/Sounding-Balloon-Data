function plot_ks(table_of_ks, h)

k_av = zeros(1,size(table_of_ks,2));

for i= 1:size(table_of_ks,2)
    k_av(i) = mean(table_of_ks(:,i));
end

plot(k_av,h, 'LineWidth', 2.5, 'LineSmoothing', 'on');
xlabel('k(z)','FontName','Arial Cyr','FontSize',20);
ylabel('Высота, метры','FontName','Arial Cyr','FontSize',20);
end