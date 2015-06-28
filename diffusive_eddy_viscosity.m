clear;

set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arial Cyr');

highestLevel = 2000;
numOfLevels = highestLevel/50;
step = highestLevel / numOfLevels;
h = step*(1:numOfLevels);
t = 35:16:4800;
height = 0:step:highestLevel - step;

for shirota = 0
[ left, right, fgr, num ] = switch_length( shirota );
table_of_ks = make_table_of_ks( highestLevel, numOfLevels, t, h, left, right );
%[ numOfTimes, matrT, matrE, matrD, interval ] = calculate_means_and_conf_interval( angle, numOfLevels );
plot_ks( table_of_ks, h );
end