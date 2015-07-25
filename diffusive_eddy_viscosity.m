clear;

set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arial Cyr');

highestLevel = 2000;
numOfLevels = highestLevel/20;
step = highestLevel / numOfLevels;
h = step*(1:numOfLevels); 
t = 35:16:4800;
height = 0:step:highestLevel - step;

   
eps = 0.5*10^-2;
smallk = 5.04;
bigk = 5.04;

for shirota = 0
[ left, right, fgr, num ] = switch_length( shirota );
%[ table_of_ks, u, v, ng, nw ] = make_table_of_ks( highestLevel, numOfLevels, t, h, left, right, smallk, bigk, eps );
[ table_of_ks, bu, bv ] = make_table_of_ks2( highestLevel, numOfLevels, t, h, left, right, smallk, bigk, eps, 28 );
%figure(2);
%plot_ks( table_of_ks, h );
end