clear;
close all;

set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arial Cyr');

highestLevel = 2000;
numOfLevels = highestLevel/20;
step = highestLevel / numOfLevels;
h = step*(1:numOfLevels); 

height = 0:step:highestLevel - step;

eps = 10^-12;
smallk = -100; % 3.5559e-06
bigk = 100.1;
heightDivisor = 1000;
alpha = 0*20*10^-3;
beta = 2;

shirota = 0;
counter = 0;
for OBSERVATION = 1:31
[ left, right, fgr, num ] = switch_length( shirota );
%[ table_of_ks, u, v, ng, nw ] = make_table_of_ks( highestLevel, numOfLevels, t, h, left, right, smallk, bigk, eps );
[ table_of_ks , fr, shir ] = make_table_of_ks2(alpha, beta, left, right, smallk, bigk, eps, heightDivisor, OBSERVATION );
print(['OBSERVATION_' int2str(OBSERVATION) '_SHIROTA_' int2str(shir) ...
    '_BETTER_' int2str(fr*1000)],'-dpng')
% figure()
% plot_ks( table_of_ks, h );
if fr > 1
    counter = counter + 1
end
end