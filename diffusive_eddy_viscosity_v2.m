clear;
close all;

set(0,'DefaultAxesFontSize',20,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontSize',12,'DefaultTextFontName','Arial Cyr');

highestLevel = 1000;
heightDivisor = 1000;
shirota = 0;
[ left, right ] = switch_length( shirota );
shir_vis = [left right highestLevel heightDivisor];

%% params for finding k(z)
eps = 10^-12;
smallk = -100;
bigk = 100.1;
fin_par = [smallk bigk eps];

alpha = 0*20*10^-3;
beta = 2;
reg_params = [alpha  beta];

%% main calculations
for OBSERVATION = 1
[ table_of_ks , fr, shir ] = make_table_of_ks_polynoms( reg_params, ...
    shir_vis, fin_par, [3 9 OBSERVATION] );

% [ table_of_ks , fr, shir ] = make_table_of_ks_splines( reg_params, ...
%     shir_vis, fin_par, [3 9 OBSERVATION] );

% print(['OBSERVATION_' int2str(OBSERVATION) '_SHIROTA_' int2str(shir) ...
%     '_BETTER_' int2str(fr*1000)],'-dpng')

end