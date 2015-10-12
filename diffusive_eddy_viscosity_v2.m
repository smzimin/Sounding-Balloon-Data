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
eps = 0.3*10^-7;
smallk = 0.1*10^-7;
bigk = 5*10^-7;
fin_par = [smallk bigk eps];

alpha = 10^-4;
beta = 2;
reg_params = [alpha  beta];

%% main calculations
for OBSERVATION = 16 
% [ table_of_ks , fr, shir ] = make_table_of_ks_polynoms( reg_params, ...
%     shir_vis, fin_par, [3 9 OBSERVATION] );
%[3 4 11 12 19 22 24 27 28]

[ table_of_ks , fr, shir ] = make_table_of_ks_splines_v2( reg_params, ...
    shir_vis, fin_par, [3 9 OBSERVATION], 0 );

% print(['OBSERVATION_' int2str(OBSERVATION) '_SHIROTA_' int2str(shir) ...
%     '_BETTER_' int2str(fr*1000)],'-dpng')
end