clear
close all
param_values = load('param_values.txt');
poms_torordland(param_values, 'endo');
% poms_torordland(nb_models, 'mid');
% poms_torordland(nb_models,'epi');