clear
close all
param_values = load('param_values.txt');
poms_torordland(param_values, 'endo');
poms_torordland(param_values, 'mid');
poms_torordland(param_values, 'epi');