clear
close all
param_values = load('uq_param_values.txt');
fid = fopen('uq_output_dir.txt');
output_dir = fgetl(fid);
if ~exist(output_dir)
    mkdir(output_dir)
end
poms_torordland(param_values, 'endo', output_dir);
poms_torordland(param_values, 'mid', output_dir);
poms_torordland(param_values, 'epi', output_dir);
quit()
