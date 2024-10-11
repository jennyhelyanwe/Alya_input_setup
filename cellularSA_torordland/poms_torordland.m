function [] = poms_torordland(param_values,celltype, outputdir)
% Population of models generation for evaluation of sensitivity analysis
% Adapted from scriptDemonstration_2_ParameterComparison.m:

%% Setting parameters
param.bcl = 800; % basic cycle length in ms
param.model = @model_ToRORd_Land; % which model is to be used
param.verbose = false; % printing numbers of beats simulated. 
params(1:size(param_values,1)) = param;
nb_models = size(param_values,1);
for ipopulation = 1:nb_models
    params(ipopulation).INa_Multiplier = param_values(ipopulation,1);
    params(ipopulation).INaL_Multiplier = param_values(ipopulation,2);
    params(ipopulation).Ito_Multiplier = param_values(ipopulation,3);
    params(ipopulation).IKr_Multiplier = param_values(ipopulation,4);
    params(ipopulation).IKs_Multiplier = param_values(ipopulation,5);
    params(ipopulation).IK1_Multiplier = param_values(ipopulation,6);
    params(ipopulation).INCX_Multiplier = param_values(ipopulation,7);
    params(ipopulation).INaK_Multiplier = param_values(ipopulation,8);
    params(ipopulation).ICaL_Multiplier = param_values(ipopulation,9);
    params(ipopulation).Jrel_Multiplier = param_values(ipopulation,10);
    params(ipopulation).Jup_Multiplier = param_values(ipopulation,11);
    params(ipopulation).ca50_Multiplier = param_values(ipopulation,12);
    params(ipopulation).kuw_Multiplier = param_values(ipopulation,13);
    params(ipopulation).kws_Multiplier = param_values(ipopulation,14);
    params(ipopulation).ksu_Multiplier = param_values(ipopulation,15);
%     params(ipopulation).stim_amp = -53 * param_values(ipopulation,16);
end
options = [];
beats = 10;
ignoreFirst = beats - 1;

%% Simulation and output extraction
% Now, the structure of parameters is used to run multiple models in a
% parallel-for loop.
figure;
hold on;
oliTraces = load('data/oliTraces.mat');
for i = 1:size(oliTraces.tracesRealigned,2)
    subplot(1,2,1), plot(oliTraces.referenceTime, oliTraces.tracesRealigned(:,i), 'k');
    hold on;
end
apd40 = zeros(size(params));
apd50 = zeros(size(params));
apd90 = zeros(size(params));
dvdt_max = zeros(size(params));
vpeak = zeros(size(params));
RMP = zeros(size(params));
Tamax = zeros(size(params));
Tamin = zeros(size(params));
TaD50 = zeros(size(params));
TaD90 = zeros(size(params));
dTadt_max = zeros(size(params));
CTD50 = zeros(size(params));
CTD90 = zeros(size(params));
CaTmax = zeros(size(params));
CaTmin = zeros(size(params));

tic
for i = 1:length(params)
%     celltype = 'endo';
    X0 = getStartingState_ToRORdLand(['m_',celltype]); % starting state - can be also m_mid or m_epi for midmyocardial or epicardial cells respectively.
    
    % Simulation and extraction of outputs
    params(i).cellType = celltype; %0 endo, 1 epi, 2 mid
    [t, X] = modelRunner_ToRORdLand(X0, options, params(i), beats, ignoreFirst);
    
    currents = getCurrentsStructure_ToRORdLand(t, X, beats, param, 0);
    time(i).value = currents.time;
    V(i).value = currents.V;
    cai(i).value = currents.Cai;
    Ta(i).value = currents.Ta;
    apd40(i) = DataReporter.getAPD(currents.time, currents.V, 0.4);
    apd50(i) = DataReporter.getAPD(currents.time, currents.V, 0.5);
    apd90(i) = DataReporter.getAPD(currents.time, currents.V, 0.9);
    dvdt_max(i) = DataReporter.getPeakDVDT(currents.time, currents.V, -1.0);
    CTD50(i) = DataReporter.getAPD(currents.time, currents.Cai, 0.5);
    CTD90(i) = DataReporter.getAPD(currents.time, currents.Cai, 0.9);
    CaTmax(i) = max(currents.Cai);
    CaTmin(i) = min(currents.Cai);
    Tamax(i) = max(currents.Ta);
    TaD50(i) = DataReporter.getAPD(currents.time, currents.Ta, 0.5);
    TaD90(i) = DataReporter.getAPD(currents.time, currents.Ta, 0.9);
    dTadt_max(i) = DataReporter.getPeakDVDT(currents.time, currents.Ta, -1.0);
    Tamin(i) = min(currents.Ta);
    vpeak(i) = max(currents.V);
    RMP(i) = min(currents.V);
end  
disp('Time per model: ')
disp(toc/nb_models)
for i = 1:length(params)
     subplot(1,2,1), plot(time(i).value, V(i).value, 'b');
     hold on;
end

tri_90_40 = apd90 - apd40;
for i = 1:length(params)
    subplot(1,2,2), plot(time(i).value, Ta(i).value, 'b');
    hold on;
end


%%
% %% Calibration based on Table 2 of Passini et al. 2019 https://doi.org/10.1111/bph.14786
% apd40_min = 85;
% apd40_max = 320;
% apd50_min = 110;
% apd50_max = 350;
% apd90_min = 180;
% apd90_max = 440;
% tri_90_40_min = 50;
% tri_90_40_max = 150;
% dvdt_max_min = 100;
% dvdt_max_max = 1000;
% vpeak_min = 10;
% vpeak_max = 55;
% rmp_min = -95;
% rmp_max = -80;
% ctd50_min = 120;
% ctd50_max = 420;
% ctd90_min = 220;
% ctd90_max = 785;
% 
% calibration_criteria_min = [apd40_min, apd50_min, apd90_min, tri_90_40_min,...
%     dvdt_max_min, vpeak_min, rmp_min, ctd50_min, ctd90_min];
% calibration_criteria_max = [apd40_max, apd50_max, apd90_max, tri_90_40_max,...
%     dvdt_max_max, vpeak_max, rmp_max, ctd50_max, ctd90_max];
% population_biomarkers = [apd40; apd50; apd90; tri_90_40; ...
%     dvdt_max; vpeak; RMP; CTD50; CTD90]';
% calibrated_population_biomarkers_index = (population_biomarkers > calibration_criteria_min) & ( population_biomarkers < calibration_criteria_max);
% calibrated_population_biomarkers_index = all(calibrated_population_biomarkers_index, 2);
% calibrated_population_V = V(calibrated_population_biomarkers_index);
% calibrated_population_Ta = Ta(calibrated_population_biomarkers_index);
% calibrated_population_t = time(calibrated_population_biomarkers_index);
% 
% % Save as table based on integer values of APD90 and APD50 (based on
% % variability) in ms. 
% APD40 = apd40(calibrated_population_biomarkers_index)';
% APD50 = apd50(calibrated_population_biomarkers_index)';
% APD90 = apd90(calibrated_population_biomarkers_index)';
% CTD50 = CTD50(calibrated_population_biomarkers_index)';
% CTD90 = CTD90(calibrated_population_biomarkers_index)';
% CaT_max = CaTmax(calibrated_population_biomarkers_index)';
% CaT_min = CaTmin(calibrated_population_biomarkers_index)';
% Ta_max = Tamax(calibrated_population_biomarkers_index)';
% Ta_min = Tamin(calibrated_population_biomarkers_index)';
T = table(apd40', apd50', apd90', CTD50', CTD90', CaTmax', CaTmin', Tamax', ...
          Tamin', TaD50', TaD90', dTadt_max', 'VariableNames', ...
          {'APD40', 'APD50', 'APD90', 'CTD50', 'CTD90', 'CaTmax', 'CaTmin', ...
          'Tamax', 'Tamin', 'TaD50', 'TaD90', 'dTadtmax'});
writetable(T,[outputdir celltype '_output.txt'],'WriteRowNames',true, 'WriteVariableNames', true)
save([outputdir 'V.mat'], 'V');
save([outputdir 'cai.mat'], 'cai');
save([outputdir 'Ta.mat'], 'Ta');
save([outputdir 'time.mat'], 'time');
end
