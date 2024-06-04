from alyaformat import AlyaFormat
from postprocessing import PostProcessing
from sensitivityanalysis_uncertaintyquantification import SAUQ
from generatefields import FieldGeneration
import os
import numpy as np
import json
import pandas as pd
import code

########################################################################################################################
# Global Settings
mesh_number = '05'
simulation_name = 'rodero_' + mesh_number + '_fine'
workdir = os.getcwd()
if 'icei' in workdir:
    system = 'jureca'
elif 'cosma' in workdir:
    system = 'cosma'
elif 'polaris' in workdir:
    system = 'polaris'
elif 'Expansion' in workdir:
    system = 'archive'
elif 'e769' in workdir:
    system = 'archer2'
else:
    system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'cosma':
    meta_data_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
elif system == 'archive':
    meta_data_dir = '/run/media/jenang/Expansion/JURECA_COSMA_download_April2024/meta_data/'
elif system == 'archer2':
    meta_data_dir = '/work/e769/e769/jennywang/Alya_pipeline/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False

########################################################################################################################
# Step 1: Initialise Alya input files writing capabilities.
simulation_root_dir = ''
if system == 'jureca':
    simulation_root_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'cosma':
    simulation_root_dir = '/snap8/scratch/dp287/dc-wang14/alya_simulations/'
elif system == 'polaris':
    simulation_root_dir = '/grand/projects/CompBioAffin/jenang/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_root_dir = './'
elif system == 'archive':
    simulation_root_dir = '/run/media/jenang/Expansion/JURECA_COSMA_download_April2024/Alya_pipeline/alya_simulations/'
elif system == 'archer2':
    simulation_root_dir = '/work/e769/e769/jennywang/Alya_pipeline/alya_simulations/'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_root_dir, verbose=verbose, job_version=system)
########################################################################################################################
# CHANGE THIS FOR DIFFERENT SAs!!!
setup_simulations = False
run_simulations = False
# Choose which groups of parameters to setup/run/evaluate
passive_mechanics = True
active_mechanics = True
cellular = True
haemodynamic = True
all_parameters_at_once = False

# Choose which groups of QoI to evaluate
evaluate_pv= True
evaluate_ecg = True
evaluate_deformation = False
evaluate_fibrework = False
evaluate_strain = False
evaluate_volume = False
fresh_qoi_evaluation = False

parameter_names = []
sa_folder_root_names = []
cellular_params = ['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup']
active_params = ['tref_scaling_myocardium', 'cal50_myocardium', 'sfkws_myocardium']
passive_params = ['pericardial_stiffness', 'Kct_myocardium', 'a_myocardium', 'af_myocardium', 'as_myocardium', 'afs_myocardium']
haemo_params = ['arterial_resistance_lv',
                                'arterial_compliance_lv',
                                'gain_error_relaxation_lv',
                                'gain_derror_relaxation_lv',
                                'ejection_pressure_threshold_lv']
if cellular:
    parameter_names = parameter_names + cellular_params
    sa_folder_root_names = sa_folder_root_names + ['sensitivity_analyses_cellular_parameters_oat']*len(cellular_params)
if active_mechanics:
    parameter_names = parameter_names + active_params
    sa_folder_root_names = sa_folder_root_names + ['sensitivity_analyses_active_mechanical_parameters_oat']*len(active_params)
if passive_mechanics:
    parameter_names = parameter_names + passive_params
    sa_folder_root_names = sa_folder_root_names + ['sensitivity_analyses_mechanical_parameters_oat']*len(passive_params)
if haemodynamic:
    parameter_names = parameter_names + haemo_params
    sa_folder_root_names = sa_folder_root_names + ['sensitivity_analyses_haemodynamics_parameters_oat']*len(haemo_params)
else:
    print('Turn on one of the SA types!')
    quit()
print('Parameters to evaluate: ', parameter_names)
print('SA root folder names: ', sa_folder_root_names)
########################################################################################################################
baseline_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
########################################################################################################################
# Run simulations
if setup_simulations or run_simulations:
    for param_i, param in enumerate(parameter_names):
        print('Dealing with parameter: ', param)
        if 'sf_' in param:
            baseline_parameter_values = np.array([simulation_dict[param][0][0]])
        elif '_lv' in param:
            baseline_parameter_values = np.array([simulation_dict[param.split('_lv')[0]][0]])
        elif '_rv' in param:
           baseline_parameter_values = np.array([simulation_dict[param.split('_rv')[0]][1]])
        elif '_myocardium' in param:
            baseline_parameter_values = np.array([simulation_dict[param.split('_myocardium')[0]][0]])
        elif '_valveplug' in param:
            baseline_parameter_values = np.array([simulation_dict[param.split('_valveplug')[0]][1]])
        else:
            baseline_parameter_values = np.array([simulation_dict[param]])
        upper_bounds = baseline_parameter_values * 2.0
        lower_bounds = baseline_parameter_values * 0.1
        sampling_method = 'uniform'
        if passive_mechanics:
            upper_bounds = baseline_parameter_values * 100.0
            lower_bounds = baseline_parameter_values * 0.01
            sampling_method = ''
        baseline_json_file = 'rodero_baseline_simulation_em.json'
        simulation_dir = ''
        baseline_dir = ''
        if system == 'jureca':
            baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
            simulation_dir = simulation_root_dir + sa_folder_root_names[param_i] + '_' + param + '/'
        elif system == 'cosma':
            baseline_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
            simulation_dir = simulation_root_dir + sa_folder_root_names[param_i] + '_' + param + '/'
        elif system == 'polaris':
            baseline_dir = '/grand/projects/CompBioAffin/jenang/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
        elif system == 'heart':
            baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
            simulation_dir = sa_folder_root_names[param_i] + '_' + param + '/'
        elif system == 'archer2':
            baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
            simulation_dir = simulation_root_dir + sa_folder_root_names[param_i] + '/'
        elif system == 'archive':
            baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
            simulation_dir = simulation_root_dir + sa_folder_root_names[param_i] + '/'
        parameter_names = np.array([param])
        sa = SAUQ(name='sa', sampling_method='uniform', n=8, parameter_names=parameter_names, baseline_json_file=baseline_json_file,
                 simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
        if setup_simulations:
            sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
        if run_simulations:
            sa.run_jobs(simulation_dir)

########################################################################################################################
# Evaluate QoIs and correlations
qoi_names = []
pv_qois = ['EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection', 'dvdt_filling', 'dpdt_max', 'EDVR', 'ESVR',
                     'PmaxR', 'SVR']
ecg_qois = ['qrs_dur_mean', 't_dur_mean', 'qt_dur_mean', 't_pe_mean', 'jt_dur_mean']
deformation_qois = ['es_ed_avpd', 'es_ed_apical_displacement', 'diff_lv_wall_thickness']
fibrework_qois = ['peak_lambda', 'min_lambda', 'peak_ta', 'diastolic_ta']
strain_qois = ['max_mid_Ecc', 'min_mid_Ecc', 'max_mid_Err', 'min_mid_Err', 'max_four_chamber_Ell', 'min_four_chamber_Ell']
if evaluate_ecg:
    qoi_names= qoi_names + ecg_qois
if evaluate_pv:
    qoi_names = qoi_names + pv_qois
if evaluate_deformation:
    qoi_names= qoi_names + deformation_qois
if evaluate_fibrework:
    qoi_names= qoi_names + fibrework_qois
if evaluate_strain:
    qoi_names= qoi_names + strain_qois
print('Setting up emtpy SA output dataframe:')
all_slopes = pd.DataFrame(columns=qoi_names, index=parameter_names)
all_intercepts = pd.DataFrame(columns=qoi_names, index=parameter_names)
all_p_values = pd.DataFrame(columns=qoi_names, index=parameter_names)
all_r_values = pd.DataFrame(columns=qoi_names, index=parameter_names)
all_ranges = pd.DataFrame(columns=qoi_names, index=parameter_names)
print(all_slopes)
for param_i, param in enumerate(parameter_names):
    sa_folder_root_name = sa_folder_root_names[param_i]
    simulation_dir = ''
    if system == 'jureca':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'cosma':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'polaris':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'heart':
        simulation_dir = sa_folder_root_name + '_' + param + '/'
    elif system == 'archive':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'archer2':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    if 'sf_' in param:
        baseline_parameter_values = np.array([simulation_dict[param][0][0]])
    elif '_lv' in param:
        baseline_parameter_values = np.array([simulation_dict[param.split('_lv')[0]][0]])
    elif '_rv' in param:
        baseline_parameter_values = np.array([simulation_dict[param.split('_rv')[0]][1]])
    elif '_myocardium' in param:
        baseline_parameter_values = np.array([simulation_dict[param.split('_myocardium')[0]][0]])
    elif '_valveplug' in param:
        baseline_parameter_values = np.array([simulation_dict[param.split('_valveplug')[0]][1]])
    else:
        baseline_parameter_values = np.array([simulation_dict[param]])
    upper_bounds = baseline_parameter_values * 2.0
    lower_bounds = baseline_parameter_values * 0.1
    parameter_names = np.array([param])  # Dummy
    baseline_dir = ''  # Dummy
    sa = SAUQ(name='sa', sampling_method='uniform', n=8, parameter_names=parameter_names, baseline_json_file=baseline_json_file,
              simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
    labels = [param +'=' + s for s in ["%.1f" % x for x in np.linspace(lower_bounds, upper_bounds, 8)]]
    beat = 1
    ####################################################################################################################
    # Can be done as the simulations are running, using raw outputs.
    if system == 'archive':
        sa.sort_simulations_archive(tag='raw')
    else:
        sa.sort_simulations(
            tag='raw')  # Collate list of finished simulations by checking the existence of particular files.
    print('Number of finished simulations: ', len(sa.finished_simulation_dirs))
    if evaluate_pv:
        # Pressure volume information
        if fresh_qoi_evaluation:
            pv_post = sa.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                       analysis_type='sa')
            sa.visualise_sa(beat=1, pv_post=pv_post, labels=labels, save_filename=simulation_dir+'/pv_post.png')
        qoi_names = pv_qois
        slopes, intercepts, p_values, r_values, ranges, params, qois = sa.analyse(filename=simulation_dir + 'pv_qois.csv', qois=qoi_names, show_healthy_ranges=False,
                                   save_filename=simulation_dir + '/pv_scatter.png')
        all_slopes.loc[param, qoi_names] = slopes[:, 0]
        all_intercepts.loc[param, qoi_names] = intercepts[:, 0]
        all_p_values.loc[param, qoi_names] = p_values[:, 0]
        all_r_values.loc[param, qoi_names] = r_values[:, 0]
        all_ranges.loc[param, qoi_names] = ranges[:, 0]
        # slopes = dict(map(lambda i, j: (i, j), qoi_names, slopes[:, 0]))
        # intercepts = dict(map(lambda i, j: (i, j), qoi_names, intercepts[:, 0]))
        # p_values = dict(map(lambda i, j: (i, j), qoi_names, p_values[:, 0]))
        # r_values = dict(map(lambda i, j: (i, j), qoi_names, r_values[:, 0]))
        # ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
        # if all_parameters_at_once:
        #     all_slopes.append(pd.DataFrame(slopes, index=[param]))
        #     all_intercepts.append(pd.DataFrame(intercepts, index=[param]))
        #     all_p_values.append(pd.DataFrame(p_values, index=[param]))
        #     all_r_values.append(pd.DataFrame(r_values, index=[param]))
        #     all_ranges.append(pd.DataFrame(ranges, index=[param]))
        # pd.DataFrame(slopes, index=[param]).to_csv(simulation_dir + '/pv_slopes.csv')
        # pd.DataFrame(intercepts, index=[param]).to_csv(simulation_dir + '/pv_intercepts.csv')
        # pd.DataFrame(p_values, index=[param]).to_csv(simulation_dir + '/pv_p_values.csv')
        # pd.DataFrame(r_values, index=[param]).to_csv(simulation_dir + '/pv_r_values.csv')
        # pd.DataFrame(ranges, index=[param]).to_csv(simulation_dir + '/pv_ranges.csv')
        # pd.DataFrame(ranges, index=[param]).to_csv(simulation_dir + '/pv_ranges.csv')
        np.savetxt(simulation_dir + '/param_'+param+'_inputs.csv', params, delimiter=',')
        pd.DataFrame(qois, columns=qoi_names).to_csv(simulation_dir + '/param_' + param + '_pv_qoi_outcomes.csv')

    # ECG information
    if evaluate_ecg:
        if fresh_qoi_evaluation:
            ecg_post = sa.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                        analysis_type='sa')
            if not ecg_post:
                print('ERROR: ECG postprocessing did not happen for some reason')
                quit()
            sa.visualise_sa(beat=1, ecg_post=ecg_post, labels=labels, save_filename=simulation_dir+'/ecg_post.png')
        qoi_names = ecg_qois
        slopes, intercepts, p_values, r_values, ranges, params, qois = sa.analyse(filename=simulation_dir + 'ecg_qois.csv', qois=qoi_names, show_healthy_ranges=False,
                                   save_filename=simulation_dir + '/ecg_scatter.png')
        all_slopes.loc[param, qoi_names] = slopes[:, 0]
        all_intercepts.loc[param, qoi_names] = intercepts[:, 0]
        all_p_values.loc[param, qoi_names] = p_values[:, 0]
        all_r_values.loc[param, qoi_names] = r_values[:, 0]
        all_ranges.loc[param, qoi_names] = ranges[:, 0]
        # slopes = dict(map(lambda i, j: (i, j), qoi_names, slopes[:, 0]))
        # intercepts = dict(map(lambda i, j: (i, j), qoi_names, intercepts[:, 0]))
        # p_values = dict(map(lambda i, j: (i, j), qoi_names, p_values[:, 0]))
        # r_values = dict(map(lambda i, j: (i, j), qoi_names, r_values[:, 0]))
        # ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
        # if all_parameters_at_once:
        #     all_slopes.append(pd.DataFrame(slopes, index=[param]))
        #     all_intercepts.append(pd.DataFrame(intercepts, index=[param]))
        #     all_p_values.append(pd.DataFrame(p_values, index=[param]))
        #     all_r_values.append(pd.DataFrame(r_values, index=[param]))
        #     all_ranges.append(pd.DataFrame(ranges, index=[param]))
        # pd.DataFrame(slopes, index=[param]).to_csv(simulation_dir + '/ecg_slopes.csv')
        # pd.DataFrame(intercepts, index=[param]).to_csv(simulation_dir + '/ecg_intercepts.csv')
        # pd.DataFrame(p_values, index=[param]).to_csv(simulation_dir + '/ecg_p_values.csv')
        # pd.DataFrame(r_values, index=[param]).to_csv(simulation_dir + '/ecg_r_values.csv')
        # pd.DataFrame(ranges, index=[param]).to_csv(simulation_dir + '/ecg_ranges.csv')
        np.savetxt(simulation_dir + '/param_' + param + '_inputs.csv', params, delimiter=',')
        pd.DataFrame(qois, columns=qoi_names).to_csv(simulation_dir + '/param_' + param + '_ecg_qoi_outcomes.csv')

    # Deformation
    if evaluate_deformation:
        if fresh_qoi_evaluation:
            deformation_post = sa.evaluate_qois(qoi_group_name='deformation', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                                analysis_type='sa')
            sa.visualise_sa(beat=1, deformation_post=deformation_post, labels=labels,
                            save_filename=simulation_dir + '/deformation_post.png')
        qoi_names = ['es_ed_avpd', 'es_ed_apical_displacement', 'diff_lv_wall_thickness']
        slopes, intercepts, p_values, r_values, ranges, params, qois = sa.analyse(filename=simulation_dir + 'deformation_qois.csv', qois=qoi_names,
                                   show_healthy_ranges=False, save_filename=simulation_dir + '/deformation_scatter.png')
        all_slopes.loc[param, qoi_names] = slopes[:, 0]
        all_intercepts.loc[param, qoi_names] = intercepts[:, 0]
        all_p_values.loc[param, qoi_names] = p_values[:, 0]
        all_r_values.loc[param, qoi_names] = r_values[:, 0]
        all_ranges.loc[param, qoi_names] = ranges[:, 0]
        # slopes = dict(map(lambda i, j: (i, j), qoi_names, slopes[:, 0]))
        # intercepts = dict(map(lambda i, j: (i, j), qoi_names, intercepts[:, 0]))
        # p_values = dict(map(lambda i, j: (i, j), qoi_names, p_values[:, 0]))
        # r_values = dict(map(lambda i, j: (i, j), qoi_names, r_values[:, 0]))
        # ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
        # if all_parameters_at_once:
        #     all_slopes.append(pd.DataFrame(slopes, index=[param]))
        #     all_intercepts.append(pd.DataFrame(intercepts, index=[param]))
        #     all_p_values.append(pd.DataFrame(p_values, index=[param]))
        #     all_r_values.append(pd.DataFrame(r_values, index=[param]))
        #     all_ranges.append(pd.DataFrame(ranges, index=[param]))
        # pd.DataFrame(slopes, index=[param]).to_csv(simulation_dir + '/deformation_slopes.csv')
        # pd.DataFrame(intercepts, index=[param]).to_csv(simulation_dir + '/deformation_intercepts.csv')
        # pd.DataFrame(p_values, index=[param]).to_csv(simulation_dir + '/deformation_p_values.csv')
        # pd.DataFrame(r_values, index=[param]).to_csv(simulation_dir + '/deformation_r_values.csv')
        # pd.DataFrame(ranges, index=[param]).to_csv(simulation_dir + '/deformation_ranges.csv')
        np.savetxt(simulation_dir + '/param_' + param + '_inputs.csv', params, delimiter=',')
        pd.DataFrame(qois, columns=qoi_names).to_csv(simulation_dir + '/param_' + param + '_deformation_qoi_outcomes.csv')

    if evaluate_volume:
        if fresh_qoi_evaluation:
            volume_post = sa.evaluate_qois(qoi_group_name='volume', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                                analysis_type='sa')
            sa.visualise_sa(beat=1, volume_post=volume_post, labels=labels)
                            # save_filename=simulation_dir + '/volume_post.png')
        qoi_names = ['percentage_volume_change']
        slopes, intercepts, p_values, r_values, ranges, params, qois = sa.analyse(filename=simulation_dir + 'volume_qois.csv', qois=qoi_names,
                                   show_healthy_ranges=False, save_filename=simulation_dir + '/volume_scatter.png')
        all_slopes.loc[param, qoi_names] = slopes[:, 0]
        all_intercepts.loc[param, qoi_names] = intercepts[:, 0]
        all_p_values.loc[param, qoi_names] = p_values[:, 0]
        all_r_values.loc[param, qoi_names] = r_values[:, 0]
        all_ranges.loc[param, qoi_names] = ranges[:, 0]
        # slopes = dict(map(lambda i, j: (i, j), qoi_names, slopes[:, 0]))
        # intercepts = dict(map(lambda i, j: (i, j), qoi_names, intercepts[:, 0]))
        # p_values = dict(map(lambda i, j: (i, j), qoi_names, p_values[:, 0]))
        # r_values = dict(map(lambda i, j: (i, j), qoi_names, r_values[:, 0]))
        # ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
        # if all_parameters_at_once:
        #     all_slopes.append(pd.DataFrame(slopes, index=[param]))
        #     all_intercepts.append(pd.DataFrame(intercepts, index=[param]))
        #     all_p_values.append(pd.DataFrame(p_values, index=[param]))
        #     all_r_values.append(pd.DataFrame(r_values, index=[param]))
        #     all_ranges.append(pd.DataFrame(ranges, index=[param]))
        # pd.DataFrame(slopes, index=[param]).to_csv(simulation_dir + '/volume_slopes.csv')
        # pd.DataFrame(intercepts, index=[param]).to_csv(simulation_dir + '/volume_intercepts.csv')
        # pd.DataFrame(p_values, index=[param]).to_csv(simulation_dir + '/volume_p_values.csv')
        # pd.DataFrame(r_values, index=[param]).to_csv(simulation_dir + '/volume_r_values.csv')
        # pd.DataFrame(ranges, index=[param]).to_csv(simulation_dir + '/volume_ranges.csv')
        np.savetxt(simulation_dir + '/param_' + param + '_inputs.csv', params, delimiter=',')
        pd.DataFrame(qois, columns=qoi_names).to_csv(simulation_dir + '/param_' + param + '_volume_qoi_outcomes.csv')


    # Fibre strain and Ta
    if evaluate_fibrework:
        if fresh_qoi_evaluation:
            fibre_work_post = sa.evaluate_qois(qoi_group_name='fibre_work', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                               analysis_type='sa')
            sa.visualise_sa(beat=1, fibre_work_post=fibre_work_post, labels=labels,
                            save_filename=simulation_dir + '/fibre_work_post.png')
        qoi_names = ['peak_lambda', 'min_lambda', 'peak_ta', 'diastolic_ta']
        slopes, intercepts, p_values, r_values, ranges, params, qois = sa.analyse(filename=simulation_dir + 'fibrework_qois.csv', qois=qoi_names,
                                   show_healthy_ranges=False, save_filename=simulation_dir + '/fibre_scatter.png')
        all_slopes.loc[param, qoi_names] = slopes[:, 0]
        all_intercepts.loc[param, qoi_names] = intercepts[:, 0]
        all_p_values.loc[param, qoi_names] = p_values[:, 0]
        all_r_values.loc[param, qoi_names] = r_values[:, 0]
        all_ranges.loc[param, qoi_names] = ranges[:, 0]
        # slopes = dict(map(lambda i, j: (i, j), qoi_names, slopes[:, 0]))
        # intercepts = dict(map(lambda i, j: (i, j), qoi_names, intercepts[:, 0]))
        # p_values = dict(map(lambda i, j: (i, j), qoi_names, p_values[:, 0]))
        # r_values = dict(map(lambda i, j: (i, j), qoi_names, r_values[:, 0]))
        # ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
        # if all_parameters_at_once:
        #     all_slopes.append(pd.DataFrame(slopes, index=[param]))
        #     all_intercepts.append(pd.DataFrame(intercepts, index=[param]))
        #     all_p_values.append(pd.DataFrame(p_values, index=[param]))
        #     all_r_values.append(pd.DataFrame(r_values, index=[param]))
        #     all_ranges.append(pd.DataFrame(ranges, index=[param]))
        # pd.DataFrame(slopes, index=[param]).to_csv(simulation_dir + '/fibre_slopes.csv')
        # pd.DataFrame(intercepts, index=[param]).to_csv(simulation_dir + '/fibre_intercepts.csv')
        # pd.DataFrame(p_values, index=[param]).to_csv(simulation_dir + '/fibre_p_values.csv')
        # pd.DataFrame(r_values, index=[param]).to_csv(simulation_dir + '/fibre_r_values.csv')
        # pd.DataFrame(ranges, index=[param]).to_csv(simulation_dir + '/fibre_ranges.csv')
        np.savetxt(simulation_dir + '/param_' + param + '_inputs.csv', params, delimiter=',')
        pd.DataFrame(qois, columns=qoi_names).to_csv(simulation_dir + '/param_' + param + '_fibre_qoi_outcomes.csv')

    # Strain information
    if evaluate_strain:
        if fresh_qoi_evaluation:
            strain_post = sa.evaluate_qois(qoi_group_name='strain', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                           analysis_type='sa')
            sa.visualise_sa(beat=1, strain_post=strain_post, labels=labels,
                            save_filename=simulation_dir + '/strain_post.png')
        qoi_names = ['max_mid_Ecc', 'min_mid_Ecc', 'max_mid_Err', 'min_mid_Err', 'max_four_chamber_Ell', 'min_four_chamber_Ell']
        slopes, intercepts, p_values, r_values, ranges, params, qois = sa.analyse(filename=simulation_dir + 'strain_qois.csv', qois=qoi_names, show_healthy_ranges=False,
                                   save_filename=simulation_dir + '/strain_scatter.png')
        all_slopes.loc[param, qoi_names] = slopes[:, 0]
        all_intercepts.loc[param, qoi_names] = intercepts[:, 0]
        all_p_values.loc[param, qoi_names] = p_values[:, 0]
        all_r_values.loc[param, qoi_names] = r_values[:, 0]
        all_ranges.loc[param, qoi_names] = ranges[:, 0]
        # slopes = dict(map(lambda i, j: (i, j), qoi_names, slopes[:, 0]))
        # intercepts = dict(map(lambda i, j: (i, j), qoi_names, intercepts[:, 0]))
        # p_values = dict(map(lambda i, j: (i, j), qoi_names, p_values[:, 0]))
        # r_values = dict(map(lambda i, j: (i, j), qoi_names, r_values[:, 0]))
        # ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
        # if all_parameters_at_once:
        #     all_slopes.append(pd.DataFrame(slopes, index=[param]))
        #     all_intercepts.append(pd.DataFrame(intercepts, index=[param]))
        #     all_p_values.append(pd.DataFrame(p_values, index=[param]))
        #     all_r_values.append(pd.DataFrame(r_values, index=[param]))
        #     all_ranges.append(pd.DataFrame(ranges, index=[param]))
        # pd.DataFrame(slopes, index=[param]).to_csv(simulation_dir + '/strain_slopes.csv')
        # pd.DataFrame(intercepts, index=[param]).to_csv(simulation_dir + '/strain_intercepts.csv')
        # pd.DataFrame(p_values, index=[param]).to_csv(simulation_dir + '/strain_p_values.csv')
        # pd.DataFrame(r_values, index=[param]).to_csv(simulation_dir + '/strain_r_values.csv')
        # pd.DataFrame(ranges, index=[param]).to_csv(simulation_dir + '/strain_ranges.csv')
        np.savetxt(simulation_dir + '/param_' + param + '_inputs.csv', params, delimiter=',')
        pd.DataFrame(qois, columns=qoi_names).to_csv(
            simulation_dir + '/param_' + param + '_strain_qoi_outcomes.csv')

#######################################################################################################################
# Concatenate all results into a single CSV file
print('Showing slopes dataframe: ')
print(all_slopes)
print('Save all results into SA_summary_OAT files...')
all_slopes.to_csv('SA_summary_OAT_slopes.csv')
all_intercepts.to_csv('SA_summary_OAT_intercepts.csv')
all_p_values.to_csv('SA_summary_OAT_p_values.csv')
all_r_values.to_csv('SA_summary_OAT_r_values.csv')
all_r_values.to_csv('SA_summary_OAT_ranges.csv')
# print('Concatenate all results into SA_summary_OAT files...')
# pd.concat(all_slopes, axis='columns').to_csv('SA_summary_OAT_slopes.csv')
# pd.concat(all_intercepts, axis='columns').to_csv('SA_summary_OAT_intercepts.csv')
# pd.concat(all_p_values, axis='columns').to_csv('SA_summary_OAT_p_values.csv')
# pd.concat(all_r_values, axis='columns').to_csv('SA_summary_OAT_r_values.csv')
# pd.concat(all_ranges, axis='columns').to_csv('SA_summary_OAT_ranges.csv')
print('Finished.')
#######################################################################################################################
# # Postprocessing for visualisation purposes only
# for param in parameter_names:
#     simulation_dir = ''
#     if system == 'jureca':
#         simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
#     elif system == 'cosma':
#         simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
#     elif system == 'polaris':
#         simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
#     elif system == 'heart':
#         simulation_dir = sa_folder_root_name + '_' + param + '/'
#     parameter_names = np.array([param]) # Dummy
#     baseline_parameter_values = 0 # Dummy
#     baseline_dir = '' # Dummy
#     sa = SAUQ(name='sa', sampling_method='uniform', n=8, parameter_names=parameter_names,
#               baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#               simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
#     sa.run_jobs_postprocess(simulation_dir)
# #######################################################################################################################
# quit()
