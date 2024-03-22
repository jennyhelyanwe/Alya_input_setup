from alyaformat import AlyaFormat
from postprocessing import PostProcessing
from sensitivityanalysis_uncertaintyquantification import SAUQ
from generatefields import FieldGeneration
import os
import numpy as np
import json

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
else:
    system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'cosma':
    meta_data_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False

########################################################################################################################
# Step 1: Initialise Alya input files writing capabilities.
if system == 'jureca':
    simulation_root_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'cosma':
    simulation_root_dir = '/snap8/scratch/dp287/dc-wang14/alya_simulations/'
elif system == 'polaris':
    simulation_root_dir = '/grand/projects/CompBioAffin/jenang/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_root_dir = './'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_root_dir, verbose=verbose)

########################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
conduction_parameter_names = np.array(['sigma_fibre', 'sigma_sheet', 'endocardial_activation_time_scaling'])
sa_folder_root_name = 'sensitivity_analyses_conduction_parameters_oat'
baseline_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
#######################################################################################################################
# Run simulations
for param in conduction_parameter_names:
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
    elif '_fibre' in param:
        baseline_parameter_values = np.array([simulation_dict[param.split('_fibre')[0]][0][0]])
    elif '_sheet' in param:
        baseline_parameter_values = np.array([simulation_dict[param.split('_sheet')[0]][0][1]])
    else:
        baseline_parameter_values = np.array([simulation_dict[param]])
    upper_bounds = baseline_parameter_values * 2.0
    lower_bounds = baseline_parameter_values * 0.1
    baseline_json_file = 'rodero_baseline_simulation_em.json'
    simulation_dir = ''
    baseline_dir = ''
    if system == 'jureca':
        baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'cosma':
        baseline_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'polaris':
        baseline_dir = '/grand/projects/CompBioAffin/jenang/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    elif system == 'heart':
        baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
        simulation_dir = sa_folder_root_name + '_' + param + '/'
    parameter_names = np.array([param])
    sa = SAUQ(name='sa', sampling_method='uniform', n=8, parameter_names=parameter_names,
              baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
              simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
    sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
    sa.run_jobs(simulation_dir)
quit()
## ########################################################################################################################
# Evaluate QoIs and correlations
for param in conduction_parameter_names:
    simulation_dir = ''
    if system == 'jureca':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'cosma':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'polaris':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'heart':
        simulation_dir = sa_folder_root_name + '_' + param + '/'
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
    sa = SAUQ(name='sa', sampling_method='uniform', n=8, parameter_names=parameter_names,
              baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
              simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
    labels = [param +'=' + s for s in ["%.1f" % x for x in np.linspace(lower_bounds, upper_bounds, 8)]]
    beat = 1
    ####################################################################################################################
    # Can be done as the simulations are running, using raw outputs.
    sa.sort_simulations(
        tag='raw')  # Collate list of finished simulations by checking the existence of particular files.
    # The lines below only need to be run once, it will write out the QoIs to the save directories as txt files.
    pv_post = sa.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                               analysis_type='sa')
    sa.visualise_sa(beat=1, pv_post=pv_post, labels=labels, save_filename=simulation_dir+'/pv_post.png')
    ecg_post = sa.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                analysis_type='sa')
    sa.visualise_sa(beat=1, ecg_post=ecg_post, labels=labels, save_filename=simulation_dir+'/ecg_post.png')
    # Analyse correlations
    qoi_names = ['EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection', 'dvdt_filling', 'dpdt_max', 'EDVR', 'ESVR',
                 'PmaxR', 'SVR']
    corrs, ranges = sa.analyse(filename=simulation_dir + 'pv_qois.csv', qois=qoi_names, show_healthy_ranges=False)
    pv_corrs = dict(map(lambda i, j: (i, j), qoi_names, corrs[:, 0]))
    pv_ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
    json.dump(pv_corrs, open(simulation_dir + '/pv_corrs.csv', 'w'))
    json.dump(pv_ranges, open(simulation_dir + '/pv_ranges.csv', 'w'))

    qoi_names = ['qt_dur_mean', 't_pe_mean', 't_peak_mean']
    corrs, ranges = sa.analyse(filename=simulation_dir + 'ecg_qois.csv', qois=qoi_names, show_healthy_ranges=False)
    ecg_corrs = dict(map(lambda i, j: (i, j), qoi_names, corrs[:, 0]))
    ecg_ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
    json.dump(ecg_corrs, open(simulation_dir + '/ecg_corrs.csv', 'w'))
    json.dump(ecg_ranges, open(simulation_dir + '/ecg_ranges.csv', 'w'))
quit()
########################################################################################################################
# Postprocessing
# for param in conduction_parameter_names:
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
########################################################################################################################


####################################################################################################################
# Needing proprocessing first
for param in conduction_parameter_names:
    simulation_dir = ''
    if system == 'jureca':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'cosma':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'polaris':
        simulation_dir = simulation_root_dir + sa_folder_root_name + '_' + param + '/'
    elif system == 'heart':
        simulation_dir = sa_folder_root_name + '_' + param + '/'
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
    sa = SAUQ(name='sa', sampling_method='uniform', n=8, parameter_names=parameter_names,
              baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
              simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
    labels = [param + '=' + s for s in ["%.1f" % x for x in np.linspace(lower_bounds, upper_bounds, 8)]]
    beat = 1
    # Visualise the simulations
    sa.sort_simulations(
        tag='postprocess')  # Collate list of finished simulations by checking the existence of particular files.
    # The lines below only need to be run once, it will write out the QoIs to the save directories as txt files.
    fibre_work_post = sa.evaluate_qois(qoi_group_name='fibre_work', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                       analysis_type='sa')
    sa.visualise_sa(beat=1, fibre_work_post=fibre_work_post, labels=labels,
                    save_filename=simulation_dir + '/fibre_work_post.png')
    deformation_post = sa.evaluate_qois(qoi_group_name='deformation', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                        analysis_type='sa')
    sa.visualise_sa(beat=1, deformation_post=deformation_post, labels=labels,
                    save_filename=simulation_dir + '/deformation_post.png')
    strain_post = sa.evaluate_qois(qoi_group_name='strain', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                   analysis_type='sa')
    sa.visualise_sa(beat=1, deformation_post=deformation_post, labels=labels,
                    save_filename=simulation_dir + '/strain_post.png')

    # Evaluate correlations
    qoi_names = ['es_ed_avpd', 'es_ed_apical_displacement', 'diff_lv_wall_thickness']
    corrs, ranges = sa.analyse(filename=simulation_dir + 'deformation_qois.csv', qois=qoi_names)
    deformation_corrs = dict(map(lambda i, j: (i, j), qoi_names, corrs[:, 0]))
    deformation_ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
    json.dump(deformation_corrs, open(simulation_dir + '/deformation_corrs.csv', 'w'))
    json.dump(deformation_ranges, open(simulation_dir + '/deformation_ranges.csv', 'w'))

    qoi_names = ['peak_lambda', 'min_lambda']
    corrs, ranges = sa.analyse(filename=simulation_dir + 'fibrework_qois.csv', qois=qoi_names)
    fibre_corrs = dict(map(lambda i, j: (i, j), qoi_names, corrs[:, 0]))
    fibre_ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
    json.dump(fibre_corrs, open(simulation_dir + '/fibre_corrs.csv', 'w'))
    json.dump(fibre_ranges, open(simulation_dir + '/fibre_ranges.csv', 'w'))

    qoi_names = ['max_median_mid_Ecc', 'max_median_mid_Err', 'max_median_four_chamber_Ell']
    corrs, ranges = sa.analyse(filename=simulation_dir + 'strain_qois.csv', qois=qoi_names)
    strain_corrs = dict(map(lambda i, j: (i, j), qoi_names, corrs[:, 0]))
    strain_ranges = dict(map(lambda i, j: (i, j), qoi_names, ranges[:, 0]))
    json.dump(strain_corrs, open(simulation_dir + '/strain_corrs.csv', 'w'))
    json.dump(strain_ranges, open(simulation_dir + '/strain_ranges.csv', 'w'))