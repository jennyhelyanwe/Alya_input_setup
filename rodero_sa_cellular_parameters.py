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
elif system == 'heart':
    simulation_root_dir = './'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_root_dir, verbose=verbose)


# Step 2: Use sampling methods to explore sensitivity analysis
sa_folder_name = 'sensitivity_analyses_cellular_parameters'
baseline_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
cellular_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup'])
baseline_parameter_values = np.array([simulation_dict['sf_gnal'][0][0],
                                      simulation_dict['sf_gkr'][0][0],
                                      simulation_dict['sf_gnak'][0][0],
                                      simulation_dict['sf_gcal'][0][0],
                                      simulation_dict['sf_jup'][0][0]])
upper_bounds = baseline_parameter_values * 2.0
lower_bounds = baseline_parameter_values * 0.5
baseline_json_file = 'rodero_baseline_simulation_em_for_sa.json'
if system == 'jureca':
    baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = simulation_root_dir + sa_folder_name + '/'
elif system == 'cosma':
    baseline_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = simulation_root_dir + sa_folder_name + '/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = sa_folder_name + '/'
sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 4, parameter_names=cellular_parameter_names,
          baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
          simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# quit()
sa.run_jobs(simulation_dir)
quit()
########################################################################################################################
# Step 3: Run Alya post-processing
if system == 'jureca':
    simulation_dir = simulation_root_dir + sa_folder_name + '/'
elif system == 'cosma':
    simulation_dir = simulation_root_dir + sa_folder_name + '/'
elif system == 'heart':
    simulation_dir = '/users/jenang/Alya_setup_SA/alya_csv/rodero_' + mesh_number + '/'
sa.run_jobs_postprocess(simulation_dir)
########################################################################################################################
# Step 4: Evaluate QoIs and write out to results file
beat = 1
sa.sort_simulations(tag='postprocess') # Collate list of finished simulations by checking the existence of particular files.
# The lines below only need to be run once, it will write out the QoIs to the save directories as txt files.
pv_post = sa.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
sa.visualise_sa(beat=1, pv_post=pv_post, labels=[])
ecg_post = sa.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
sa.visualise_sa(beat=1, ecg_post=ecg_post, labels=[])
fibre_work_post = sa.evaluate_qois(qoi_group_name='fibre_work', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
sa.visualise_sa(beat=1, fibre_work_post=fibre_work_post)
deformation_post = sa.evaluate_qois(qoi_group_name='deformation', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
sa.visualise_sa(beat=1, deformation_post=deformation_post, labels=[])

# Visualise scatter plots
sa.analyse(filename=simulation_dir+'pv_qois.csv', qois = ['EDVL', 'ESVL', 'PmaxL', 'LVEF', 'dvdt_ejection', 'dvdt_filling', 'dpdt_max'])
sa.analyse(filename=simulation_dir+'ecg_qois.csv', qois=['qt_dur_mean', 't_pe_mean', 't_peak_mean'])
sa.analyse(filename=simulation_dir+'deformation_qois.csv', qois=['es_ed_avpd', 'es_ed_apical_displacement'])
sa.analyse_qoi_vs_qoi(filename_1=simulation_dir+'pv_qois.csv', filename_2=simulation_dir+'ecg_qois.csv',
                      qois_1=['EDVL', 'ESVL', 'PmaxL', 'LVEF'], qois_2=['qt_dur_mean', 't_pe_mean'])

########################################################################################################################
# Step 5: Evaluate Sobol indices and plot results
sa_figures_directory = simulation_dir
sa.analyse(sa_figures_directory, qois=['EDVL', 'LVEF', 'PmaxL', 'SVL'])

########################################################################################################################
# Step 6: Illustrate sensitivity
