from alyaformat import AlyaFormat
from postprocessing import PostProcessing
from sensitivityanalysis import SA
import os
import numpy as np

########################################################################################################################
# Global Settings
mesh_number = '05'
simulation_name = 'rodero_' + mesh_number + '_fine'
workdir = os.getcwd()
if 'icei' in workdir:
    system = 'jureca'
else:
    system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'

verbose = False

########################################################################################################################
# Step 1: Initialise Alya input files writing capabilities.
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_dir = './'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)

########################################################################################################################
# # Step 2: Use sampling methods to explore sensitivity analysis
# # Parameters to change: INaL, IKr, INaK, ICaL, Jrel, Jup, Ca50, kws, LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s,
# cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
# baseline_parameter_values = np.array([1,1,1,1,1,0.805,1])
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
# elif system == 'heart':
#     baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = sa_folder_name + '/'
# sa = SA(name='sa', sampling_method='saltelli', n=2 ** 4, parameter_names=cell_parameter_names,
#        baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#        simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup()
# # sa.run_jobs(simulation_dir+sa_folder_name +'/')
# # quit()

#######################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
sa_folder_name = 'sensitivity_analyses_cellular_parameters'
cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
baseline_parameter_values = np.array([1,1,1,1,1,0.805,1])
upper_bounds = baseline_parameter_values * 2.0
lower_bounds = baseline_parameter_values * 0.5
baseline_json_file = 'rodero_baseline_simulation_em.json'
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = sa_folder_name + '/'
# sa = SA(name='sa', sampling_method='saltelli', n=2 ** 4, parameter_names=cell_parameter_names,
#        baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#        simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# sa.run_jobs(simulation_dir)
# quit()

#######################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
# cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
sa_folder_name = 'sensitivity_analyses_mechanical_parameters'
haemodynamic_mechanical_parameter_names = np.array(['pericardial_stiffness', 'Kct_myocardium', 'a_myocardium', 'af_myocardium', 'afs_myocardium', 'cal50_myocardium', 'tref_scaling_myocardium'])
baseline_parameter_values = np.array([1200000,4000000,4500,15600,3000,0.805,15])
upper_bounds = baseline_parameter_values * 2.0
lower_bounds = baseline_parameter_values * 0.5
baseline_json_file = 'rodero_baseline_simulation_em.json'
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = sa_folder_name + '/'
sa = SA(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=haemodynamic_mechanical_parameter_names,
       baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
       simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# sa.run_jobs(simulation_dir)
# quit()
########################################################################################################################
# Step 3: Run Alya post-processing
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
elif system == 'heart':
    simulation_dir = '/users/jenang/Alya_setup_SA/alya_csv/rodero_' + mesh_number + '/'
sa.run_jobs_postprocess(simulation_dir)
quit()
########################################################################################################################
# Step 4: Evaluate QoIs and write out to results file
beat = 1
pv_post = sa.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
ecg_post = sa.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
fibre_work_post = sa.evaluate_qois(qoi_group_name='fibre_work', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
deformation_post = sa.evaluate_qois(qoi_group_name='deformation', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
sa.analyse(simulation_dir+'pv_qois.csv', qois = ['EDVL', 'LVEF', 'PmaxL', 'SVL'])
sa.analyse(simulation_dir+'deformation_qois.csv', qois=['max_basal_ab_displacement', 'min_basal_ab_displacement', 'max_apical_ab_displacement', 'min_apical_ab_displacement'])
quit()
########################################################################################################################
# Step 5: Evaluate Sobol indices and plot results
sa_figures_directory = simulation_dir
sa.analyse(sa_figures_directory, qois=['EDVL', 'LVEF', 'PmaxL', 'SVL'])