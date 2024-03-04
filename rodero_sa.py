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
# sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 4, parameter_names=cell_parameter_names,
#        baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#        simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup()
# # sa.run_jobs(simulation_dir+sa_folder_name +'/')
# # quit()

#######################################################################################################################
# # Step 2: Use sampling methods to explore sensitivity analysis
# # Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
# sa_folder_name = 'sensitivity_analyses_cellular_parameters'
# cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
# baseline_parameter_values = np.array([1,1,1,1,1,0.805,1])
# upper_bounds = baseline_parameter_values * 2.0
# lower_bounds = baseline_parameter_values * 0.5
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# baseline_dir = ''
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
# elif system == 'heart':
#     baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = sa_folder_name + '/'
# # sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 4, parameter_names=cell_parameter_names,
# #        baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
# #        simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# # sa.run_jobs(simulation_dir)
# # quit()

######################################################################################################################
# # Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
# # cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
# sa_folder_name = 'sensitivity_analyses_mechanical_parameters'
# mechanical_parameter_names = np.array(['pericardial_stiffness', 'Kct_myocardium', 'a_myocardium', 'af_myocardium', 'afs_myocardium', 'cal50_myocardium', 'tref_scaling_myocardium'])
# baseline_parameter_values = np.array([1200000,4000000,4500,15600,3000,0.805,15])
# upper_bounds = baseline_parameter_values * 2.0
# lower_bounds = baseline_parameter_values * 0.5
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
# elif system == 'heart':
#     baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = sa_folder_name + '/'
# sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=mechanical_parameter_names,
#           baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#           simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# # sa.run_jobs(simulation_dir)
# # quit()
########################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
# cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
sa_folder_name = 'sensitivity_analyses_literature_baseline_mechanical_parameters'
baseline_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
mechanical_parameter_names = np.array(['pericardial_stiffness', 'Kct_myocardium', 'a_myocardium', 'af_myocardium', 'afs_myocardium', 'tref_scaling_myocardium'])
baseline_parameter_values = np.array([simulation_dict['pericardial_stiffness'],
                                      simulation_dict['Kct'][0],
                                      simulation_dict['a'][0],
                                      simulation_dict['af'][0],
                                      simulation_dict['afs'][0],
                                      simulation_dict['tref_scaling'][0]])
upper_bounds = baseline_parameter_values * 2.0
lower_bounds = baseline_parameter_values * 0.5
baseline_json_file = 'rodero_baseline_simulation_em.json'
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = sa_folder_name + '/'
sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 3, parameter_names=mechanical_parameter_names,
          baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
          simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# sa.run_jobs(simulation_dir)
# quit()

######################################################################################################################
# # Step 2: Explore fibre angles
# fields = FieldGeneration(name=simulation_name, geometric_data_dir=geometric_data_dir,
#                 personalisation_data_dir=personalisation_data_dir, verbose=verbose)
# fibre_sa_folder_name = 'sensitivity_analyses_transmural_fibres'
# sa_folder_name = fibre_sa_folder_name
# parameter_names = ''
# baseline_parameter_values = []
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# baseline_dir = ''
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + fibre_sa_folder_name + '/'
# elif system == 'heart':
#     baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = fibre_sa_folder_name + '/'
# doste_fibre_directory = geometric_data_dir + 'CR05_orig_files/'
# fibre_filenames = [doste_fibre_directory+'Long_Fibers0_0_20_CR05_orig.vtk', doste_fibre_directory+'Long_Green_Fibers0_0_0_CR05_orig.vtk',
#                doste_fibre_directory+'Long_Green_inv_Fibers0_0_0_CR05_orig.vtk', doste_fibre_directory+'Long_12090_Fibers0_0_0_CR05_orig.vtk',
#                doste_fibre_directory+'Long_9090_Fibers0_0_0_CR05_orig.vtk']
# sheet_filenames = [doste_fibre_directory+'Sheet_Fibers0_0_20_CR05_orig.vtk', doste_fibre_directory+'Sheet_Green_Fibers0_0_0_CR05_orig.vtk',
#                doste_fibre_directory+'Sheet_Green_inv_Fibers0_0_0_CR05_orig.vtk', doste_fibre_directory+'Sheet_12090_Fibers0_0_0_CR05_orig.vtk',
#                doste_fibre_directory+'Sheet_9090_Fibers0_0_0_CR05_orig.vtk']
# normal_filenames = [doste_fibre_directory+'Normal_Fibers0_0_20_CR05_orig.vtk', doste_fibre_directory+'Normal_Green_Fibers0_0_0_CR05_orig.vtk',
#                doste_fibre_directory+'Normal_Green_inv_Fibers0_0_0_CR05_orig.vtk', doste_fibre_directory+'Normal_12090_Fibers0_0_0_CR05_orig.vtk',
#                doste_fibre_directory+'Normal_9090_Fibers0_0_0_CR05_orig.vtk']
# map_filename = doste_fibre_directory+'doste_to_rodero_node_map.txt'
# sa = SAUQ(name='fibre_sa', sampling_method='fibre_test', n=0, parameter_names=np.array(['']),
#           baseline_parameter_values=np.array([]), baseline_json_file=baseline_json_file,
#           simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup_fibre_sa(fields=fields, fibre_filenames=fibre_filenames, sheet_filenames=sheet_filenames,
# #                   normal_filenames=normal_filenames, baseline_json_file=baseline_json_file, map_filename=map_filename)
# # quit()
# # sa.run_jobs(simulation_dir=simulation_dir, start_id=0)
# # quit()
#
# # sa.run_jobs_postprocess(simulation_dir=simulation_dir)
# # quit()
######################################################################################################################
# # Step 2: Use sampling methods to explore sensitivity analysis
# sa_folder_name = 'sensitivity_analyses_avpd'
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# simulation_json_file = baseline_json_file
# simulation_dict = json.load(open(simulation_json_file, 'r'))
# mechanical_parameter_names = np.array(['tref_normal_scaling_lv'])
# baseline_parameter_values = np.array([simulation_dict['tref_normal_scaling'][0]])
# upper_bounds = baseline_parameter_values * 2.0
# lower_bounds = baseline_parameter_values * 0.5
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
# elif system == 'heart':
#     baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = sa_folder_name + '/'
# sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=mechanical_parameter_names,
#           baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#           simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# # sa.run_jobs(simulation_dir)
# # quit()
#######################################################################################################################
# # Step 2: Use sampling methods to explore sensitivity analysis
# sa_folder_name = 'sensitivity_analyses_haemodynamic_parameters'
# # sa_folder_name = 'sensitivity_analyses_haemodynamic_parameters_v2'
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# simulation_json_file = baseline_json_file
# simulation_dict = json.load(open(simulation_json_file, 'r'))
# # First attempt; trying to get PmaxL in range, and dvdt in range for filling and relaxation
# haemodynamic_parameter_names = np.array(['arterial_compliance_lv',
#                                          'arterial_resistance_lv',
#                                          'gain_derror_relaxation_lv',
#                                          'gain_error_relaxation_lv',
#                                          'ejection_pressure_threshold_lv'])
# baseline_parameter_values = np.array([simulation_dict['arterial_compliance'][0],
#                                       simulation_dict['arterial_resistance'][0],
#                                       simulation_dict['gain_derror_relaxation'][0],
#                                       simulation_dict['gain_error_relaxation'][0],
#                                       simulation_dict['ejection_pressure_threshold'][0]])
# # Second attempt: narrow down parameters that give good PmaxL, trying to get slower dvdt ejection, and faster dvdt filling
# # haemodynamic_parameter_names = np.array(['gain_error_contraction_lv',
# #                                          'gain_derror_contraction_lv'])
# # baseline_parameter_values = np.array([simulation_dict['gain_error_contraction'][0],
# #                                       simulation_dict['gain_derror_contraction'][0]])
#
# upper_bounds = baseline_parameter_values * 2.0
# lower_bounds = baseline_parameter_values * 0.5
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + sa_folder_name + '/'
# elif system == 'heart':
#     baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = sa_folder_name + '/'
# sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=haemodynamic_parameter_names,
#           baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#           simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# # sa.run_jobs(simulation_dir)
# # quit()
########################################################################################################################
# Step 3: Run Alya post-processing
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + fibre_sa_folder_name + '/'
elif system == 'heart':
    simulation_dir = '/users/jenang/Alya_setup_SA/alya_csv/rodero_' + mesh_number + '/'
# sa.run_jobs_postprocess(simulation_dir)
# quit()
########################################################################################################################
# Step 4: Evaluate QoIs and write out to results file
beat = 1
sa.sort_simulations(tag='postprocess') # Collate list of finished simulations by checking the existence of particular files.
# The lines below only need to be run once, it will write out the QoIs to the save directories as txt files.
pv_post = sa.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='fibre_sa')
sa.visualise_sa(beat=1, pv_post=pv_post, labels=['Original Doste', 'Greenbaum','Inverted Greenbaum', '12090','9090'])
ecg_post = sa.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='fibre_sa')
sa.visualise_sa(beat=1, ecg_post=ecg_post, labels=['Original Doste', 'Greenbaum','Inverted Greenbaum', '12090','9090'])
# fibre_work_post = sa.evaluate_qois(qoi_group_name='fibre_work', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
# sa.visualise_sa(beat=1, fibre_work_post=fibre_work_post)
deformation_post = sa.evaluate_qois(qoi_group_name='deformation', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='fibre_sa')
sa.visualise_sa(beat=1, deformation_post=deformation_post, labels=['Original Doste', 'Greenbaum','Inverted Greenbaum', '12090','9090'])

# Visualise scatter plots
# sa.analyse(filename=simulation_dir+'pv_qois.csv', qois = ['EDVL', 'ESVL', 'PmaxL', 'LVEF', 'dvdt_ejection', 'dvdt_filling', 'dpdt_max'])
# sa.analyse(filename=simulation_dir+'ecg_qois.csv', qois=['qt_dur_mean', 't_pe_mean', 't_peak_mean'])
# sa.analyse(filename=simulation_dir+'deformation_qois.csv', qois=['es_ed_avpd', 'es_ed_apical_displacement'])
# sa.analyse_qoi_vs_qoi(filename_1=simulation_dir+'pv_qois.csv', filename_2=simulation_dir+'ecg_qois.csv',
#                       qois_1=['EDVL', 'ESVL', 'PmaxL', 'LVEF'], qois_2=['qt_dur_mean', 't_pe_mean'])
quit()
########################################################################################################################
# Step 5: Evaluate Sobol indices and plot results
sa_figures_directory = simulation_dir
sa.analyse(sa_figures_directory, qois=['EDVL', 'LVEF', 'PmaxL', 'SVL'])

