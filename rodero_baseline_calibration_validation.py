from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
# from calibratecv import CalibrateCV
from alyaformat import AlyaFormat, run_job, run_job_postprocess
from postprocessing import PostProcessing
import numpy as np
import os
import json
from sensitivityanalysis_uncertaintyquantification import SAUQ

########################################################################################################################
# Global Settings
mesh_number = '05'
simulation_name = 'rodero_' + mesh_number + '_fine'
workdir = os.getcwd()
if 'icei' in workdir:
    system = 'jureca'
elif 'cosma' in workdir:
    system = 'cosma'
elif 'e769' in workdir:
    system = 'archer2'
elif 'Expansion' in workdir:
    system = 'archive'
else:
    system = 'heart'

if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'cosma':
    meta_data_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/meta_data/'
elif system == 'archer2':
    meta_data_dir = '/work/e769/e769/jennywang/Alya_pipeline/meta_data/'
elif system == 'archive':
    meta_data_dir = '/run/media/jenang/Expansion/JURECA_COSMA_download_April2024/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
print(geometric_data_dir)
clinical_data_dir = meta_data_dir + 'clinical_data/'

verbose = False
mesh_preprocess = False
calibrate_cv = False
generate_fields_original_doste = False
generate_fields_Green_fibres = False
generate_fields_12090_fibres = False
generate_fields_slices_and_local_bases = False
setup_em_alya_literature_parameters_files = False
setup_em_alya_files = False
setup_ep_alya_files = False
run_alya_baseline_simulation = False
run_alya_baseline_postprocessing = False
decide_calibration_parameter_ranges = False
setup_calibration_alya_simulations = False
run_alya_calibration_simulations = True
evaluate_calibration_sa = False
setup_validation_alya_simulations = False
run_alya_validation_simulations = False
run_alya_validation_postprocessing = False
evaluate_validation_biomarkers = False

#######################################################################################################################
# Step 1: Save input mesh into CSV format, as prescribed in myformat.py
vtk_dir = ''
if system == 'heart':
    vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number + '/'
elif system == 'jureca':
    vtk_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/geometric_data/vtk/'
elif system == 'archer2':
    vtk_dir = meta_data_dir + '/geometric_data/vtk/'
vtk_name = mesh_number + '_bivent_only'
simulation_json_file = 'rodero_baseline_simulation_ep.json'
if mesh_preprocess:
    mesh = MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir, verbose=verbose)
    mesh.read_geometry_from_vtk_rodero(save=False)
    mesh.generate_boundary_data_rodero(save=True)
    mesh.check_fields_for_qrs_inference()
    mesh.check_fields_for_twave_inference()

########################################################################################################################
# Step 2: Run QRS and T inference and write personalised results to personalisation_data_dir

########################################################################################################################
# Step 3: Calibrate conductivities to personalisation conduction velocity
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
# if calibrate_cv:
#     calibration_dir = meta_data_dir + 'calibration_dir/'
#     simulation_json_file = 'rodero_baseline_simulation_ep.json'
#     calibrate = CalibrateCV(name=simulation_name, geometric_data_dir=geometric_data_dir,
#                             calibration_dir=calibration_dir, simulation_json_file=simulation_json_file,
#                             personalisation_data_dir=personalisation_data_dir, verbose=verbose)

########################################################################################################################
# Step 4: Generate fields for Alya simulation
electrode_data_filename = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_electrode_xyz.csv'
fields = FieldGeneration(name=simulation_name, geometric_data_dir=geometric_data_dir,
                personalisation_data_dir=personalisation_data_dir, verbose=verbose)

if generate_fields_slices_and_local_bases:
    fields.generate_orthogonal_local_bases(save=False)
    fields.generate_short_long_axes_slices(save=True)

doste_fibre_directory = geometric_data_dir + 'CR05_orig_files/'
if generate_fields_original_doste:
    map = fields.map_doste_nodes_to_rodero_nodes(fibre_vtk_filename=doste_fibre_directory+'Long_Fibers0_0_20_CR05_orig.vtk',
                                                 map_filename=doste_fibre_directory+'doste_to_rodero_node_map.txt')
    fields.read_doste_fibre_fields_vtk(fibre_vtk_filename=doste_fibre_directory+'Long_Fibers0_0_20_CR05_orig.vtk',
                                       sheet_vtk_filename=doste_fibre_directory+'Sheet_Fibers0_0_20_CR05_orig.vtk',
                                       normal_vtk_filename=doste_fibre_directory+'Normal_Fibers0_0_20_CR05_orig.vtk',
                                       save=True, map=map)
if generate_fields_Green_fibres:
    map = fields.map_doste_nodes_to_rodero_nodes(
        fibre_vtk_filename=doste_fibre_directory + 'Long_Green_Fibers0_0_0_CR05_orig.vtk',
        map_filename=doste_fibre_directory + 'doste_to_rodero_node_map.txt')
    fields.read_doste_fibre_fields_vtk(fibre_vtk_filename=doste_fibre_directory+'Long_Green_Fibers0_0_0_CR05_orig.vtk',
                                       sheet_vtk_filename=doste_fibre_directory+'Sheet_Green_Fibers0_0_0_CR05_orig.vtk',
                                       normal_vtk_filename=doste_fibre_directory+'Normal_Green_Fibers0_0_0_CR05_orig.vtk',
                                       save=True, map=map)

if generate_fields_12090_fibres:
    map = fields.map_doste_nodes_to_rodero_nodes(
        fibre_vtk_filename=doste_fibre_directory + 'Long_Green_Fibers0_0_0_CR05_orig.vtk',
        map_filename=doste_fibre_directory + 'doste_to_rodero_node_map.txt')
    fields.read_doste_fibre_fields_vtk(fibre_vtk_filename=doste_fibre_directory+'Long_12090_Fibers0_0_0_CR05_orig.vtk',
                                       sheet_vtk_filename=doste_fibre_directory+'Sheet_12090_Fibers0_0_0_CR05_orig.vtk',
                                       normal_vtk_filename=doste_fibre_directory+'Normal_12090_Fibers0_0_0_CR05_orig.vtk',
                                       save=True, map=map)

########################################################################################################################
# Step 5: Write Alya input files according to simulation protocol saved in .json file.
simulation_dir = ''
if system == 'jureca':
    simulation_root_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'cosma':
    simulation_root_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_root_dir = './'
elif system == 'archer2':
    simulation_root_dir = '/work/e769/e769/jennywang/Alya_pipeline/alya_simulations/'
elif system == 'archive':
    simulation_root_dir = '/run/media/jenang/Expansion/JURECA_COSMA_download_April2024/Alya_pipeline/alya_simulations/'

alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_root_dir, job_version=system, verbose=verbose)

# Sanity check:
# if not system == 'jureca' and not system == 'cosma' and not system == 'archer2':
#     alya.visual_sanity_check(simulation_json_file=simulation_json_file)
if setup_em_alya_files:
    simulation_json_file = 'rodero_baseline_simulation_em.json'
    alya.do(simulation_json_file=simulation_json_file)
if setup_ep_alya_files:
    simulation_json_file = 'rodero_baseline_simulation_ep.json'
    alya.do(simulation_json_file=simulation_json_file)

if setup_em_alya_literature_parameters_files:
    # simulation_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
    simulation_json_file = 'rodero_baseline_simulation_em_literature_parameters_high_sfkws.json'
    alya.do(simulation_json_file=simulation_json_file)

########################################################################################################################
# Step 5: Run baseline simulation
if run_alya_baseline_simulation:
    run_job(alya.output_dir)
if run_alya_baseline_postprocessing:
    run_job_postprocess(alya.output_dir)

########################################################################################################################
# Step 6: Postprocess
if decide_calibration_parameter_ranges:
    print('Evaluating simulated biomarkers')
    simulation_json_file = 'rodero_baseline_simulation_em.json'
    alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[0] + '_literature_parameters_' + simulation_name + '/'
    # alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[
    #     0] + '_' + simulation_name + '_mec_baseline/'
    pp = PostProcessing(alya=alya, simulation_json_file=simulation_json_file,
                        alya_output_dir=alya_output_dir, protocol='', verbose=verbose)
    beat = 1
    pp.evaluate_pv_biomarkers(beat=beat)
    pp.evaluate_ecg_biomarkers(beat=beat, show_landmarks=False)
    # pp.evaluate_deformation_biomarkers(beat=beat)
    # pp.evaluate_fibre_work_biomarkers(beat=beat)
    # pp.evaluate_strain_biomarkers(beat=beat)
    # pp.visualise_qoi_comparisons(qoi_names = ['qrs_dur_mean', 'qt_dur_mean', 't_pe_mean', 'EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection',
    #                  'dvdt_filling', 'dpdt_max'], save_figure=alya_output_dir+'/qoi_evaluation.png')
    # Rank QoIs according to importance for matching - based on ??

    # Use OAT SA results and evaluation of biomarkers to assign new ranges for calibration SA
    ranked_qoi_names = ['LVEF', 'PmaxL', 'SVL', 'EDVL', 'ESVL', 'dpdt_max',
                        'dvdt_ejection', 'dvdt_filling']
    baseline_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
    simulation_json_file = baseline_json_file
    simulation_dict = json.load(open(simulation_json_file, 'r'))
    pp.calculate_calibration_sa_parameter_ranges(ranked_qoi_names=ranked_qoi_names,
                                                 oat_sa_slopes='SA_summary_OAT_slopes.csv',
                                                 oat_sa_p_values='SA_summary_OAT_p_values.csv',
                                                 oat_sa_r_values='SA_summary_OAT_r_values.csv',
                                                 oat_sa_ranges='SA_summary_OAT_ranges.csv',
                                                 oat_sa_intercepts='SA_summary_OAT_intercepts.csv',
                                                 simulation_dict=simulation_dict, strategy='one_qoi_at_a_time', qoi_input='LVEF')
    # pp.visualise_calibration_sa_parameter_ranges(oat_sa_slopes='SA_summary_OAT_slopes.csv',
    #                                              oat_sa_p_values='SA_summary_OAT_p_values.csv',
    #                                              oat_sa_r_values='SA_summary_OAT_r_values.csv',
    #                                              oat_sa_ranges='SA_summary_OAT_ranges.csv',
    #                                              oat_sa_intercepts='SA_summary_OAT_intercepts.csv',
    #                                              simulation_root_dir=simulation_root_dir)
    quit()
# Using sfkws and Kct ranges of {"sfkws_myocardium": [2.6391168230682993, 4.0779357886520335], "Kct_myocardium": [500000.0, 5000000.0]}
# Gave an averaged LVEF of ~30%, but an averaged peak LVP of ~22.5 kPa. The rational thing to do now is to reduce the resistance
# to get within healthy ranges.
# Compliance ranges:
# upper_range = (143000 - 225000 -39524877.19 * 0.000150) / -39524877.19 = (-82000 -  5928)/-39524877.19 = 0.0022
# lower_range = (153000 - 225000 - 5928 ) / -39524877.19 = 77928/  39524877.19 = 0.00197

# Third iteration. LVEF now sitting at around 42 %, PmaxL around 13.2 kPa, which is a bit too low, and the ejection phase is too linear.
# We need one more parameter to push the LVEF higher, and allow PmaxL to go back up.
# The next parameter with highest range in LVEF is ICaL. So, time to ramp that up to the 1.8 to 2.0 range.
# lower_range = 1.8
# upper_range = 2.0
# Adjust compliance:
# [ 0.0016, 0.002]

# Fifth iteration:
# "sfkws_myocardium": [4, 7], "cal50_myocardium": [0.5, 0.7], "Kct_myocardium": [50000.0,  500000.0], "arterial_resistance_lv":  [600, 999]}
# Max LVEF: 46.37
# Parameters:
# sfkws: 4.468
# cal50: 0.518
# Kct: 429687.5
# R_LV: 712.218

# Sixth iteration:
# Set: sfkws: 4.5
# Set: cal50: 0.52
# Set: Kct: 450,000
# Search: R_LV: [600, 999], ejection_threshold: [7, 10] kPa
#

# Max LVEF: 42%,
# Seventh iteration: add compliance
# Set: sfkws: 5, set cal50: 0.5, set Kct: 450,000 - use same baseline as the sixth iteration
# Search: R_LV: [700, 999], ejection threshold: [7, 10], C_LV: [0.000065, 0.00015]


########################################################################################################################
# Step 7: Use OAT SA results and evaluation of biomarkers to assign new ranges for calibration SA
# calibration_folder_name = 'calibration_simulations_third_iteration'
calibration_folder_name = 'calibration_simulations_seventh_iteration'
baseline_json_file = 'rodero_baseline_simulation_baseline_sixth_iteration.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
# perturbed_parameters = json.load(open('calibration_sa_ranges_third_iteration.json', 'r'))
perturbed_parameters = json.load(open('calibration_sa_ranges_seventh_iteration.json', 'r'))
perturbed_parameters_name = np.array(list(perturbed_parameters.keys()))
print(perturbed_parameters_name)
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + calibration_folder_name + '/'
elif system == 'cosma':
    baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = simulation_root_dir + calibration_folder_name + '/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = calibration_folder_name + '/'
elif system == 'archer2':
    baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_literature_parameters_rodero_05_fine/'
    simulation_dir = simulation_root_dir + calibration_folder_name + '/'
upper_bounds = []
lower_bounds = []
for param in perturbed_parameters_name:
    lower_bounds.append(perturbed_parameters[param][0])
    upper_bounds.append(perturbed_parameters[param][1])
calibration = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 4 , parameter_names=perturbed_parameters_name,
                   baseline_json_file=baseline_json_file, simulation_dir=simulation_dir, alya_format=alya,
                   baseline_dir=baseline_dir, verbose=verbose)
if setup_calibration_alya_simulations:
    calibration.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
if run_alya_calibration_simulations:
    # calibration.run_jobs(simulation_dir, start_id=0) # Maximum job submission in archer2 is 128 for QoS:taskfarm.
    calibration.run_jobs(simulation_dir, start_id=55)
# calibration.run_jobs(simulation_dir, start_id=0, end_id=128) # Maximum job submission in archer2 is 128 for QoS:taskfarm.
if evaluate_calibration_sa:
    if system == 'archive':
        calibration.sort_simulations_archive(tag='raw')
    else:
        calibration.sort_simulations(
            tag='raw')  # Collate list of finished simulations by checking the existence of particular files.
    # calibration.visualise_finished_parameter_sets(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
    evaluate_pv = True
    evaluate_ecg = False
    evaluate_deformation = False
    evaluate_fibrework = False
    evaluate_strain = False
    beat = 1
    if evaluate_pv:
        pv_post = calibration.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                             analysis_type='sa')
        calibration.visualise_sa(beat=1, pv_post=pv_post, save_filename=simulation_dir + '/pv_post.png', show=True)
        # calibration.visualise_sa(beat=1, pv_post=pv_post, show=True)
        # qoi_names = ['EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection', 'dvdt_filling', 'dpdt_max', 'EDVR',
        #              'ESVR', 'PmaxR', 'SVR']
        qoi_names = ['LVEF', 'ESVL', 'PmaxL', 'SVL', 'dvdt_ejection', 'dvdt_filling', 'dpdt_max']
        # slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
        #     filename=simulation_dir + 'pv_qois.csv', qois=qoi_names, show_healthy_ranges=False,
        #     save_filename=simulation_dir + '/pv_scatter.png')
        slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
            filename=simulation_dir + 'pv_qois.csv', qois=qoi_names, show_healthy_ranges=False)
    elif evaluate_ecg:
        ecg_post = calibration.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                             analysis_type='sa')
        calibration.visualise_sa(beat=1, pv_post=ecg_post, save_filename=simulation_dir + '/ecg_post.png')
        qoi_names = ['qrs_dur_mean', 't_dur_mean', 'qt_dur_mean', 't_pe_mean', 'jt_dur_mean']
        slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
            filename=simulation_dir + 'ecg_qois.csv', qois=qoi_names, show_healthy_ranges=False,
            save_filename=simulation_dir + '/ecg_scatter.png')
    elif evaluate_deformation:
        deformation_post = calibration.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                             analysis_type='sa')
        calibration.visualise_sa(beat=1, pv_post=deformation_post, save_filename=simulation_dir + '/ecg_post.png')
        qoi_names = ['es_ed_avpd', 'es_ed_apical_displacement', 'diff_lv_wall_thickness']
        slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
            filename=simulation_dir + 'deformation_qois.csv', qois=qoi_names, show_healthy_ranges=False,
            save_filename=simulation_dir + '/deformation_scatter.png')

    elif evaluate_strain:
        strain_post = calibration.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat,
                                                     qoi_save_dir=simulation_dir,
                                                     analysis_type='sa')
        calibration.visualise_sa(beat=1, pv_post=strain_post, save_filename=simulation_dir + '/ecg_post.png')
        qoi_names = ['max_mid_Ecc', 'min_mid_Ecc', 'max_mid_Err', 'min_mid_Err', 'max_four_chamber_Ell',
                     'min_four_chamber_Ell']
        slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
            filename=simulation_dir + 'strain_qois.csv', qois=qoi_names, show_healthy_ranges=False,
            save_filename=simulation_dir + '/strain_scatter.png')

quit()
########################################################################################################################
# Step 6: Validation experiments - volume perturbation
validation_folder_name = 'volume_perturbation_validation_experiments'
baseline_json_file = 'rodero_baseline_simulation_em.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
perturbed_parameters_name = np.array(['end_diastole_p_lv'])
baseline_parameter_values = np.array([simulation_dict['end_diastole_p'][0]])
upper_bounds = baseline_parameter_values * 4.0
lower_bounds = baseline_parameter_values * 0.8
baseline_dir = ''
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + validation_folder_name + '/'
elif system == 'cosma':
    baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_rodero_05_fine_mec_baseline/'
    simulation_dir = simulation_root_dir + validation_folder_name + '/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = validation_folder_name + '/'
elif system == 'archer2':
    baseline_dir = simulation_root_dir + 'rodero_baseline_simulation_em_rodero_05_fine_mec_baseline/'
    simulation_dir = simulation_root_dir + validation_folder_name + '/'
experiment = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=perturbed_parameters_name,
                  baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
                  simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
if setup_validation_alya_simulations:
    experiment.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
if run_alya_validation_simulations:
    experiment.run_jobs(simulation_dir)
if run_alya_validation_postprocessing:
    experiment.run_jobs_postprocess(simulation_dir)
if evaluate_validation_biomarkers:
    beat = 1
    experiment.sort_simulations(
        tag='raw')  # Collate list of finished simulations by checking the existence of particular files.
    pv_post = experiment.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir, analysis_type='sa')
    ed_es_pvr_post = experiment.visualise_ed_es_pvr_biomarkers(beat=beat, pv_post=pv_post)
    ecg_post = experiment.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                analysis_type='sa')
    experiment.visualise_sa(beat=1, ecg_post=ecg_post)