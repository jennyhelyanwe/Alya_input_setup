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
else:
    system = 'heart'

if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'cosma':
    meta_data_dir = '/cosma8/data/dp287/dc-wang14/Alya_pipeline/meta_data/'
elif system == 'archer2':
    meta_data_dir = '/work/e769/e769/jennywang/Alya_pipeline/meta_data/'
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
evaluate_simulated_biomarkers = True
calculate_sa_range_based_on_oat_results = True
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
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_root_dir, job_version=system, verbose=verbose)

# Sanity check:
if not system == 'jureca' and not system == 'cosma' and not system == 'archer2':
    alya.visual_sanity_check(simulation_json_file=simulation_json_file)
if setup_em_alya_files:
    simulation_json_file = 'rodero_baseline_simulation_em.json'
    alya.do(simulation_json_file=simulation_json_file)
if setup_ep_alya_files:
    simulation_json_file = 'rodero_baseline_simulation_ep.json'
    alya.do(simulation_json_file=simulation_json_file)

if setup_em_alya_literature_parameters_files:
    simulation_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
    alya.do(simulation_json_file=simulation_json_file)

########################################################################################################################
# Step 5: Run baseline simulation
if run_alya_baseline_simulation:
    run_job(alya.output_dir)
if run_alya_baseline_postprocessing:
    run_job_postprocess(alya.output_dir)

########################################################################################################################
# Step 6: Postprocess
if evaluate_simulated_biomarkers:
    print('Evaluating simulated biomarkers')
    simulation_json_file = 'rodero_baseline_simulation_em.json'
    alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[0] + '_literature_parameters_' + simulation_name + '/'
    # alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[
    #     0] + '_' + simulation_name + '_mec_baseline/'
    pp = PostProcessing(alya=alya, simulation_json_file=simulation_json_file,
                        alya_output_dir=alya_output_dir, protocol='postprocess', verbose=verbose)
    beat = 1
    pp.evaluate_pv_biomarkers(beat=beat)
    pp.evaluate_ecg_biomarkers(beat=beat, show_landmarks=False)
    # pp.evaluate_deformation_biomarkers(beat=beat)
    # pp.evaluate_fibre_work_biomarkers(beat=beat)
    # pp.evaluate_strain_biomarkers(beat=beat)
    # pp.visualise_qoi_comparisons(qoi_names = ['qrs_dur_mean', 'qt_dur_mean', 't_pe_mean', 'EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection',
    #                  'dvdt_filling', 'dpdt_max'], save_figure=alya_output_dir+'/qoi_evaluation.png')
    # Rank QoIs according to importance for matching - based on ??
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
                                                 simulation_dict=simulation_dict)

    # pp.visualise_calibration_comparisons_global(beat=beat)
    # pp.visualise_calibration_comparisons_strain()
    # pp.compare_ecg_with_clinical_ranges(beat=beat)
    # pp.compare_deformation_with_clinical_ranges(beat=beat)

########################################################################################################################
# Step 7: Use OAT SA results and evaluation of biomarkers to assign new ranges for calibration SA
if calculate_sa_range_based_on_oat_results:
    print ('Calculating new SA ranges based on OAT results and evaluation of previous baseline')

########################################################################################################################
# Step 6: Validation experiments - volume perturbation
validation_folder_name = 'volume_perturbation_validation_experiments'
baseline_json_file = 'rodero_baseline_simulation_em.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
perturbed_parameter_name = np.array(['end_diastole_p_lv'])
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
experiment = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=perturbed_parameter_name,
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