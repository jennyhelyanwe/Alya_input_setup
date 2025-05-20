from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
# from calibratecv import CalibrateCV
from alyaformat import AlyaFormat, run_job, run_job_postprocess
from postprocessing import PostProcessing
import numpy as np
import os
import json

########################################################################################################################
# Global Settings
vtu_names = ['1032146', '1037406', '1347962']
# vtk_name = '1097970'
vtu_name = vtu_names[0] + '_combo_fibres'
mesh_number = vtu_names[0]
simulation_name = mesh_number
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

max_cores_used = 3

if system == 'archer2':
    meta_data_dir = '/work/e769/e769/jennywang/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/UKB_'+mesh_number+'_HSmith/'
print('Geometric data directory: ', geometric_data_dir)
clinical_data_dir = meta_data_dir + 'clinical_data/'

verbose = False
mesh_preprocess = True
calibrate_cv = False
generate_fields_original_doste = False
generate_fields_Green_fibres = False
generate_fields_12090_fibres = False
generate_fields_slices_and_local_bases = False
setup_em_alya_literature_parameters_files = False
setup_em_alya_files = True
setup_ep_alya_files = False
run_alya_baseline_simulation = False
run_alya_baseline_postprocessing = False
evaluate_calibrated_baseline = False
setup_calibration_alya_simulations = False
run_alya_calibration_simulations = False
evaluate_calibration_sa = False
setup_validation_alya_simulations = False
run_alya_validation_simulations = False
run_alya_validation_postprocessing = False
evaluate_validation_biomarkers = False

#######################################################################################################################
# Step 1: Save input mesh into CSV format, as prescribed in myformat.py
vtk_dir = geometric_data_dir
# Meshes from Hannah Smith 20 May 2025

if mesh_preprocess:
    mesh = MeshPreprocessing(vtk_name=vtu_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
                             max_cores_used=max_cores_used, verbose=verbose)
    # mesh.read_geometry_from_vtu_hannah_smith(vtu_name=vtu_name,save=True)
    # mesh.add_cap_to_mesh()

    # mesh.generate_boundary_data_UKB(save=True)
    # mesh.check_fields_for_qrs_inference()
    # mesh.check_fields_for_twave_inference()

########################################################################################################################
# Step 2: Run QRS and T inference and write personalised results to personalisation_data_dir

########################################################################################################################
# Step 3: Calibrate conductivities to personalisation conduction velocity
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/UKB_'+mesh_number+'/'
# if calibrate_cv:
#     calibration_dir = meta_data_dir + 'calibration_dir/'
#     simulation_json_file = 'UKB_baseline_simulation_ep.json'
#     calibrate = CalibrateCV(name=simulation_name, geometric_data_dir=geometric_data_dir,
#                             calibration_dir=calibration_dir, simulation_json_file=simulation_json_file,
#                             personalisation_data_dir=personalisation_data_dir, verbose=verbose)

########################################################################################################################
# Step 4: Generate fields for Alya simulation
electrode_data_filename = meta_data_dir + 'geometric_data/UKB_'+mesh_number+'/UKB_'+mesh_number+'_electrode_xyz.csv'
fields = FieldGeneration(name=simulation_name, geometric_data_dir=geometric_data_dir,
                personalisation_data_dir=personalisation_data_dir, max_cores_used=max_cores_used, verbose=verbose)

fields.generate_infarct_borderzone_Wallman(visualise=True)
quit()
# fields.generate_infarct_borderzone()
# fields.generate_celltype(save=True)
#
# if generate_fields_slices_and_local_bases:
#     fields.generate_orthogonal_local_bases(save=False)
#     fields.generate_short_long_axes_slices(save=True)

########################################################################################################################
# Step 5: Write Alya input files according to simulation protocol saved in .json file.
simulation_root_dir = ''
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
                  simulation_dir = simulation_root_dir, job_version=system, max_cores_used=max_cores_used,
                  verbose=verbose)

# Sanity check:
# if not system == 'jureca' and not system == 'cosma' and not system == 'archer2':
#     alya.visual_sanity_check(simulation_json_file=simulation_json_file)
if setup_em_alya_files:
    simulation_json_file = 'UKB_'+mesh_number + '_healthy_baseline_simulation_em.json'
    alya.do(simulation_json_file=simulation_json_file)
quit()
if setup_em_alya_literature_parameters_files:
    simulation_json_file = 'UKB_baseline_simulation_em_literature_parameters.json'
    alya.do(simulation_json_file=simulation_json_file)

########################################################################################################################
# Step 5: Run baseline simulation
if run_alya_baseline_simulation:
    run_job(alya.output_dir)
if run_alya_baseline_postprocessing:
    run_job_postprocess(alya.output_dir)

########################################################################################################################
# Step 6: Postprocess
if evaluate_calibrated_baseline:
    print('Evaluating simulated biomarkers')
    simulation_json_file = 'UKB_baseline_simulation_baseline.json'
    alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[0] + '_' + simulation_name + '/'
    # simulation_json_file = 'UKB_1097970_acute_bz1_baseline_simulation_em.json'
    # alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[0] + '_literature_parameters_' + simulation_name + '/'
    # alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[
    #     0] + '_' + simulation_name + '_mec_baseline/'
    pp = PostProcessing(alya=alya, simulation_json_file=simulation_json_file,
                        alya_output_dir=alya_output_dir, protocol='raw', verbose=verbose)
    beat = 1
    pp.read_ecg_pv()
    pp.evaluate_pv_biomarkers(beat=beat)
    pp.evaluate_ecg_biomarkers(beat=beat, show_landmarks=False)
    # pp.evaluate_deformation_biomarkers(beat=beat)
    # pp.evaluate_fibre_work_biomarkers(beat=beat)
    # pp.evaluate_strain_biomarkers(beat=beat)
    pp.visualise_calibration_comparisons_global(beat=1, save_filename=alya_output_dir + '/calibration_result_pvecg.png')
    # pp.visualise_qoi_comparisons(qoi_names = ['qrs_dur_mean', 'qt_dur_mean', 't_pe_mean', 'EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection',
    #                  'dvdt_filling', 'dpdt_max'], save_figure=alya_output_dir+'/qoi_evaluation.png')
    # Rank QoIs according to importance for matching - based on ??

    # # Use OAT SA results and evaluation of biomarkers to assign new ranges for calibration SA
    # ranked_qoi_names = ['LVEF', 'PmaxL', 'SVL', 'EDVL', 'ESVL', 'dpdt_max',
    #                     'dvdt_ejection', 'dvdt_filling']
    # baseline_json_file = 'UKB_baseline_simulation_em_literature_parameters.json'
    # simulation_json_file = baseline_json_file
    # simulation_dict = json.load(open(simulation_json_file, 'r'))
    # pp.calculate_calibration_sa_parameter_ranges(ranked_qoi_names=ranked_qoi_names,
    #                                              oat_sa_slopes='SA_summary_OAT_slopes.csv',
    #                                              oat_sa_p_values='SA_summary_OAT_p_values.csv',
    #                                              oat_sa_r_values='SA_summary_OAT_r_values.csv',
    #                                              oat_sa_ranges='SA_summary_OAT_ranges.csv',
    #                                              oat_sa_intercepts='SA_summary_OAT_intercepts.csv',
    #                                              simulation_dict=simulation_dict, strategy='one_qoi_at_a_time', qoi_input='LVEF')
    # # pp.visualise_calibration_sa_parameter_ranges(oat_sa_slopes='SA_summary_OAT_slopes.csv',
    # #                                              oat_sa_p_values='SA_summary_OAT_p_values.csv',
    # #                                              oat_sa_r_values='SA_summary_OAT_r_values.csv',
    # #                                              oat_sa_ranges='SA_summary_OAT_ranges.csv',
    # #                                              oat_sa_intercepts='SA_summary_OAT_intercepts.csv',
    # #                                              simulation_root_dir=simulation_root_dir)
    quit()

########################################################################################################################
# Step 7: Use OAT SA results and evaluation of biomarkers to assign new ranges for calibration SA
# calibration_folder_name = 'calibration_simulations_third_iteration'
calibration_iteration = '1'
calibration_folder_name = 'UKB_'+mesh_number+'_calibration_simulations_' + calibration_iteration + '_iteration'
baseline_json_file = 'UKB_'+mesh_number+'_baseline_simulation_baseline_' + calibration_iteration + '_iteration.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
perturbed_parameters = json.load(open('UKB_' + mesh_number + '_calibration_sa_ranges_' + calibration_iteration + '_iteration.json', 'r'))
perturbed_parameters_name = np.array(list(perturbed_parameters.keys()))
print(perturbed_parameters_name)
if system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/UKB_baseline_simulation_em_literature_parameters_UKB_05_fine/'
    simulation_dir = calibration_folder_name + '/'
elif system == 'archer2':
    baseline_dir = simulation_root_dir + 'UKB_' + mesh_number + '/'
    simulation_dir = simulation_root_dir + calibration_folder_name + '/'
upper_bounds = []
lower_bounds = []
for param in perturbed_parameters_name:
    lower_bounds.append(perturbed_parameters[param][0])
    upper_bounds.append(perturbed_parameters[param][1])
calibration = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 3 , parameter_names=perturbed_parameters_name,
                   baseline_json_file=baseline_json_file, simulation_dir=simulation_dir, alya_format=alya,
                   baseline_dir=baseline_dir, verbose=verbose)
# calibration = SAUQ(name='sa', sampling_method='uniform', n=8 , parameter_names=perturbed_parameters_name,
#                    baseline_json_file=baseline_json_file, simulation_dir=simulation_dir, alya_format=alya,
#                    baseline_dir=baseline_dir, verbose=verbose)
if setup_calibration_alya_simulations:
    calibration.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
if run_alya_calibration_simulations:
    # calibration.run_jobs(simulation_dir, start_id=0) # Maximum job submission in archer2 is 128 for QoS:taskfarm.
    calibration.run_jobs(simulation_dir, start_id=0)

if evaluate_calibration_sa:
    calibration.sort_simulations(tag='raw')  # Collate list of finished simulations by checking the existence of particular files.
    # calibration.visualise_finished_parameter_sets(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
    evaluate_pv = True
    evaluate_ecg = False
    evaluate_deformation = False
    evaluate_fibrework = False
    evaluate_strain = False
    beat = 1
    labels = []
    if len(perturbed_parameters_name) == 1:
        param = perturbed_parameters_name[0]
        for param_i in range(8):
            filename = simulation_dir + 'sa_' + str(param_i) + '.json'
            simulation_dict = json.load(open(filename, 'r'))
            if 'sf_' in param:
                labels.append(param + '=' + str(simulation_dict[param][0][0]))
            elif '_lv' in param:
                labels.append(param + '=' + str(simulation_dict[param.split('_lv')[0]][0]))
            elif '_rv' in param:
                labels.append(param + '=' + str(simulation_dict[param.split('_rv')[0]][1]))
            elif '_myocardium' in param:
                labels.append(param + '=' + str(simulation_dict[param.split('_myocardium')[0]][0]))
            elif '_valveplug' in param:
                labels.append(param + '=' + str(simulation_dict[param.split('_valveplug')[0]][1]))
            else:
                labels.append(param + '=' + str(simulation_dict[param]))
        print(labels)
    if evaluate_pv:
        pv_post = calibration.evaluate_qois(qoi_group_name='pv', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                             analysis_type='sa')
        calibration.visualise_sa(beat=1, pv_post=pv_post, save_filename=simulation_dir + '/pv_post', labels=labels, show=True, highlight_max_lvef=True)
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
        calibration.visualise_sa(beat=1, pv_post=ecg_post, save_filename=simulation_dir + '/ecg_post')
        qoi_names = ['qrs_dur_mean', 't_dur_mean', 'qt_dur_mean', 't_pe_mean', 'jt_dur_mean']
        slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
            filename=simulation_dir + 'ecg_qois.csv', qois=qoi_names, show_healthy_ranges=False,
            save_filename=simulation_dir + '/ecg_scatter.png')
    elif evaluate_deformation:
        deformation_post = calibration.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat, qoi_save_dir=simulation_dir,
                                             analysis_type='sa')
        calibration.visualise_sa(beat=1, pv_post=deformation_post, save_filename=simulation_dir + '/deformation_post')
        qoi_names = ['es_ed_avpd', 'es_ed_apical_displacement', 'diff_lv_wall_thickness', 'percentage_volume_change']
        slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
            filename=simulation_dir + 'deformation_qois.csv', qois=qoi_names, show_healthy_ranges=False,
            save_filename=simulation_dir + '/deformation_scatter.png')

    elif evaluate_strain:
        strain_post = calibration.evaluate_qois(qoi_group_name='ecg', alya=alya, beat=beat,
                                                     qoi_save_dir=simulation_dir,
                                                     analysis_type='sa')
        calibration.visualise_sa(beat=1, pv_post=strain_post, save_filename=simulation_dir + '/strain_post')
        qoi_names = ['max_mid_Ecc', 'min_mid_Ecc', 'max_mid_Err', 'min_mid_Err', 'max_four_chamber_Ell',
                     'min_four_chamber_Ell']
        slopes, intercepts, p_values, r_values, ranges, params, qois = calibration.analyse(
            filename=simulation_dir + 'strain_qois.csv', qois=qoi_names, show_healthy_ranges=False,
            save_filename=simulation_dir + '/strain_scatter.png')

quit()
########################################################################################################################
# Step 6: Validation experiments - volume perturbation
validation_folder_name = 'volume_perturbation_validation_experiments'
baseline_json_file = 'UKB_1097970_acute_bz1_baseline_simulation_em.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
perturbed_parameters_name = np.array(['end_diastole_p_lv'])
baseline_parameter_values = np.array([simulation_dict['end_diastole_p'][0]])
upper_bounds = baseline_parameter_values * 4.0
lower_bounds = baseline_parameter_values * 0.8
baseline_dir = ''
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/UKB_baseline_simulation_em_UKB_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + validation_folder_name + '/'
elif system == 'cosma':
    baseline_dir = simulation_root_dir + 'UKB_baseline_simulation_em_UKB_05_fine_mec_baseline/'
    simulation_dir = simulation_root_dir + validation_folder_name + '/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/UKB_baseline_simulation_em_UKB_05_fine/'
    simulation_dir = validation_folder_name + '/'
elif system == 'archer2':
    baseline_dir = simulation_root_dir + 'UKB_baseline_simulation_em_UKB_05_fine_mec_baseline/'
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