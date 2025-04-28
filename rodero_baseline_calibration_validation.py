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

max_cores_used = 3
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
evaluate_calibrated_baseline = False
setup_calibration_alya_simulations = False
run_alya_calibration_simulations = False
evaluate_calibration_sa = False
visualise_calibration_result = True
setup_validation_alya_simulations = False
run_alya_validation_simulations = False
run_alya_validation_postprocessing = False
evaluate_validation_biomarkers = False

iteration = '1st'

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
    mesh = MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
                             max_cores_used=max_cores_used, verbose=verbose)
    # mesh.read_geometry_from_vtk_rodero(save=False)
    # mesh.generate_boundary_data_rodero(save=True)
    # mesh.check_fields_for_qrs_inference()
    # mesh.check_fields_for_twave_inference()
    mesh.generate_lower_resolution_mesh()

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
                personalisation_data_dir=personalisation_data_dir, max_cores_used=max_cores_used, verbose=verbose)

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
                  simulation_dir = simulation_root_dir, job_version=system, max_cores_used=max_cores_used, verbose=verbose)

# Sanity check:
# if not system == 'jureca' and not system == 'cosma' and not system == 'archer2':
#     alya.visual_sanity_check(simulation_json_file=simulation_json_file)
if setup_em_alya_files:
    simulation_json_file = 'rodero_baseline_simulation_baseline_' + iteration + '_iteration.json'
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
if evaluate_calibrated_baseline:
    print('Evaluating simulated biomarkers')
    simulation_json_file = 'rodero_baseline_simulation_baseline_' + iteration + '_iteration.json'
    alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[0] + '_' + simulation_name + '/'
    # simulation_json_file = 'rodero_baseline_simulation_em.json'
    # alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[0] + '_literature_parameters_' + simulation_name + '/'
    # alya_output_dir = simulation_root_dir + simulation_json_file.split('/')[-1].split('.')[
    #     0] + '_' + simulation_name + '_mec_baseline/'
    pp = PostProcessing(alya=alya, simulation_json_file=simulation_json_file,
                        alya_output_dir=alya_output_dir, protocol='raw', verbose=verbose)
    beat = 1
    pp.read_ecg_pv()
    # pp.evaluate_pv_biomarkers(beat=beat)
    # pp.evaluate_ecg_biomarkers(beat=beat, show_landmarks=False)
    pp.evaluate_deformation_biomarkers(beat=beat)
    # pp.evaluate_fibre_work_biomarkers(beat=beat)
    # pp.evaluate_strain_biomarkers(beat=beat)
    pp.visualise_calibration_comparisons_global(beat=1, save_filename=alya_output_dir + '/calibration_result_global.png')
    # pp.visualise_qoi_comparisons(qoi_names = ['qrs_dur_mean', 'qt_dur_mean', 't_pe_mean', 'EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection',
    #                  'dvdt_filling', 'dpdt_max'], save_figure=alya_output_dir+'/qoi_evaluation.png')
    # Rank QoIs according to importance for matching - based on ??

    # # Use OAT SA results and evaluation of biomarkers to assign new ranges for calibration SA
    # ranked_qoi_names = ['LVEF', 'PmaxL', 'SVL', 'EDVL', 'ESVL', 'dpdt_max',
    #                     'dvdt_ejection', 'dvdt_filling']
    # baseline_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
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

# Fourth iteration: instead of adjusting ICaL, which could bring on arrhythmic effects, it's better to alter the cross bridge cycling
# dynamics and compressibility.
# {"sfkws_myocardium": [4, 7], "Kct_myocardium": [100000.0,  1000000.0]}

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
# Max LVEF: 42%,

# Seventh iteration: add compliance
# Set: sfkws: 5, set cal50: 0.5, set Kct: 450,000
# Search: R_LV: [700, 999], ejection threshold: [7, 10], C_LV: [0.000065, 0.00015]

# I need to lower R_LV even more... this may not be very justified, but the peak pressure is almost double what it should be.
# From the scatter plot, it isn't clear which parameter has the strongest effect on LVEF, there seems to be a lot of interaction
# between these parameters... Perhaps best to fix the ejection threshold and compliance and only tune R at this point...

# Eighth iteration: only resistance
# Set: 5, set cal50: 0.5, set Kct: 450,000, ejection threshold: 7 kPa, C: 0.00010
# Search: R_LV: [400, 700]

# Ninth iteration: Even with really low RV_LV: 450, the peak pressure is still too high (20 kPa). So, I'll now leave RV_LV at 500, and
# do a calibration of sf_kws down from 5 to 1, and see if we can get a higher LVEF...

# Tenth iteration: Sf_kws doesn't do much to reduce the peak pressure...
# Going to just drastically reduce the R... this could be justified because the Windkessel model is the least realistic model
# out of all the options that we can change to alter LVEF.

# 11th iteration: Lowering the R is working. Based on the linear regression for the LVEF scatter plots, R= 200 is going to
# give LVEF: 52%, and ESP of 14.2 kPa, which would fall perfectly in range. I'll run this now as a single simulation.
# Then, next step, I will need to tune the PFR and PERs. PER is strongly dependent on R, so we may need a separate way of
# getting good LVEF while having slow PERs... This will be tricky.

# 12th iteration: So, the current problem is that PER is too high, and PFR is too low.
# Decreasing sf_kws will decrease PER, but it comes at the expense of decreasing LVEF and also decreasing peak pressure.
# turns out, gain contraction has no effect on the PER or dPdt. Investigating the effect of conduction velocity on this...

# 13th iteration: tune gain relaxation to get higher PFR. Bearing in mind that this may need to be retuned when we get
# slower ejection speed...
# {"gain_error_relaxation_lv":  [0.005, 0.05]}
# gain_error_relaxation_lv=0.035 is going to give a dVdt filling of roughly 300, which will sit right in the middle of the
# healthy range.

# 14th iteration: Here comes the delicate balance. What we want is:
# LVEF, ESP, dPdt, dvdt ejection, dvdt filling
# currently, dVdt ejection, dPdt max, and percentage volume change are all too high.
# The parameters that affect dPdt are: skws, Tref, and cal50, GCaL. Of these, sf_kws has the strongest effect, and is most promising.
# The parameters that affect dVdt ejection: Tref, cal50, sf_kws, Kct (higher it is, slower the ejection), GCaL. The dVdt ejection is
# strongly correlated with ESV.
# The set of parameter that increases LVEF or decreases ESV: Tref, cal50, sf_kws, Kct, R_LV, GCaL, Jup (the only parameter that can decrease ESV without changing the dVdt is Jup).
# The parameters that strongly affect peak pressure are: Tref, sf_kws, pericardial stiffness, Kct, af, R_LV, C_LV (affects peak pressure without changing anything else),
# ejection threshold, GCaL, Jup

# So, to decrease dPdt and dVdt ejection, I will decrease sf_Kws to below 1.0. This is going to cause a decrease in peak pressure and LVEF.
# To get the LVEF back up, I will decrease Jup to 0.5 and increase GCaL to 2, and potentially change C_LV to get peak pressure back up.

# So, 14th iteration:
# Set: sf_jup: 0.5, sf_gcal: 2, gain_error_relaxation_lv: 0.035
# Search: sf_kws: [0.2, 1.0] to get good dPdt and dVdt ejection.

# So that didn't work. sf)gcal at 2 and sf_jup at 0.5 is too much, the diastolic phase is stunted and the mechanics
# fails before contraction happens.

# Perhaps better to do a SA for this and search within the sf_gcal and sf_jup space
# 15th iteration:
# search: sf_kws: [0.2, 1.0], sf_jup: [0.5, 1.0]


########################################################################################################################
# Step 7: Use OAT SA results and evaluation of biomarkers to assign new ranges for calibration SA
# calibration_folder_name = 'calibration_simulations_third_iteration'
calibration_folder_name = 'calibration_simulations_' + iteration + '_iteration'
baseline_json_file = 'rodero_baseline_simulation_em_literature_parameters.json'
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
# perturbed_parameters = json.load(open('calibration_sa_ranges_third_iteration.json', 'r'))
perturbed_parameters = json.load(open('calibration_sa_ranges_' + iteration + '_iteration.json', 'r'))
perturbed_parameters_name = np.array(list(perturbed_parameters.keys()))
if system == 'heart':
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
calibration = SAUQ(name='sa', sampling_method='lhs', n=70 , parameter_names=perturbed_parameters_name,
                   baseline_json_file=baseline_json_file, simulation_dir=simulation_dir, alya_format=alya,
                   baseline_dir=baseline_dir, max_cores_used=max_cores_used, verbose=verbose)
# calibration = SAUQ(name='sa', sampling_method='uniform', n=8 , parameter_names=perturbed_parameters_name,
#                    baseline_json_file=baseline_json_file, simulation_dir=simulation_dir, alya_format=alya,
#                    baseline_dir=baseline_dir, verbose=verbose)
calibrated_simulation_dir = simulation_root_dir + 'calibrated_baseline_' + iteration + '_iteration_rodero_'+mesh_number + '/'
calibrated_json_filename = simulation_dir + 'calibrated_baseline_' + iteration + '_iteration_rodero_'+mesh_number + '.json'
if setup_calibration_alya_simulations:
    calibration.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
if run_alya_calibration_simulations:
    # calibration.run_jobs(simulation_dir, start_id=0) # Maximum job submission in archer2 is 128 for QoS:taskfarm.
    calibration.run_jobs(simulation_dir, start_id=0)
# calibration.run_jobs(simulation_dir, start_id=0, end_id=128) # Maximum job submission in archer2 is 128 for QoS:taskfarm.
if evaluate_calibration_sa:
    print('Evaluate calibration sampling')
    if system == 'archive':
        calibration.sort_simulations_archive(tag='raw')
    else:
        calibration.sort_simulations(
            tag='raw')  # Collate list of finished simulations by checking the existence of particular files.
    # calibration.visualise_finished_parameter_sets(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
    evaluate_pv = False
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


    calibration_qoi_names = ['LVEF','PmaxL','ESVL','EDVL', 'dvdt_ejection', 'dvdt_filling']
    calibration_qoi_weights = [5, 5, 2, 0.8, 0.8, 2]
    best_simulation_id = calibration.select_best_calibration_result(pv_qois_filename=simulation_dir+'pv_qois.csv',
                                               qoi_names=calibration_qoi_names,
                                               qoi_weights=calibration_qoi_weights,
                                               image_save_dir=simulation_dir,
                                               calibrated_simulation_dir=calibrated_simulation_dir,
                                               calibrated_json_filename=calibrated_json_filename, show=True)

if visualise_calibration_result:
    # Visualise calibrated model against healthy ranges
    alya_output_dir = calibrated_simulation_dir
    simulation_json_file = calibrated_json_filename
    pp = PostProcessing(alya=alya, simulation_json_file=simulation_json_file,
                        alya_output_dir=alya_output_dir, protocol='raw', max_cores_used=max_cores_used,
                        verbose=verbose)
    beat = 1
    # pp.evaluate_strain_biomarkers(beat=beat)
    # pp.evaluate_fibre_work_biomarkers(beat=beat)
    # pp.visualise_calibration_comparisons_strain()
    pp.read_ecg_pv()
    # pp.evaluate_deformation_biomarkers(beat=beat)
    pp.visualise_calibration_comparisons_global(beat=beat,
                                                save_filename=simulation_dir + '/calibration_result.png',
                                                show=True)

########################################################################################################################
# Step 6: Validation experiments - volume perturbation
print('Volume perturbation experiment')
validation_folder_name = 'volume_perturbation_validation_experiments'
baseline_json_file = calibrated_json_filename
simulation_json_file = baseline_json_file
simulation_dict = json.load(open(simulation_json_file, 'r'))
perturbed_parameters_name = np.array(['end_diastole_p_lv'])
baseline_parameter_values = np.array([simulation_dict['end_diastole_p'][0]])
upper_bounds = baseline_parameter_values * 5.0
lower_bounds = baseline_parameter_values * 0.5
baseline_dir = calibrated_simulation_dir
simulation_dir = simulation_root_dir + validation_folder_name + '/'
experiment = SAUQ(name='sa', sampling_method='uniform', n=16, parameter_names=perturbed_parameters_name,
                  baseline_json_file=baseline_json_file, simulation_dir=simulation_dir, alya_format=alya,
                  baseline_dir=baseline_dir, max_cores_used=max_cores_used, verbose=verbose)
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