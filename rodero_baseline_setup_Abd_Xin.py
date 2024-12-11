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
system = 'jureca'
max_cores_used = 32
meta_data_dir = '/data/Exchange/new_rodero_baseline_for_Abdallah_Xin_11Dec2024/Alya_pipeline/meta_data/'
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
setup_em_alya_files = True
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
simulation_root_dir = '/data/Exchange/new_rodero_baseline_for_Abdallah_Xin_11Dec2024/Alya_pipeline/alya_simulations/'

alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_root_dir, job_version=system, max_cores_used=max_cores_used, verbose=verbose)

if setup_em_alya_files:
    output_dir = simulation_root_dir + '/test/'
    simulation_json_file = simulation_root_dir + '/calibrated_baseline_' + iteration + '_iteration_rodero_'+mesh_number+'.json'
    alya.do(simulation_json_file=simulation_json_file, output_dir=output_dir)
if setup_ep_alya_files:
    simulation_json_file = 'rodero_baseline_simulation_ep.json'
    alya.do(simulation_json_file=simulation_json_file)

if setup_em_alya_literature_parameters_files:
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
    simulation_json_file = 'calibrated_baseline_' + iteration + '_iteration_rodero_'+mesh_number+'.json'
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
