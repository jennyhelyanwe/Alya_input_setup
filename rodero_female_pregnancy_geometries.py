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
evaluate_simulated_biomarkers = False
setup_validation_alya_simulations = False
run_alya_validation_simulations = False
run_alya_validation_postprocessing = False
evaluate_validation_biomarkers = True

#######################################################################################################################
# Step 1: Save input mesh into CSV format, as prescribed in myformat.py
vtk_dir = ''
if system == 'heart':
    vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number + '/'
elif system == 'jureca':
    vtk_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/geometric_data/vtk/'
vtk_name = mesh_number + '_bivent_only'
simulation_json_file = 'rodero_baseline_simulation_ep.json'
if mesh_preprocess:
    mesh = MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir, verbose=verbose)
    mesh.read_geometry_from_vtk_rodero(save=False)
    mesh.generate_boundary_data_rodero(save=True)
    mesh.check_fields_for_qrs_inference()
    mesh.check_fields_for_twave_inference()

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

