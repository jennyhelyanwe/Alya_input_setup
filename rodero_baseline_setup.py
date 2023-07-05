from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from calibratecv import CalibrateCV
from alyaformat import AlyaFormat
from sensitivityanalysis import SA, run
from postprocessing import PostProcessing
import numpy as np
import os
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

#######################################################################################################################
# Step 1: Save input mesh into CSV format, as prescribed in myformat.py
if system == 'heart':
    vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number + '/'
elif system == 'jureca':
    vtk_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/geometric_data/vtk/'
vtk_name = mesh_number + '_bivent_only'
simulation_json_file = 'rodero_baseline_simulation_ep.json'
# MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
#                   simulation_json_file=simulation_json_file, four_chamber=False, verbose=verbose)

########################################################################################################################
# Step 2: Run QRS and T inference and write personalised results to personalisation_data_dir

########################################################################################################################
# Step 3: Calibrate conductivities to personalisation conduction velocity
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
# calibration_dir = meta_data_dir + 'calibration_dir/'
# simulation_json_file = 'rodero_baseline_simulation_ep.json'
# calibrate = CalibrateCV(name=simulation_name, geometric_data_dir=geometric_data_dir,
#                         calibration_dir=calibration_dir, simulation_json_file=simulation_json_file,
#                         personalisation_data_dir=personalisation_data_dir, verbose=verbose)

########################################################################################################################
# Step 4: Generate fields for Alya simulation
electrode_data_filename = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_electrode_xyz.csv'
# FieldGeneration(name=simulation_name, geometric_data_dir=geometric_data_dir, electrode_data_filename=electrode_data_filename,
#                 personalisation_data_dir=personalisation_data_dir, endo_mid_divide=0.3, mid_epi_divide=0.7, verbose=verbose)

########################################################################################################################
# Step 5: Write Alya input files according to simulation protocol saved in .json file.
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_dir = './'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)

# Sanity check:
if not system == 'jureca':
    alya.visual_sanity_check()

simulation_json_file = 'rodero_baseline_simulation_em.json'
alya.do(simulation_json_file=simulation_json_file)
simulation_json_file = 'rodero_baseline_simulation_ep.json'
alya.do(simulation_json_file=simulation_json_file)

########################################################################################################################
# Step 5: Run baseline simulation
# run_job('/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine')
# rub_job_postprocess('/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine')

########################################################################################################################
# Step 6: Postprocess
simulation_json_file = 'rodero_baseline_simulation_em.json'
pp = PostProcessing(name=simulation_name, geometric_data_dir=geometric_data_dir, simulation_json_file=simulation_json_file,
                    alya_output_dir=simulation_dir, verbose=verbose)
pp.evaluate_ecg_pv_biomarkers(beat=0)
# pp.save_sa_results(sa_results_dir+'sa_results.csv')