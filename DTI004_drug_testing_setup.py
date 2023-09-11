from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from alyaformat import AlyaFormat
from populationdrugtest import PopulationDrugTest
from postprocessing import PostProcessing
import numpy as np
import os
########################################################################################################################
# Global Settings
mesh_number = '05'
simulation_name = 'DTI004'
workdir = os.getcwd()
if 'icei' in workdir:
    system = 'jureca'
else:
    system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/'+simulation_name+'/'+simulation_name+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = True

#######################################################################################################################
# # Step 1: Save input mesh into CSV format, as prescribed in myformat.py
# if system == 'heart':
#     vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number + '/'
# elif system == 'jureca':
#     vtk_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/geometric_data/vtk/'
# vtk_name = mesh_number + '_bivent_only'
# simulation_json_file = 'rodero_baseline_simulation_ep.json'
# MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
#                   simulation_json_file=simulation_json_file, four_chamber=False, verbose=verbose)
#
# Step 4: Generate fields for Alya simulation
electrode_data_filename = meta_data_dir + 'geometric_data/'+simulation_name+'/'+simulation_name+'_electrode_xyz.csv'
personalisation_data_dir = meta_data_dir + 'personalisation_data/'+simulation_name+'/twave_sf_IKs_GKs5_GKr0.6_tjca60/drug_testing_population/'
fields = FieldGeneration(name=simulation_name+'_fine', geometric_data_dir=geometric_data_dir,
                         personalisation_data_dir=personalisation_data_dir, verbose=verbose)
# fields.generate_celltype(read_celltype_filename=personalisation_data_dir+simulation_name+'_fine_nodefield_personalisation-biomarker_0.csv')
# fields.generate_electrode_locations(electrode_data_filename=electrode_data_filename, save=True)

########################################################################################################################
# Step 5: Set up Alya input files according to simulation protocol saved in .json file.
simulation_dir = ''
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_dir = './'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)

########################################################################################################################
# Step 6: Set up population based drug tests
baseline_json_file = simulation_name+'_baseline_simulation_ep.json'
baseline_dir = ''
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'+simulation_name+'_baseline_simulation_ep_'+simulation_name+'_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/population_drug_test/'+simulation_name+'/'

paths = os.listdir(personalisation_data_dir)
population_size = 0
for path in paths:
    if 'personalisation-biomarker' in path:
        population_size = population_size + 1
pdt = PopulationDrugTest(name=simulation_name, personalised_population_dir=personalisation_data_dir, population_size=population_size,
                          alya_format=alya, verbose=verbose)
drug_doses = [{'sf_gkr':0.6}, {'sf_gkr':0.5}, {'sf_gkr':0.4}, {'sf_gkr':0.34}, {'sf_gkr':0.3}, {'sf_gkr':0.27}, {'sf_gkr':0.25}]
pdt.setup_drug_test(drug_name='dofetilide', drug_doses=drug_doses, simulation_dir=simulation_dir,population_sf_names=['sf_IKs'],
                    baseline_json_file=baseline_json_file, baseline_dir=baseline_dir)
#
pdt.run_jobs(simulation_dir=simulation_dir, start_id=0)
# pdt.visualise_drug_effect_ecgs(beat=1, drug_name='dofetilide', drug_doses=drug_doses, simulation_dir=simulation_dir)
# pdt.evaluate_drug_effect_ecgs(drug_name='dofetilide', drug_doses=drug_doses, simulation_dir=simulation_dir)