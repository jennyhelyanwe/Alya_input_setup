from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from alyaformat import AlyaFormat
from populationdrugtest import PopulationDrugTest
from postprocessing import PostProcessing
import numpy as np
from calibratecv import CalibrateCV
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
    personalisation_data_dir = meta_data_dir + '/personalisation_data/' + simulation_name + '/twave_sf_IKs_GKs5_GKr0.6_tjca60/translation_to_monodomain/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/'+simulation_name+'/'+simulation_name+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False


#######################################################################################################################
# # Calibrate CV
# if system == 'heart':
#     calibration_dir = '/users/jenang/Alya_setup_SA/calibration_dir/'
#     alya_executable_path = '/data/Personalisation_projects/alya-compbiomed2_ionic_sf_fields/Executables/unix/Alya.x'
#     personalisation_data_dir = meta_data_dir + '/results/personalisation_data/DTI004/qrs/'
# calibrate = CalibrateCV(name=simulation_name + '_fine', geometric_data_dir=geometric_data_dir, calibration_dir=calibration_dir,
#                         alya_executable_path=alya_executable_path, personalisation_data_dir=personalisation_data_dir, verbose=verbose)
# quit()
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


fields = FieldGeneration(name=simulation_name+'_fine', geometric_data_dir=geometric_data_dir,
                         personalisation_data_dir=personalisation_data_dir, verbose=verbose)

# fields.generate_celltype(read_celltype_filename=personalisation_data_dir+simulation_name+'_fine_nodefield_personalisation_biomarker_best.csv')
# fields.generate_electrode_locations(electrode_data_filename=electrode_data_filename, save=True)
# fields.generate_stimuli(lat_filename=personalisation_data_dir+'/'+simulation_name+'_fine_nodefield_personalisation-lat.csv')
# fields.generate_ionic_scaling_factors(
#     read_biomarker_filename=personalisation_data_dir + '/' + simulation_name + '_fine_nodefield_personalisation_biomarker_best.csv', save=True)

########################################################################################################################
# Step 5: Simulate best match model.
simulation_dir = ''
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/twave_best_match/'
alya = AlyaFormat(name=simulation_name + '_fine', geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)

baseline_json_file = simulation_name+'_baseline_simulation_ep.json'
baseline_dir = ''
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/twave_best_match/'+simulation_name+'_baseline_simulation_ep_'+simulation_name+'_fine/'

#fields.generate_ionic_scaling_factors(
#    read_biomarker_filename=personalisation_data_dir + '/' + simulation_name + '_fine_nodefield_personalisation_biomarker_best.csv', save=True)
#alya.do(simulation_json_file=baseline_json_file, SA_flag=False, drug_flag=False,
#        best_match_biomarker_file=personalisation_data_dir+simulation_name+'_fine_nodefield_personalisation_biomarker_best.csv',
#        baseline_dir=baseline_dir)
# alya.do(simulation_json_file=baseline_json_file, SA_flag=False, drug_flag=False, baseline_dir=baseline_dir)
#
# baseline_json_file = simulation_name+'_baseline_simulation_ep_wideAPD.json'
# baseline_dir = ''
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/twave_best_match/'+simulation_name+'_baseline_simulation_ep_'+simulation_name+'_fine_with_epi/'
#
# fields.generate_ionic_scaling_factors(
#     read_biomarker_filename=personalisation_data_dir + '/' + simulation_name + '_fine_nodefield_personalisation_biomarker_best.csv', save=True)
# alya.do(simulation_json_file=baseline_json_file, SA_flag=False, drug_flag=False,
#         best_match_biomarker_file=personalisation_data_dir+simulation_name+'_fine_nodefield_personalisation_biomarker_best.csv',
#         baseline_dir=baseline_dir)


#quit()
########################################################################################################################
# Step 6: Set up population based drug tests
personalisation_data_dir = meta_data_dir + 'personalisation_data/'+simulation_name+'/twave_sf_IKs_GKs5_GKr0.6_tjca60/drug_testing_population/'
simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/population_drug_test/'+simulation_name+'/'
alya = AlyaFormat(name=simulation_name + '_fine', geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)
paths = os.listdir(personalisation_data_dir)
population_size = 0
for path in paths:
    if 'personalisation-biomarker' in path:
        population_size = population_size + 1
pdt = PopulationDrugTest(name=simulation_name + '_fine', personalised_population_dir=personalisation_data_dir, population_size=population_size,
                          alya_format=alya, verbose=verbose)
drug_doses = [{'sf_gkr':0.6}, {'sf_gkr':0.5}, {'sf_gkr':0.4}, {'sf_gkr':0.34}, {'sf_gkr':0.3}, {'sf_gkr':0.27}, {'sf_gkr':0.25}]
setup_and_run = False
postprocess = False
visualise = False
evaluate_ecg = True
download = False
if setup_and_run:
    pdt.setup_drug_test(drug_name='dofetilide', drug_doses=drug_doses, simulation_dir=simulation_dir,population_sf_names=['sf_IKs'],
                        baseline_json_file=baseline_json_file, baseline_dir=baseline_dir)
    pdt.run_jobs(simulation_dir=simulation_dir, start_id=0)
elif postprocess:
    pdt.run_jobs_postprocess(simulation_dir=simulation_dir)
elif visualise:
    pdt.visualise_drug_effect_ecgs(beat=1, drug_name='dofetilide', drug_doses=drug_doses, simulation_dir=simulation_dir)
elif download:
    dir_for_download = '/p/project/icei-prace-2022-0003/wang1/Personalisation_projects/meta_data/results/personalisation_data/DTI004/twave_sf_IKs_GKs5_GKr0.6_tjca60/only_endo/monodomain_drug_test_results/'
    # simulation_dir + '/monodomain_drug_test_results/'
    pdt.extract_vms_for_download(simulation_dir=simulation_dir, dir_for_download=dir_for_download, drug_doses=drug_doses, drug_name='dofetilide')
elif evaluate_ecg:
    pdt.evaluate_drug_effect_ecgs(drug_name='dofetilide', drug_doses=drug_doses, simulation_dir=simulation_dir)
