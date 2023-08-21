from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from calibratecv import CalibrateCV
from alyaformat import AlyaFormat
from postprocessing import PostProcessing
from sensitivityanalysis import SA
import os
import numpy as np

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
# # Step 1: Run Alya simulation and post-processing
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
elif system == 'heart':
    simulation_dir = '/users/jenang/Alya_setup_SA/alya_csv/rodero_' + mesh_number + '/'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)

########################################################################################################################
# Step 2: Read in post-processed Alya files
simulation_json_file = 'rodero_baseline_simulation_em.json'
pp = PostProcessing(alya=alya, simulation_json_file=simulation_json_file,
                    alya_output_dir=simulation_dir, verbose=verbose)
beat = 1
# pp.evaluate_ecg_pv_biomarkers(beat=beat)
# pp.read_csv_fields('EPSXX', 'scalar')
# pp.read_csv_fields('EPSYY', 'scalar')
# pp.read_csv_fields('EPSZZ', 'scalar')
# pp.read_csv_fields('EPSXY', 'scalar')
# pp.read_csv_fields('EPSYZ', 'scalar')
# pp.read_csv_fields('EPSXZ', 'scalar')
pp.read_csv_fields('DISPL', 'vector')
pp.evaluate_basal_displacement()
quit()
pp.evaluate_ecg_pv_biomarkers(beat=0)
# pp.save_sa_results(sa_results_dir+'sa_results.csv')
########################################################################################################################

cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
baseline_parameter_values = np.array([1,1,1,1,1,0.805,1])
baseline_json_file = 'rodero_baseline_simulation_em.json'
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_dir = './'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = 'sensitivity_analyses/'
# sa = SA(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=cell_parameter_names,
#        baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#        simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)

# sa.analyse(sa_results_dir+'/results.txt')



