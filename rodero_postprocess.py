from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from calibratecv import CalibrateCV
from alyaformat import AlyaFormat
from postprocessing import PostProcessing

########################################################################################################################
# Global Settings
mesh_number = '05'
simulation_name = 'rodero_' + mesh_number + '_fine'
system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False

########################################################################################################################
# Step 1: Run Alya simulation and post-processing


########################################################################################################################
# Step 2: Read in post-processed Alya files
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_ep_rodero_'+mesh_number+'_fine/'
    sa_results_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/sensitivity_analyses/'
elif system == 'heart':
    simulation_dir = '/users/jenang/Alya_setup_SA/alya_csv/rodero_' + mesh_number + '/'
    sa_results_dir = '/users/jenang/Alya_setup_SA/sensitivity_analyses/'
simulation_json_file = 'rodero_baseline_simulation_ep.json'
pp = PostProcessing(name=simulation_name, simulation_json_file=simulation_json_file,  geometric_data_dir=geometric_data_dir,
                    alyacsv_dir=simulation_dir+'results_csv/', results_dir=sa_results_dir, verbose=verbose)


