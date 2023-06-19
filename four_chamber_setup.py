import numpy as np
from meshpreprocessing import MeshPreprocessing

# Global Settings
mesh_number = '05'
simulation_name = 'rodero_' + mesh_number + '_fine'
system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine_four_chamber/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False

#######################################################################################################################
# Step 1: Save input mesh into CSV format, as prescribed in myformat.py
if system == 'heart':
    vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number + '_four_chamber/'
elif system == 'jureca':
    vtk_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/geometric_data/vtk/'
# vtk_name = mesh_number + '_bivent_only'
vtk_name = mesh_number
simulation_json_file = 'rodero_baseline_simulation_ep.json'
MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
                  simulation_json_file=simulation_json_file, four_chamber=True, verbose=verbose)



