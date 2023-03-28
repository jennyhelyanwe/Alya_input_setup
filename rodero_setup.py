from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from alyaformat import AlyaFormat

########################################################################################################################
# Global Settings
mesh_number = '05'
simulation_name = 'rodero_' + mesh_number + '_fine'
geometric_data_dir = '/data/Personalisation_projects/meta_data/geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
verbose = False

#######################################################################################################################
# Step 1: Save input mesh into CSV format, as prescribed in myformat.py
vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number + '/'
vtk_name = mesh_number + '_bivent_only'
MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
                  verbose=verbose)

########################################################################################################################
# Step 2: Run QRS and T inference and write personalised results to personalisation_data_dir


########################################################################################################################
# Step 3: Generate fields for Alya simulation
personalisation_data_dir = '/data/Personalisation_projects/meta_data/results/personalisation_data/rodero_'+mesh_number+'/'
electrode_data_filename = '/data/Personalisation_projects/meta_data/geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_electrode_xyz.csv'
FieldGeneration(name=simulation_name, geometric_data_dir=geometric_data_dir, electrode_data_filename=electrode_data_filename,
                personalisation_data_dir=personalisation_data_dir, verbose=verbose)

########################################################################################################################
# Step 4: Write Alya input files according to simulation protocol saved in .json file.
simulation_json_file = 'rodero_baseline_simulation.json'
AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir, personalisation_dir=personalisation_data_dir,
           simulation_json_file=simulation_json_file, verbose=verbose)


