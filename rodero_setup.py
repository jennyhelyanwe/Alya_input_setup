from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from alyaformat import AlyaFormat


mesh_number = '01'
vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number + '/'
geometric_data_dir = '/data/Personalisation_projects/meta_data/geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
boundary_data_dir = geometric_data_dir + 'boundary_data/'
field_data_dir = geometric_data_dir + 'field_data/'
personalisation_data_dir = '/data/Personalisation_projects/meta_data/results/personalisation_data/rodero_'+mesh_number+'/'
vtk_name = mesh_number + '_bivent_only'
output_name = 'rodero_' + mesh_number + '_fine'
simulation_json_file = 'rodero_baseline_simulation.json'
verbose = False
# MeshPreprocessing(vtk_name=vtk_name, name=output_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
#                                boundary_data_dir=boundary_data_dir, field_data_dir=field_data_dir, verbose=verbose)

# FieldGeneration(name=output_name, geometric_data_dir=geometric_data_dir, boundary_data_dir=boundary_data_dir,
#                 field_data_dir=field_data_dir, personalisation_data_dir=personalisation_data_dir, verbose=verbose)

AlyaFormat(name=output_name, geometric_data_dir=geometric_data_dir, boundary_data_dir=boundary_data_dir,
           field_data_dir=field_data_dir, personalisation_dir=personalisation_data_dir, simulation_json_file=simulation_json_file, verbose=verbose)

# alya = AlyaFormat(name=output_name, geometric_data_dir=geometric_data_dir, boundary_data_dir=boundary_data_dir, field_data_dir=field_data_dir,
#                   simulation_json_file='rodero_baseline_simulation.json')
