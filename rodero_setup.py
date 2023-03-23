from vtk_to_csv_ensight import *
from alyaformat import AlyaFormat


mesh_number = '05'
vtk_dir = '/users/jenang/RoderoNiedererMeshHealthy/' + mesh_number
csv_dir = '/data/Personalisation_projects/meta_data/geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
vtk_name = mesh_number + '_bivent_only'
output_name = 'rodero_' + mesh_number + '_fine'

# vtk_reader = VTKtoCSVEnsight(vtk_name=vtk_name, output_name=output_name, input_dir=vtk_dir, output_dir=csv_dir)

alya = AlyaFormat(name=output_name, geometry_and_fields_input_dir=csv_dir,
                  simulation_json_file='rodero_baseline_simulation.json')
