from vtk_to_csv_ensight import *
from alyaformat import AlyaFormat

vtk_dir = '/users/jenang/RoderoNidererMeshHealthy/01'
csv_dir = '/data/Personalisation_projects/meta_data/geometric_data/rodero_01/rodero_01_fine/'
vtk_name = '01_bivent_only'
output_name = 'rodero_01_fine'

vtk_reader = VTKtoCSVEnsight(vtk_name=vtk_name, output_name=output_name, input_dir=vtk_dir, output_dir=csv_dir)

alya = AlyaFormat(name=output_name, geometry_and_fields_input_dir=csv_dir,
                  simulation_json_file='rodero_01_baseline_simulation.json')
