from vtk_to_csv import *
from alyaformat import AlyaFormat

vtk_dir = '/users/jenang/RoderoNidererMeshHealthy/01'
csv_dir = '/users/jenang/Alya_setup_SA/csv/'
# csv_dir = '/data/Personalisation_projects/meta_data/geometric_data/rodero_01/'
vtk_name = 'bivent_base_only'
output_name = 'rodero_01'

vtk_reader = VTKtoCSV(vtk_name=vtk_name, output_name=output_name, input_dir=vtk_dir, output_dir=csv_dir)
alya = AlyaFormat(name=output_name, geometry_and_fields_input_dir=csv_dir,
                  simulation_json_file='rodero_01_baseline_simulation.json')
