from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration

mesh_numbers = [ '02', '03', '04', '06', '07', '08', '09', '10'] # Process all other rodero meshes except for 05 which is already ready.
max_cores_used = 10
verbose = True
for mesh_number in mesh_numbers:
    simulation_name = 'rodero_' + mesh_number + '_fine'
    # Step 1: Save input mesh into CSV format, as prescribed in myformat.py
    meta_data_dir = '/work/e769/e769/jennywang/Alya_pipeline/meta_data/'
    vtk_dir = meta_data_dir + '/geometric_data/rodero_vtk/' + mesh_number + '/'
    vtk_name = mesh_number + '_bivent_only'
    simulation_json_file = 'rodero_baseline_simulation_ep.json'
    geometric_data_dir = meta_data_dir + 'geometric_data/rodero_' + mesh_number + '/rodero_' + mesh_number + '_fine/'
    mesh = MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
                             max_cores_used=max_cores_used, verbose=verbose)
    mesh.read_geometry_from_vtk_rodero(save=False)
    mesh.generate_boundary_data_rodero(surface_label_json=vtk_dir+'/surface_labels.json', save=True)
    # mesh.generate_lower_resolution_mesh(simulation_name='rodero_' + mesh_number + '_coarse')

