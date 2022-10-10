import readfunctions
import writefunctions
import visualisationfunctions
import argparse
import os, sys

false = ['False', 'false', 'F', 'f']
true = ['True', 'true', 'T', 't']
parser = argparse.ArgumentParser()
parser.add_argument("--root_dir", help='Root directory, default is ./', default='./')
parser.add_argument("--name", help='Mesh name', default='heart')
parser.add_argument("--input_dir", help='Directory where the input mesh can be found, default ./', default='./')
parser.add_argument("--save_dir", help='Directory where the saved arrays can be found, default save/', default='save/')
parser.add_argument("--output_dir", help='Directory where the Alya files are written to, default alya/', default='alya/')
parser.add_argument("--input_format", help='Format of input mesh file: [cobiveco, vtk] default cobiveco', default='cobiveco')
parser.add_argument("--refresh", help='Either ecg or pv, or single cell [true, false] default false', default=False)
parser.add_argument("--visualisation_dir", help='Visualisation file directory, default visualisation/', default='visualisation/')
parser.add_argument("--visualisation_format", help='Visualisation type: [ensight, MATLAB, vtk] default ensight', default='ensight')
parser.add_argument("--visualise", help='Toggle whether or not to write visualisation files to check, [true, false]  default false', default=False)
parser.add_argument("--run_type", help='Set which subroutines to call using keyword: [full, mesh, surfaces]. Default full', default='full')
parser.add_argument("--activation_type", help='Toggle whether to evaluate activation from root node locations or from existing map', default='rootnodes')
args = parser.parse_args()
name = args.name
refresh = args.refresh
visualise = args.visualise
if refresh in true:
    refresh = True
else:
    refresh = False
if visualise in true:
    visualise = True
else:
    visualise = False
root_dir = args.root_dir
input_format = args.input_format
input_dir = root_dir + args.input_dir
save_dir = root_dir + args.save_dir
output_dir = root_dir + args.output_dir
visualisation_format = args.visualisation_format
visualisation_dir = args.visualisation_dir
run_type = args.run_type
activation_type = args.activation_type
if not os.path.exists(input_dir):
    exit('Input directory does not exist: '+input_dir)
if not os.path.exists(save_dir):
    os.mkdir(save_dir)
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
if visualise:
    if not os.path.exists(visualisation_dir):
        os.mkdir(visualisation_dir)
print('##### Alya setup #####')
print('# Configuration:')
print('Run type:\t'+run_type)
print('Name:\t'+name)
print('Activation type:\t'+activation_type)
print('Input format:\t'+input_format)
print('Input directory:\t'+input_dir)
print('Save directory:\t'+save_dir)
print('Output directory:\t'+output_dir)
if visualise:
    print('Visualisation format:\t'+visualisation_format)
    print('Visualisation directory:\t'+visualisation_dir)
else:
    print('Visualisation:\tOFF')
print('#######################')

i = readfunctions.I(name, input_dir, input_format, save_dir, refresh)
o = writefunctions.O(i,name, output_dir)
v = visualisationfunctions.V(i, name, visualisation_dir, visualisation_format)


if run_type == 'mesh':
    i._read_mesh()
    o._write_alya_mesh()
    if visualise:
        if visualisation_format == 'ensight':
            v._ensight_export_mesh(i.nodes, i.elems)
            v._ensight_export_vector_per_node('TRANSV', i.transmural_vector)
            v._ensight_export_vector_per_node('LONGIV', i.longitudinal_vector)
            v._ensight_export_vector_per_node('ROTATV', i.rotational_vector)
            v._ensight_export_scalar_per_node('TRANS', i.transmural_coordinate)
            v._ensight_export_scalar_per_node('LONGI', i.longitudinal_coordinate)
            v._ensight_export_scalar_per_node('ROTAT', i.rotational_coordinate)
elif run_type == 'surfaces':
    i._read_mesh()
    i._read_surfaces()
    o._write_alya_surface_boundary()
    if visualise:
        if visualisation_format == 'ensight':
            v._ensight_export_scalar_per_node('FACELABEL', i.facesnodes_label)
elif run_type == 'celltype':
    _read_mesh()
    _read_wall_axes()
    _write_alya_celltype()
    if visualise:
        if visualisation_format == 'ensight':
            v._ensight_export_scalar_per_node('CELLTYPE', i.celltype)
elif run_type == 'activation':
    if activation_type == 'readin':
        _read_endocardial_activation()
    elif activation_type == 'rootnodes':
        _evaluate_dijkstra_endocardial_activation()
    else:
        exit('Activation type: '+activation_type+' not recognised')
    _write_alya_stimulus()
    if visualise:
        if visualisation_format == 'ensight':
            v._ensight_export_scalar_per_node('STIMULI', i.stimuli)
elif run_type == 'full':
    i._read_mesh()
    o._write_alya_mesh()
    if visualise:
        if visualisation_format == 'ensight':
            v._ensight_export_mesh(i.nodes, i.elems)
            v._ensight_export_vector_per_node('TRANSV', i.transmural_vector)
            v._ensight_export_vector_per_node('LONGIV', i.longitudinal_vector)
            v._ensight_export_vector_per_node('ROTATV', i.rotational_vector)
            v._ensight_export_scalar_per_node('TRANS', i.transmural_coordinate)
            v._ensight_export_scalar_per_node('LONGI', i.longitudinal_coordinate)
            v._ensight_export_scalar_per_node('ROTAT', i.rotational_coordinate)
    i._read_surfaces()
    o._write_alya_surface_boundary()
    if visualise:
        if visualisation_format == 'ensight':
            v._ensight_export_scalar_per_node('FACELABEL', i.facesnodes_label)
    o._write_alya_celltype()
    if visualise:
        if visualisation_format == 'ensight':
            v._ensight_export_scalar_per_node('CELLTYPE', i.celltype)
else:
    exit('Run type not found: '+run_type)
