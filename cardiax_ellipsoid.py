from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from alyaformat import AlyaFormat, run_job, run_job_postprocess
from populationdrugtest import PopulationDrugTest
from postprocessing import PostProcessing
import numpy as np
import os

simulation_name = 'ho_LV_mesh'
workdir = os.getcwd()
meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/'+simulation_name+'/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False

vtk_dir = geometric_data_dir
vtk_name = 'ho_LV_mesh'
#######################################################################################################################
# Step 1: Save input mesh into CSV format, as prescribed in myformat.py
meshpreprocessing = True
if meshpreprocessing:
    mesh = MeshPreprocessing(vtk_name=vtk_name, name=simulation_name, input_dir=vtk_dir, geometric_data_dir=geometric_data_dir,
                   verbose=verbose)
    mesh.geometry.lv_endocardium = 1 # Change from biventricular default
    mesh.geometry.epicardium = 2 # Change from biventricular default
    mesh.geometry.lid = 3 # Change form biventricular default
    mesh.read_geometry_from_vtk_cardiax_ellipsoid(save=False)
    mesh.generate_boundary_data_cardiax_ellipsoid(save=True)

########################################################################################################################
# Step 2: Generate fields for Alya simulation
fieldgeneration = False
if fieldgeneration:
    fields = FieldGeneration(name=simulation_name, geometric_data_dir=geometric_data_dir,
                             personalisation_data_dir='', verbose=verbose)
    fields.read_cardiax_ellipsoid_fibre_fields_vtk(fibre_vtk_filename=vtk_dir+'ho_LV_mesh_fibre.vtk')
    fields.generate_prestress_field_cardiax()
    fields.read_cavity_landmarks(save=True)

########################################################################################################################
# Step 5: Write Alya input files according to simulation protocol saved in .json file.
alyasetup = True
simulation_dir = ''
personalisation_data_dir = ''
system='jureca'
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_dir = './'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)
alya.geometry.lv_endocardium = 1 # Change from biventricular default
alya.geometry.epicardium = 2 # Change from biventricular default
alya.geometry.lid = 3 # Change form biventricular default
if alyasetup:
    simulation_json_file = 'cardiax_ellipsoid_mech.json'
    alya.do(simulation_json_file=simulation_json_file)

run = False
if run:
    run_job(alya.output_dir)