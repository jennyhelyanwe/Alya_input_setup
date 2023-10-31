from meshpreprocessing import MeshPreprocessing
from generatefields import FieldGeneration
from alyaformat import AlyaFormat
from populationdrugtest import PopulationDrugTest
from postprocessing import PostProcessing
import numpy as np
import os

simulation_name = 'ho_LV_mesh'
workdir = os.getcwd()
meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/'+simulation_name+'/'+simulation_name+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False

vtk_dir = meta_data_dir + simulation_name + '/'

fields = FieldGeneration(name=simulation_name, geometric_data_dir=geometric_data_dir,
                         personalisation_data_dir='', verbose=verbose)
fields.generate_


