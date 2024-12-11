from alyaformat import AlyaFormat
from postprocessing import PostProcessing
from sensitivityanalysis_uncertaintyquantification import SAUQ
from generatefields import FieldGeneration
import os
import numpy as np
import json

########################################################################################################################
# Global Settings
mesh_number = '3D'
simulation_name = 'cube_free_contraction'
verbose = False

## Set up Alya format
simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/sensitivity_analysis_cube_free_contraction/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=None,
                  personalisation_dir=None, clinical_data_dir=None,
                  simulation_dir=simulation_dir, verbose=verbose)

## Sample ICaL, Jup, Ca50, kws, and INaL parameters
baseline_json_file = './cube_free_contraction_em.json'
simulation_dict = json.load(open(baseline_json_file, 'r'))
parameter_names = np.array(['sf_gcal', 'sf_jup', 'cal50_myocardium', 'sfkws_myocardium', 'sf_gnal'])
baseline_parameter_values = np.array([simulation_dict['sf_gcal'][0][0], simulation_dict['sf_jup'][0][0], simulation_dict['cal50'][0],
                                      simulation_dict['sfkws'][0], simulation_dict['sf_gnal'][0][0]])
print(baseline_parameter_values)
upper_bounds = baseline_parameter_values * 2.0
lower_bounds = baseline_parameter_values * 0.5
baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/free_contraction/baseline/'
simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/sensitivity_analysis_cube_free_contraction/'
sa = SAUQ(name='sa', sampling_method='saltelli', n=2 ** 3, parameter_names=parameter_names,
          baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
          simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# sa.run_jobs(simulation_dir=simulation_dir, start_id=0)

# sa.run_jobs_postprocess(simulation_dir=simulation_dir)

#
beat = 1
print('Sorting simulations...')
sa.sort_simulations(tag='cube_postprocess')
print('Evaluating x displacements...')
cube_deformation_ta_post = sa.evaluate_qois(qoi_group_name='cube_deformation_ta', alya=alya, beat=beat,
                                         qoi_save_dir=simulation_dir, analysis_type='sa')
print('Visualising sa...')
sa.visualise_sa(beat=beat, cube_deformation_ta_post=cube_deformation_ta_post)

sa.analyse(filename=simulation_dir+'cube_deformation_qois.csv', qois=['peak_displ_x', 'rise_dxdt', 'decay_dxdt', 'peak_ta', 'rise_dtadt', 'decay_dtadt'])