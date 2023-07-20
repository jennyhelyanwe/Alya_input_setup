from alyaformat import AlyaFormat
from postprocessing import PostProcessing
from sensitivityanalysis import SA
import os
import numpy as np

########################################################################################################################
# Global Settings
mesh_number = '05'
simulation_name = 'rodero_' + mesh_number + '_fine'
workdir = os.getcwd()
if 'icei' in workdir:
    system = 'jureca'
else:
    system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
clinical_data_dir = meta_data_dir + 'clinical_data/'
verbose = False

########################################################################################################################
# Step 1: Initialise Alya input files writing capabilities.
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/machine_learning_mechanics/'
elif system == 'heart':
    simulation_dir = './'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)

simulation_json_file = 'rodero_baseline_simulation_em_diastole_only.json'
# alya.do(simulation_json_file=simulation_json_file)
# quit()
#######################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
cell_parameter_names = np.array(['a_myocardium', 'end_diastole_p_lv', 'end_diastole_p_rv', 'pericardial_stiffness'])
baseline_parameter_values = np.array([6100,15000,0.1,1500000])
lower_bounds = [0.3050*10000, 0.5*10000, 0.5*10000, 10*10000]
upper_bounds = [6.1*10000, 3*10000, 3*10000, 200*10000]
baseline_json_file = simulation_json_file
baseline_dir = ''
simulation_dir = ''
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/machine_learning_mechanics/rodero_baseline_simulation_em_diastole_only_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/machine_learning_mechanics/'
elif system == 'heart':
    print('Not implemented for heart! ')
    quit()
sa = SA(name='sa', sampling_method='gridsampling', n=4, parameter_names=cell_parameter_names,
       baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
       simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
# quit()
sa.run_jobs(simulation_dir, start_id=0)
quit()
########################################################################################################################
# Step 3: Run Alya post-processing
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/machine_learning_mechanics'
elif system == 'heart':
    print('Not implemented for heart! ')
    quit()
# sa.run_jobs_postprocess(simulation_dir)
# quit()
########################################################################################################################
# Step 4: Downsample displacement data and save in appropriate output folder.

