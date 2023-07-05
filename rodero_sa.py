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
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    simulation_dir = './'
personalisation_data_dir = meta_data_dir + 'results/personalisation_data/rodero_'+mesh_number+'/'
alya = AlyaFormat(name=simulation_name, geometric_data_dir=geometric_data_dir,
                  personalisation_dir=personalisation_data_dir, clinical_data_dir=clinical_data_dir,
                  simulation_dir = simulation_dir, verbose=verbose)

########################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: INaL, IKr, INaK, ICaL, Jrel, Jup, Ca50, kws, LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s,
cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
baseline_parameter_values = np.array([1,1,1,1,1,0.805,1])
baseline_json_file = 'rodero_baseline_simulation_em.json'
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = 'sensitivity_analyses/'
sa = SA(name='sa', sampling_method='saltelli', n=2 ** 2, parameter_names=cell_parameter_names,
       baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
       simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup()
# sa.run_jobs(simulation_dir+'sensitivity_analyses/')
# quit()

#######################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
baseline_parameter_values = np.array([1,1,1,1,1,0.805,1])
baseline_json_file = 'rodero_baseline_simulation_em.json'
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/sensitivity_analyses_256_samples/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
    simulation_dir = 'sensitivity_analyses/'
sa = SA(name='sa', sampling_method='saltelli', n=2 ** 4, parameter_names=cell_parameter_names,
       baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
       simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# sa.setup()
sa.run_jobs(simulation_dir+'sensitivity_analyses/')
quit()
########################################################################################################################
# Step 3: Run Alya post-processing
if system == 'jureca':
    simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/sensitivity_analyses/'
elif system == 'heart':
    simulation_dir = '/users/jenang/Alya_setup_SA/alya_csv/rodero_' + mesh_number + '/'
# sa.run_jobs_postprocess(simulation_dir)

########################################################################################################################
# Step 4: Evaluate QoIs and write out to results file
sa.evaluate_qois(alya=alya, beat=0, qoi_save_dir=simulation_dir)
sa.analyse()
quit()
########################################################################################################################
# Step 5: Evaluate Sobol indices and plot results
sa_figures_directory = simulation_dir
sa.analyse(sa_figures_directory)