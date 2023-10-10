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
# # Step 2: Use sampling methods to explore sensitivity analysis
# # Parameters to change: INaL, IKr, INaK, ICaL, Jrel, Jup, Ca50, kws, LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s,
# cell_parameter_names = np.array(['sf_gnal', 'sf_gkr', 'sf_gnak', 'sf_gcal', 'sf_jup', 'cal50', 'sfkws'])
# baseline_parameter_values = np.array([1,1,1,1,1,0.805,1])
# baseline_json_file = 'rodero_baseline_simulation_em.json'
# if system == 'jureca':
#     baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/'
# elif system == 'heart':
#     baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#     simulation_dir = sa_folder_name + '/'
# sa = SA(name='sa', sampling_method='saltelli', n=2 ** 4, parameter_names=cell_parameter_names,
#        baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#        simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
# # sa.setup()
# # sa.run_jobs(simulation_dir+sa_folder_name +'/')
# # quit()

#######################################################################################################################
# Step 2: Use sampling methods to explore sensitivity analysis
# Parameters to change: LVR, LVE, EDP, Tascaling, af, kepi, sigma_f, sigma_s, Kct
uq_cell_parameters = np.array([ 'cal50_myocardium', 'sf_gcal', 'sfkws_myocardium', 'sf_jup', 'sf_gnal'])
baseline_parameter_values = np.array([0.805, 1,1,1,1])
baseline_json_file = 'rodero_baseline_simulation_em.json'
if system == 'jureca':
    baseline_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/rodero_baseline_simulation_em_rodero_05_fine/'
elif system == 'heart':
    baseline_dir = '/users/jenang/Alya_setup_SA/rodero_baseline_simulation_em_rodero_05_fine/'
#
for i, parameter_name in enumerate(uq_cell_parameters):
    uq_folder_name = 'uq_' + parameter_name
    if system == 'jureca':
        simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + uq_folder_name + '/'
    elif system == 'heart':
        simulation_dir = uq_folder_name + '/'
    sa = SA(name='uq', sampling_method='range', n=3, parameter_names=np.array([parameter_name]),
            baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
            simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
    upper_bounds = baseline_parameter_values[i] * 2.0
    lower_bounds = baseline_parameter_values[i] * 0.5
    sa.setup(upper_bounds=upper_bounds, lower_bounds=lower_bounds)
    sa.run_jobs(simulation_dir=simulation_dir)
quit()

########################################################################################################################
# Step 3: Run Alya post-processing
# for i, parameter_name in enumerate(uq_cell_parameters):
#     uq_folder_name = 'uq_' + parameter_name
#     if system == 'jureca':
#         simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + uq_folder_name + '/'
#     elif system == 'heart':
#         simulation_dir = uq_folder_name + '/'
#     sa = SA(name='uq', sampling_method='range', n=3, parameter_names=np.array([parameter_name]),
#             baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
#             simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
#     sa.run_jobs_postprocess(simulation_dir=simulation_dir)
# quit()
########################################################################################################################
# Step 4: Evaluate QoIs and write out to results file
print('Plotting UQ curves')
for i, parameter_name in enumerate(uq_cell_parameters):
    print(parameter_name)
    uq_folder_name = 'uq_' + parameter_name
    if system == 'jureca':
        simulation_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/' + uq_folder_name + '/'
    elif system == 'heart':
        simulation_dir = uq_folder_name + '/'
    sa = SA(name='uq', sampling_method='range', n=3, parameter_names=np.array([parameter_name]),
            baseline_parameter_values=baseline_parameter_values, baseline_json_file=baseline_json_file,
            simulation_dir=simulation_dir, alya_format=alya, baseline_dir=baseline_dir, verbose=verbose)
    # sa.evaluate_qois(alya=alya, beat=0, qoi_save_dir=simulation_dir)
    sa.visualise_uq(alya=alya, beat=1, parameter_name=parameter_name)
    # sa.analyse(simulation_dir+'all_qois.csv')
quit()

########################################################################################################################
# Step 5: Evaluate Sobol indices and plot results
sa_figures_directory = simulation_dir
sa.analyse(sa_figures_directory)