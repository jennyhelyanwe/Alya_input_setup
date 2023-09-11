import numpy as np
import json
import os
from alyaformat import AlyaFormat
import pymp, multiprocessing
import matplotlib
matplotlib.use('tkagg')
from matplotlib.gridspec import GridSpec
from matplotlib import pyplot as plt
import seaborn as sns
from postprocessing import PostProcessing, mapIndices
import pandas as pd
from ECGPV_visualisation import ECGPV_visualisation

class PopulationDrugTest:
    def __init__(self, name, personalised_population_dir, population_size, alya_format, verbose):
        self.name = name
        self.population_size = population_size
        self.personalised_population_dir = personalised_population_dir

        self.alya_format = alya_format
        self.verbose = verbose
        self.all_simulation_dirs = []

    def setup_drug_test(self, drug_name, drug_doses, simulation_dir, population_sf_names,  baseline_json_file, baseline_dir):
        if not os.path.exists(simulation_dir):
            os.mkdir(simulation_dir)
        simulation_dict = json.load(open(baseline_json_file, 'r'))
        # self.population_sf_filenames = []
        # for population_i in range(self.population_size):
        #     data = pd.read_csv(self.personalised_population_dir + self.name + '_fine_nodefield_personalisation-biomarker_' + str(
        #         population_i) + '.csv')
        #     for sf_name in population_sf_names:
        #         filename = self.personalised_population_dir + self.name + '_populationid_' + str(population_i) + '.' + sf_name
        #         self.population_sf_filenames.append(filename)
        #         with open(filename, 'w') as f:
        #             for i in range(data[sf_name].values.shape[0]):
        #                 f.write(str(data[sf_name].index[i] + 1) + '\t' + str(data[sf_name].values[i]) + '\n')
        ## Set up Alya simulations for baseline population
        baseline_population_dir = simulation_dir + 'baseline_population/'
        if not os.path.exists(baseline_population_dir):
            os.mkdir(baseline_population_dir)
        for population_i in range(self.population_size):
            json_filename = baseline_population_dir + 'populationid_' + str(population_i) + '_baseline.json'
            with open(json_filename, 'w') as f:
                json.dump(simulation_dict, f)
            self.alya_format.simulation_dir = baseline_population_dir
            self.alya_format.do(simulation_json_file=json_filename, drug_flag=True, baseline_dir=baseline_dir)
            self.all_simulation_dirs.append(self.alya_format.output_dir)
        ## Set up Alya simulations for drug doses
        for dose_i in range(len(drug_doses)):
            # For each dose of the drug run a population of models
            drug_dict = drug_doses[dose_i]
            dose_dir = simulation_dir+'dose_'+str(dose_i) + '/'
            if not os.path.exists(dose_dir):
                os.mkdir(dose_dir)
            with open(dose_dir + 'drug_sf.json', 'w') as f:
                json.dump(drug_dict, f)
            sf_keys = list(drug_dict.keys())
            for population_i in range(self.population_size):
                for sf_key in sf_keys:
                    assert 'sf_' in sf_key
                    simulation_dict[sf_key][0][0] = drug_dict[sf_key]
                    simulation_dict[sf_key][0][1] = drug_dict[sf_key]
                    simulation_dict[sf_key][0][2] = drug_dict[sf_key]
                json_filename = dose_dir + 'populationid_'+str(population_i)+'_'+drug_name+'_dose_' + str(dose_i) + '.json'
                with open(json_filename, 'w') as f:
                    json.dump(simulation_dict, f)
                self.alya_format.simulation_dir = dose_dir
                # print('Writing out Alya file to '+alya_format.output_dir)
                self.alya_format.do(simulation_json_file=json_filename,
                                    drug_flag=True, baseline_dir=baseline_dir)
                # Retrieve scaling factor for this member of the population
                for sf_name in population_sf_names:
                    cmd = 'cp '+self.personalised_population_dir+self.name+'_populationid_'+str(population_i)+'.' + sf_name + ' '+self.alya_format.output_dir+'heart.'+sf_name
                    os.system(cmd)
                self.all_simulation_dirs.append(self.alya_format.output_dir)
        with open(simulation_dir+'/all_simulation_dirs.txt', 'w') as f:
            for i in range(len(self.all_simulation_dirs)):
                f.write(self.all_simulation_dirs[i]+'\n')


    def run_jobs(self, simulation_dir, start_id=0):
        with open(simulation_dir+'/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        for simulation_i in range(start_id, len(all_simulation_dirs)):
            cmd = 'cd '+all_simulation_dirs[simulation_i].split()[0]
            os.system('cd '+all_simulation_dirs[simulation_i].split()[0]+'; pwd ; sbatch run_job.cmd')


    def run_jobs_postprocess(self, simulation_dir):
        with open(simulation_dir + '/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        for simulation_i in range(len(all_simulation_dirs)):
            cmd = 'cd ' + all_simulation_dirs[simulation_i].split()[0]
            os.system('cd ' + all_simulation_dirs[simulation_i].split()[0] + '; pwd ; sbatch run_job_postprocess.cmd')


    def visualise_drug_effect_ecgs(self, beat, drug_name, drug_doses, simulation_dir):
        CL = 1.0
        vis = ECGPV_visualisation(CL)
        fig = plt.figure(tight_layout=True, figsize=(20, 10))
        gs = GridSpec(2,4)
        axes = []
        for i in [0,1]:
            for j in [0,1,2,3]:
                axes.append(fig.add_subplot(gs[i,j]))
        # fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        nb_leads = 8
        # t_padding = 0.2
        # v_padding = 0.5
        # scale = 2
        # fig_size = [np.ceil(CL*2+t_padding)*0.5/0.2*scale, 6*scale]
        # fig = plt.figure(tight_layout=True, figsize=fig_size)
        # gs = GridSpec(1,1)
        # ax = fig.add_subplot(gs[2,4])

        lead_names = ['Is', 'IIs', 'V1s', 'V2s', 'V3s', 'V4s', 'V5s', 'V6s']

        # baseline_dir ='/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/DTI004_baseline_simulation_ep_DTI004_fine/'
        # clinical_ecg = np.loadtxt(baseline_dir + 'DTI004_clinical_full_ecg.txt', delimiter=',')
        # max_leads = np.amax(clinical_ecg)
        # clinical_ecg = clinical_ecg[:, 100:] / max_leads
        # Plot baseline popluation
        for population_i in range(self.population_size):
            current_sim_dir = simulation_dir + 'baseline_population/'
            if os.path.exists(current_sim_dir + '/heart.exm.vin'):
                ecgs = vis._read_ECG(current_sim_dir + '/heart')
                for lead_i in range(nb_leads):
                    axes[lead_i].plot(ecgs['ts'][beat - 1],
                                              ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_all_leads'],
                                              'k', linewidth=0.5, label='Baseline')
                else:
                    axes[lead_i].plot(ecgs['ts'][beat - 1],
                                      ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_all_leads'],
                                      'k', linewidth=0.5)

        # colours = ['g', 'b', 'r', 'm', 'c', 'y']
        colours = ['#fde725', '#90d743', '#35b779', '#21918c', '#31688e', '#443983', '#440154']
        for dose_i in range(len(drug_doses)):
            # For each dose of the drug run a population of models
            drug_dict = drug_doses[dose_i]
            dose_dir = simulation_dir+'dose_'+str(dose_i) + '/'
            if dose_i == 0:
                lines = []
                legends = []
            for population_i in range(self.population_size):
                self.alya_format.simulation_dir = dose_dir
                current_sim_dir = dose_dir + 'populationid_'+str(population_i) + '_' + drug_name + \
                                  '_dose_' + str(dose_i) + '_' + self.name + '/'

                if os.path.exists(current_sim_dir+'/heart.exm.vin'):
                    ecgs = vis._read_ECG(current_sim_dir + '/heart')
                    for lead_i in range(nb_leads):
                        if population_i == 0:
                            axes[lead_i].plot(ecgs['ts'][beat - 1],
                                              ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_all_leads'],
                                              colours[dose_i], linewidth=0.5, label='Dose :'+str(dose_i))
                        else:
                            axes[lead_i].plot(ecgs['ts'][beat-1], ecgs[lead_names[lead_i]][beat-1]/ecgs['max_all_leads'],
                                          colours[dose_i], linewidth=0.5)
                        axes[lead_i].set_title(lead_names[lead_i], fontsize=20)
                        axes[lead_i].set_ylim([-1.5, 1.5])
                        for tick in axes[lead_i].xaxis.get_major_ticks():
                            tick.label1.set_fontsize(14)
                        for tick in axes[lead_i].yaxis.get_major_ticks():
                            tick.label1.set_fontsize(14)
                        # if (dose_i==0) & (lead_i == 0) & (population_i == 0):
                        #     lines.append(line)
                        #     legends.append('Dose: '+str(dose_i))
        axes[7].legend()
        plt.show()

    def evaluate_drug_effect_ecgs(self, beat, drug_name, drug_doses, simulation_dir):
        CL = 1.0
        vis = ECGPV_visualisation(CL)
        for dose_i in range(len(drug_doses)):
            # For each dose of the drug run a population of models
            drug_dict = drug_doses[dose_i]
            dose_dir = simulation_dir+'dose_'+str(dose_i) + '/'
            if dose_i == 0:
                lines = []
                legends = []
            for population_i in range(self.population_size):
                self.alya_format.simulation_dir = dose_dir
                current_sim_dir = dose_dir + 'populationid_'+str(population_i) + '_' + drug_name + \
                                  '_dose_' + str(dose_i) + '_' + self.name + '/'
                if os.path.exists(current_sim_dir+'/heart.exm.vin'):
                    ecgs = vis._read_ECG(current_sim_dir + '/heart')

