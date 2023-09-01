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
        self.population_sf_filenames = []
        for population_i in range(self.population_size):
            data = pd.read_csv(self.personalised_population_dir + self.name + '_fine_nodefield_personalisation-biomarker_' + str(
                population_i) + '.csv')
            for sf_name in population_sf_names:
                filename = self.personalised_population_dir + self.name + '_populationid_' + str(population_i) + '.' + sf_name
                self.population_sf_filenames.append(filename)
                with open(filename, 'w') as f:
                    for i in range(data[sf_name].values.shape[0]):
                        f.write(str(data[sf_name].index[i] + 1) + '\t' + str(data[sf_name].values[i]) + '\n')
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
                current_simulation_dir = dose_dir + 'populationid_'+str(population_i)+'_'+drug_name+'_dose_' + str(dose_i)
                with open(json_filename, 'w') as f:
                    json.dump(simulation_dict, f)
                if not os.path.exists(current_simulation_dir):
                    os.mkdir(current_simulation_dir)
                self.alya_format.simulation_dir = dose_dir
                print('Writing out Alya file to '+current_simulation_dir)
                self.alya_format.do(simulation_json_file=json_filename,
                                    drug_flag=True, baseline_dir=baseline_dir)
                # Retrieve scaling factor for this member of the population
                for sf_name in self.population_sf_names:
                    os.system('cp '+self.personalised_population_dir+self.name+'_populationid_'+str(population_i)+'.' + sf_name + ' '+current_simulation_dir + 'heart.'+sf_name)
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
