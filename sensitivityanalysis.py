import scipy
from SALib.sample import saltelli
import numpy as np
import json
import os
from alyaformat import AlyaFormat

class SA:
    def __init__(self, name, sampling_method, n, parameter_names, baseline_parameter_values, baseline_json_file,
                 simulation_dir, alya_format, verbose):
        self.name = name
        self.sampling_method = sampling_method
        self.n = n # This is a parameter that could either be the number of samples (for LHS) or is related to number of samples (Saltelli)
        self.parameter_names = parameter_names
        self.baseline_parameter_values = baseline_parameter_values
        self.number_of_parameters = parameter_names.shape[0]
        self.parameter_set = None
        self.simulation_dir = simulation_dir
        self.alya_format = alya_format
        self.verbose = verbose
        self.all_simulation_dirs = []
        if not os.path.exists(self.simulation_dir):
            os.mkdir(self.simulation_dir)
        self.generate_parameter_set()
        self.generate_alya_simulation_json(baseline_json_file)


    def generate_parameter_set(self):
        upper_bounds = self.baseline_parameter_values * 2.0
        lower_bounds = self.baseline_parameter_values * 0.5
        if self.sampling_method == 'lhs':
            sampler = scipy.stats.qms.LatinHypercube(d=self.number_of_parameters)
            sample = sampler.random(n=self.n)
            self.parameter_set = scipy.qmc.scale(sample, lower_bounds, upper_bounds)
        elif self.sampling_method == 'saltelli':
            ranges = np.vstack((lower_bounds, upper_bounds)).transpose()
            problem = {
                'num_vars' : self.number_of_parameters,
                'names' : self.parameter_names,
                'bounds' : ranges
            }
            self.parameter_set = saltelli.sample(problem, self.n)

    def generate_alya_simulation_json(self, baseline_json_file):
        baseline_simulation_dict = json.load(open(baseline_json_file, 'r'))
        for sample_i in range(self.parameter_set.shape[0]):
            sample_simulation_dict = baseline_simulation_dict
            for variable_i in range(self.parameter_names.shape[0]):
                sample_simulation_dict[self.parameter_names[variable_i]][0] = self.parameter_set[sample_i, variable_i]
            with open(self.simulation_dir + self.name + '_' + str(sample_i) + '.json', 'w') as f:
                json.dump(sample_simulation_dict, f)
            self.alya_format.output_dir = self.simulation_dir + self.name + '_' + str(sample_i) + '/'
            print('Writing to: ' + self.alya_format.output_dir)
            self.alya_format.do(simulation_json_file=self.simulation_dir + self.name + '_' + str(sample_i) + '.json')
            self.all_simulation_dirs.append(self.alya_format.output_dir)

    def run(self):
        os.system('cd '+self.simulation_dir)
        for simulation_i in range(len(self.all_simulation_dirs)):
            os.system('cd '+self.all_simulation_dirs[simulation_i])
            os.system('sbatch run_job.cmd')
            os.system('cd ../')



