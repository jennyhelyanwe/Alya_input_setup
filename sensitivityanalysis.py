import scipy
from SALib.sample import saltelli
from SALib import ProblemSpec
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
from postprocessing import PostProcessing
import pandas as pd
class SA:
    def __init__(self, name, sampling_method, n, parameter_names, baseline_parameter_values, baseline_json_file, baseline_dir,
                 simulation_dir, alya_format, verbose):
        self.name = name
        self.sampling_method = sampling_method
        self.n = n # This is a parameter that could either be the number of samples (for LHS) or is related to number of samples (Saltelli)
        self.parameter_names = parameter_names
        self.baseline_parameter_values = baseline_parameter_values
        self.baseline_dir = baseline_dir
        self.number_of_parameters = parameter_names.shape[0]
        self.parameter_set = None
        self.simulation_dir = simulation_dir
        self.alya_format = alya_format
        self.verbose = verbose
        self.all_simulation_dirs = []
        if not os.path.exists(self.simulation_dir):
            os.mkdir(self.simulation_dir)
        self.baseline_json_file = baseline_json_file

    def setup(self):
        self.generate_parameter_set()
        self.generate_alya_simulation_json(self.baseline_json_file)

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
            self.parameter_set = saltelli.sample(problem, self.n, calc_second_order=True)

    def generate_alya_simulation_json(self, baseline_json_file):
        baseline_simulation_dict = json.load(open(baseline_json_file, 'r'))
        for sample_i in range(self.parameter_set.shape[0]):
            sample_simulation_dict = baseline_simulation_dict
            for variable_i in range(self.parameter_names.shape[0]):
                if 'sf_' in self.parameter_names[variable_i]:
                    sample_simulation_dict[self.parameter_names[variable_i]][0][0] = self.parameter_set[
                        sample_i, variable_i]
                    sample_simulation_dict[self.parameter_names[variable_i]][0][1] = self.parameter_set[
                        sample_i, variable_i]
                    sample_simulation_dict[self.parameter_names[variable_i]][0][2] = self.parameter_set[
                        sample_i, variable_i]
                else:
                    sample_simulation_dict[self.parameter_names[variable_i]][0] = self.parameter_set[sample_i, variable_i]
            with open(self.simulation_dir + self.name + '_' + str(sample_i) + '.json', 'w') as f:
                json.dump(sample_simulation_dict, f)
            #self.alya_format.output_dir = self.simulation_dir + self.name + '_' + str(sample_i) + '/'
            self.alya_format.simulation_dir = self.simulation_dir
            self.alya_format.do(simulation_json_file=self.simulation_dir+self.name + '_' + str(sample_i) + '.json',
                                SA_flag=True, baseline_dir=self.baseline_dir)
            self.all_simulation_dirs.append(self.alya_format.output_dir)
        with open(self.simulation_dir+'/all_simulation_dirs.txt', 'w') as f:
            for i in range(len(self.all_simulation_dirs)):
                f.write(self.all_simulation_dirs[i]+'\n')
        # baseline_simulation_dict = json.load(open(baseline_json_file, 'r'))
        # threadsNum = multiprocessing.cpu_count()
        # threadsNum = 10
        # all_simulation_dirs = pymp.shared.array(self.parameter_set.shape[0])
        # parameter_names = pymp.shared.array(self.parameter_names.shape)
        # parameter_names = self.parameter_names
        # parameter_set = pymp.shared.array(self.parameter_set.shape)
        # parameter_set = self.parameter_set
        # manager = multiprocessing.Manager()
        # sample_simulation_dict = manager.dict()
        # sample_simulation_dict = baseline_simulation_dict
        # with pymp.Parallel(min(threadsNum, self.parameter_set.shape[0])) as p1:
        #     for sample_i in p1.range(parameter_set.shape[0]):
        #         for variable_i in range(parameter_names.shape[0]):
        #             if 'sf_' in self.parameter_names[variable_i]:
        #                 sample_simulation_dict[self.parameter_names[variable_i]][0][0] = self.parameter_set[
        #                                 sample_i, variable_i]
        #                 sample_simulation_dict[self.parameter_names[variable_i]][0][1] = self.parameter_set[
        #                     sample_i, variable_i]
        #                 sample_simulation_dict[self.parameter_names[variable_i]][0][2] = self.parameter_set[
        #                     sample_i, variable_i]
        #             else:
        #                 sample_simulation_dict[parameter_names[variable_i]][0] = parameter_set[sample_i, variable_i]
        #         with open(self.simulation_dir + self.name + '_' + str(sample_i) + '.json', 'w') as f:
        #             json.dump(sample_simulation_dict, f)
        #         self.alya_format.output_dir = self.simulation_dir + self.name + '_' + str(sample_i) + '/'
        #         print('Writing to: ' + self.alya_format.output_dir)
        #         self.alya_format.do(simulation_json_file=self.simulation_dir + self.name + '_' + str(sample_i) + '.json', parallel_flag=True)
        #         print('finished alya do! ')
        #         all_simulation_dirs[sample_i] = self.alya_format.output_dir

    def evaluate_qois(self, alya, beat, qoi_save_dir):
        with open(self.simulation_dir+'/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        finished_simulation_dirs = []
        for simulation_i in range(len(all_simulation_dirs)):
            with open(all_simulation_dirs[simulation_i].split()[0]+'/heart.post.alyafil', 'r') as f:
                if len(f.readlines()) == 1000:
                    finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])
        all_simulation_dirs_shared = pymp.shared.array((len(all_simulation_dirs)), dtype=str)
        threadsNum = multiprocessing.cpu_count()
        # with pymp.Parallel(min(threadsNum, len(all_simulation_dirs))) as p1:
        #     for simulation_i in p1.range(len(all_simulation_dirs)):
        if False:
            for simulation_i in range(len(finished_simulation_dirs)):
                alya_output_dir = finished_simulation_dirs[simulation_i]
                json_file = self.simulation_dir + 'sa_' + str(simulation_i) + '.json'
                if os.path.exists(alya_output_dir+'/results_csv/'):
                    post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                          alya_output_dir=alya_output_dir, verbose=self.verbose)
                    post.evaluate_ecg_pv_biomarkers(beat=beat)
                    post.save_qoi(filename=qoi_save_dir+'qoi_'+str(simulation_i)+'.csv')
        qois = pd.DataFrame({})
        for simulation_i in range(len(finished_simulation_dirs)):
            filename = qoi_save_dir+'qoi_'+str(simulation_i)+'.csv'
            qoi = json.load(open(filename, 'r'))
            qois = qois.append(pd.DataFrame(qoi, index=[simulation_i]))
        qois.to_csv(qoi_save_dir+'all_qois.csv')
        self.qois_db =qois


    def analyse(self, filename):
        self.qois_db = pd.read_csv(filename, index_col=False)
        names = self.parameter_names
        # selected_qois = self.qois_db.columns.values.tolist() # All QoIs
        selected_qois = ['t_pe_mean', 't_peak_mean', 'qt_dur_mean', 't_polarity_mean', 'EDVL', 'LVEF', 'PmaxL', 'SVL']
        sp = ProblemSpec({'num_vars': len(names),
                          'names': names,
                          'bounds': [[0.5, 2]] * len(names),
                          'outputs': selected_qois
                          })
        qoi_names = sp.get('outputs')
        Y = self.qois_db[selected_qois].values
        ## Need to get parameter values!
        X = np.zeros((Y.shape[0], len(names)))
        for simulation_i in range(Y.shape[0]):
            filename = self.simulation_dir+'sa_'+str(simulation_i)+'.json'
            dict = json.load(open(filename, 'r'))
            for param_i in range(len(names)):
                if 'sf_' in names[param_i]:
                    X[simulation_i,param_i] = dict[names[param_i]][0][0]
                else:
                    X[simulation_i,param_i] = dict[names[param_i]][0]

        # Scatter plots with correlation coefficients
        if len(Y.shape) == 1:
            num_qois = 1
        else:
            num_qois = Y.shape[1]

        ################################################################################################################
        # fig = plt.figure(tight_layout=True, figsize=(18, 10))
        # gs = GridSpec(num_qois, X.shape[1])
        # for qoi_i in range(num_qois):
        #     for param_j in range(X.shape[1]):
        #         ax = fig.add_subplot(gs[qoi_i, param_j])
        #         x = X[:,param_j]
        #         if num_qois == 1:
        #             y = Y
        #         else:
        #             y = Y[:,qoi_i]
        #         sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s':1})
        #         ax.text(x=np.amin(x), y=np.amax(y), va='top', ha='left',
        #                                   s='p=%.2f' % (np.corrcoef(x,y)[0,1]))
        #         if qoi_i == num_qois-1:
        #             ax.set_xlabel(names[param_j])
        #         if param_j == 0:
        #             ax.set_ylabel(qoi_names[qoi_i])
        # plt.savefig('scatter_qois_vs_parameters.png')
        ################################################################################################################
        # fig = plt.figure(tight_layout=True, figsize=(18, 10))
        # gs = GridSpec(num_qois, num_qois)
        # for qoi_i in range(num_qois):
        #     for qoi_j in range(num_qois):
        #         if qoi_i >= qoi_j:
        #             ax = fig.add_subplot(gs[qoi_i, qoi_j])
        #             if num_qois == 1:
        #                 y = Y
        #             else:
        #                 y = Y[:, qoi_i]
        #             x = Y[:, qoi_j]
        #             sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s': 1})
        #             ax.text(x=np.amin(x), y=np.amax(y), va='top', ha='left',
        #                     s='p=%.2f' % (np.corrcoef(x, y)[0, 1]))
        #         if qoi_i == num_qois - 1:
        #             ax.set_xlabel(qoi_names[qoi_j])
        #         if qoi_j == 0:
        #             ax.set_ylabel(qoi_names[qoi_i])
        # plt.show()
        # plt.savefig('scatter_qois_vs_qois.png')
        ####################################################################################################################
        # Tornado plot of sensitivity indices https://seaborn.pydata.org/examples/part_whole_bars.html
        sp.set_results(Y)
        Si = sp.analyze_sobol(print_to_console=False, calc_second_order=False)
        fig = plt.figure(tight_layout=True, figsize=(15,6))
        gs = GridSpec(1, num_qois)
        data = Si.to_df()
        sns.set_theme(style="whitegrid")
        for qoi_i in range(num_qois):
            ax = fig.add_subplot(gs[0, qoi_i])
            st_s1_data = pd.concat([data[qoi_i][0], data[qoi_i][1]], axis=1)
            sorted_data = st_s1_data.reindex(st_s1_data.abs().sort_values('ST', ascending=False).index)
            names = []
            for row in sorted_data.index:
                names.append(row)
            sns.set_color_codes("pastel")
            sns.barplot(data=sorted_data, x='ST', y=names, label='ST', color='b')
            sns.set_color_codes("muted")
            sns.barplot(data=sorted_data, x='S1', y=names, label='S1', color='b')

            ax.set( ylabel="",
                   xlabel=qoi_names[qoi_i])
            if qoi_i == num_qois-1:
                ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', frameon=True)
        plt.show()
        plt.savefig('ST_S1_tornado.png')



    def run_jobs(self, simulation_dir):
        with open(simulation_dir+'/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        for simulation_i in range(len(all_simulation_dirs)):
            cmd = 'cd '+all_simulation_dirs[simulation_i].split()[0]
            os.system('cd '+all_simulation_dirs[simulation_i].split()[0]+'; pwd ; sbatch run_job.cmd')


    def run_jobs_postprocess(self, simulation_dir):
        with open(simulation_dir + '/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        for simulation_i in range(len(all_simulation_dirs)):
            cmd = 'cd ' + all_simulation_dirs[simulation_i].split()[0]
            os.system('cd ' + all_simulation_dirs[simulation_i].split()[0] + '; pwd ; sbatch run_job_postprocess.cmd')

        

