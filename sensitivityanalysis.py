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
from postprocessing import PostProcessing, mapIndices
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

    def setup(self, upper_bounds, lower_bounds):
        self.generate_parameter_set( upper_bounds, lower_bounds)
        print('Number of samples/simulations to set up: ', str(self.parameter_set.shape[0]))
        print('Number of parameters included: ', str(self.parameter_set.shape[1]))
        self.generate_alya_simulation_json(self.baseline_json_file)

    def generate_parameter_set(self, upper_bounds, lower_bounds):
        # upper_bounds = self.baseline_parameter_values * 2.0
        # lower_bounds = self.baseline_parameter_values * 0.5
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
        elif self.sampling_method == 'gridsampling':
            self.parameter_set = np.zeros((self.n**self.number_of_parameters, self.number_of_parameters))
            a = np.linspace(lower_bounds[0], upper_bounds[0], self.n)
            b = np.linspace(lower_bounds[1], upper_bounds[1], self.n)
            c = np.linspace(lower_bounds[2], upper_bounds[2], self.n)
            d = np.linspace(lower_bounds[3], upper_bounds[3], self.n)
            i = 0
            for a_i in a:
                for b_i in b:
                    for c_i in c:
                        for d_i in d:
                            self.parameter_set[i, :] = [a_i, b_i, c_i, d_i]
                            i = i + 1

        elif self.sampling_method == 'range':
            self.parameter_set = np.zeros((3, 1))
            self.parameter_set[0, 0] = lower_bounds
            self.parameter_set[2, 0] = upper_bounds
            self.parameter_set[1, 0] = (lower_bounds + upper_bounds) / 2

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
                elif '_lv' in self.parameter_names[variable_i]:
                    variable_name = self.parameter_names[variable_i].split('_lv')[0]
                    sample_simulation_dict[variable_name][0] = self.parameter_set[sample_i, variable_i]
                elif '_rv' in self.parameter_names[variable_i]:
                    variable_name = self.parameter_names[variable_i].split('_rv')[0]
                    sample_simulation_dict[variable_name][1] = self.parameter_set[sample_i, variable_i]
                elif '_myocardium' in self.parameter_names[variable_i]:
                    variable_name = self.parameter_names[variable_i].split('_myocardium')[0]
                    sample_simulation_dict[variable_name][0] = self.parameter_set[sample_i, variable_i]
                elif '_valveplug' in self.parameter_names[variable_i]:
                    variable_name = self.parameter_names[variable_i].split('_valveplug')[0]
                    sample_simulation_dict[variable_name][0] = self.parameter_set[sample_i, variable_i]
                else:
                    sample_simulation_dict[self.parameter_names[variable_i]]= self.parameter_set[sample_i, variable_i]
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

    def evaluate_qois(self, qoi_group_name, alya, beat, qoi_save_dir, analysis_type):
        postp_objects = []
        for simulation_i in range(len(self.finished_simulation_dirs)):
            alya_output_dir = self.finished_simulation_dirs[simulation_i]
            json_file = self.simulation_dir + analysis_type + '_' + str(simulation_i) + '.json'
            if os.path.exists(alya_output_dir + '/results_csv/timeset_1.csv'):
                post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                      alya_output_dir=alya_output_dir, verbose=self.verbose)
                if qoi_group_name == 'ecg':
                    post.evaluate_ecg_biomarkers(beat=beat)
                    post.save_qoi(filename=qoi_save_dir + 'ecg_qoi_' + str(simulation_i) + '.csv')
                elif qoi_group_name == 'pv':
                    post.evaluate_pv_biomarkers(beat=beat)
                    post.save_qoi(filename=qoi_save_dir + 'pv_qoi_' + str(simulation_i) + '.csv')
                elif qoi_group_name == 'deformation':
                    post.evaluate_deformation_biomarkers(beat=beat)
                    post.save_qoi(filename=qoi_save_dir + 'deformation_qoi_' + str(simulation_i) + '.csv')
                elif qoi_group_name == 'fibre_work':
                    post.evaluate_fibre_work_biomarkers(beat=beat)
                    post.save_qoi(filename=qoi_save_dir + 'fibrework_qoi_' + str(simulation_i) + '.csv')
                postp_objects.append(post)

        # Gather all simulated QoIs to a single CSV file for SA or UQ analysis later.
        qois = pd.DataFrame({})
        if qoi_group_name == 'ecg':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'ecg_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = qois.append(pd.DataFrame(qoi, index=[simulation_i]))
            qois.to_csv(qoi_save_dir + 'ecg_qois.csv')
        elif qoi_group_name == 'pv':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'pv_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = qois.append(pd.DataFrame(qoi, index=[simulation_i]))
            qois.to_csv(qoi_save_dir + 'pv_qois.csv')
        elif qoi_group_name == 'deformation':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'deformation_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = qois.append(pd.DataFrame(qoi, index=[simulation_i]))
            qois.to_csv(qoi_save_dir + 'deformation_qois.csv')
        elif qoi_group_name == 'fibrework':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'fibrework_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = qois.append(pd.DataFrame(qoi, index=[simulation_i]))
            qois.to_csv(qoi_save_dir + 'fibrework_qois.csv')
        self.qois_db = qois
        return postp_objects



    def sort_simulations(self):
        with open(self.simulation_dir + '/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        self.finished_simulation_dirs = []
        for simulation_i in range(len(all_simulation_dirs)):
            if os.path.exists(all_simulation_dirs[simulation_i].split()[0] + '/results_csv/timeset_1.csv'):
                self.finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])
            # if os.path.exists(all_simulation_dirs[simulation_i].split()[0] + '/heart.post.alyafil'):
            #     with open(all_simulation_dirs[simulation_i].split()[0] + '/heart.post.alyafil', 'r') as f:
            #         if len(f.readlines()) >= 1000:
            #             self.finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])


    def visualise_uq(self, beat, parameter_name, ecg_post=None, pv_post=None, deformation_post=None, fibre_work_post=None):
        with open(self.simulation_dir+'/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        finished_simulation_dirs = []
        for simulation_i in range(len(all_simulation_dirs)):
            if os.path.exists(all_simulation_dirs[simulation_i].split()[0]+'/heart.post.alyafil'):
                with open(all_simulation_dirs[simulation_i].split()[0]+'/heart.post.alyafil', 'r') as f:
                    if len(f.readlines()) >= 1000:
                        finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])
        if len(finished_simulation_dirs) == 3:
            print('Visualising UQ for ', parameter_name)
            alya_output_dir = finished_simulation_dirs[0]
            json_file = self.simulation_dir + 'uq_0.json'
            # post = PostProcessing(alya=alya, simulation_json_file=json_file,
            #                       alya_output_dir=alya_output_dir, verbose=self.verbose)
            # post.read_ecg_pv()
            pvt_lower = pv_post[0].pvs['ts'][beat-1]
            vl_lower = pv_post[0].pvs['vls'][beat-1]
            pl_lower = pv_post[0].pvs['pls'][beat-1]/10000
            ecgt_lower = ecg_post[0].ecgs['ts'][beat-1]
            V3_lower = ecg_post[0].ecgs['V3s'][beat-1]

            pvt_baseline = pv_post[1].pvs['ts'][beat - 1]
            vl_baseline = pv_post[1].pvs['vls'][beat - 1]
            pl_baseline = pv_post[1].pvs['pls'][beat - 1]/10000
            ecgt_baseline = ecg_post[1].ecgs['ts'][beat - 1]
            V3_baseline = ecg_post[1].ecgs['V3s'][beat - 1]

            pvt_upper = pv_post[2].pvs['ts'][beat - 1]
            vl_upper = pv_post[2].pvs['vls'][beat - 1]
            pl_upper = pv_post[2].pvs['pls'][beat - 1]/10000
            ecgt_upper = ecg_post[2].ecgs['ts'][beat - 1]
            V3_upper = ecg_post[2].ecgs['V3s'][beat - 1]

            # Plot PV and ECG UQ
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
            fig.suptitle(parameter_name)
            gs = GridSpec(2,4)
            ax_vol = fig.add_subplot(gs[0,0])
            ax_vol.set_title('LVV(t)')
            ax_vol.plot(pvt_lower, vl_lower, 'k', linestyle='dotted')
            ax_vol.plot(pvt_baseline, vl_baseline, 'k', linestyle='dashed')
            ax_vol.plot(pvt_upper, vl_upper, 'k', linestyle='solid')
            ax_vol.fill(np.append(pvt_lower, pvt_upper[::-1]), np.append(vl_lower, vl_upper[::-1]), alpha=0.3,
                        edgecolor=None, color='k')

            ax_p = fig.add_subplot(gs[0, 1])
            ax_p.set_title('LVP(t)')
            ax_p.plot(pvt_lower, pl_lower, 'k', linestyle='dotted')
            ax_p.plot(pvt_baseline, pl_baseline, 'k', linestyle='dashed')
            ax_p.plot(pvt_upper, pl_upper, 'k', linestyle='solid')
            ax_p.fill(np.append(pvt_lower, pvt_upper[::-1]), np.append(pl_lower, pl_upper[::-1]), alpha=0.3,
                        edgecolor=None, color='k')

            ax_pv = fig.add_subplot(gs[0, 2])
            ax_pv.set_title('LV PV loop')
            ax_pv.plot(vl_lower, pl_lower, 'k', linestyle='dotted')
            ax_pv.plot(vl_baseline, pl_baseline, 'k', linestyle='dashed')
            ax_pv.plot(vl_upper, pl_upper, 'k', linestyle='solid')
            ax_pv.fill(np.append(vl_lower, vl_upper[::-1]), np.append(pl_lower, pl_upper[::-1]), alpha=0.3,
                      edgecolor=None, color='k')

            ax_ecg = fig.add_subplot(gs[0, 3])
            ax_ecg.set_title('ECG V3')
            ax_ecg.plot(ecgt_lower, V3_lower, 'k', linestyle='dotted')
            ax_ecg.plot(ecgt_baseline, V3_baseline, 'k', linestyle='dashed')
            ax_ecg.plot(ecgt_upper, V3_upper, 'k', linestyle='solid')
            ax_ecg.fill(np.append(ecgt_lower, ecgt_upper[::-1]), np.append(V3_lower, V3_upper[::-1]), alpha=0.3,
                      edgecolor=None, color='k')
            if deformation_post:
                deformation_t_lower = deformation_post[0].deformation_transients['deformation_t']
                avpd_lower = deformation_post[0].deformation_transients['avpd']
                wall_thickness_lower = deformation_post[0].deformation_transients['lv_wall_thickness']
                apical_d_lower = deformation_post[0].deformation_transients['apical_displacement']

                deformation_t_baseline = deformation_post[1].deformation_transients['deformation_t']
                avpd_baseline = deformation_post[1].deformation_transients['avpd']
                wall_thickness_baseline = deformation_post[1].deformation_transients['lv_wall_thickness']
                apical_d_baseline = deformation_post[1].deformation_transients['apical_displacement']

                deformation_t_upper = deformation_post[2].deformation_transients['deformation_t']
                avpd_upper = deformation_post[2].deformation_transients['avpd']
                wall_thickness_upper = deformation_post[2].deformation_transients['lv_wall_thickness']
                apical_d_upper = deformation_post[2].deformation_transients['apical_displacement']

                ax_vol = fig.add_subplot(gs[1, 0])
                ax_vol.set_title('Deformations')
                ax_vol.plot(deformation_t_lower, avpd_lower, 'b', deformation_t_lower, wall_thickness_lower, 'g', deformation_t_lower, apical_d_lower, 'm', linestyle='dotted')
                ax_vol.plot(deformation_t_baseline, avpd_baseline, 'b', deformation_t_baseline, wall_thickness_baseline, 'g', deformation_t_baseline, apical_d_baseline, 'm', linestyle='dashed')
                ax_vol.plot(deformation_t_upper, avpd_upper, 'b', label='AVPD', linestyle='solid')
                ax_vol.plot(deformation_t_upper, wall_thickness_upper, 'g', label='Wall thickness', linestyle='solid')
                ax_vol.plot(deformation_t_upper, apical_d_upper, 'm', label='Apical displacement', linestyle='solid')
                ax_vol.fill(np.append(deformation_t_lower, deformation_t_upper[::-1]), np.append(avpd_lower, avpd_upper[::-1]), alpha=0.3,
                            edgecolor=None, color='b')
                ax_vol.fill(np.append(deformation_t_lower, deformation_t_upper[::-1]), np.append(wall_thickness_lower, wall_thickness_upper[::-1]), alpha=0.3,
                            edgecolor=None, color='g')
                ax_vol.fill(np.append(deformation_t_lower, deformation_t_upper[::-1]),
                            np.append(apical_d_lower, apical_d_upper[::-1]), alpha=0.3,
                            edgecolor=None, color='m')
                ax_vol.legend()

            if fibre_work_post:
                fibrework_t_lower = fibre_work_post[0].fibre_work['fibrework_t']
                apical_lambda_lower = fibre_work_post[0].fibre_work['apical_lambda']
                midv_lambda_lower = fibre_work_post[0].fibre_work['midv_lambda']
                basal_lambda_lower = fibre_work_post[0].fibre_work['basal_lambda']

                fibrework_t_baseline = fibre_work_post[1].fibre_work['fibrework_t']
                apical_lambda_baseline = fibre_work_post[1].fibre_work['apical_lambda']
                midv_lambda_baseline = fibre_work_post[1].fibre_work['midv_lambda']
                basal_lambda_baseline = fibre_work_post[1].fibre_work['basal_lambda']

                fibrework_t_upper = fibre_work_post[2].fibre_work['fibrework_t']
                apical_lambda_upper = fibre_work_post[2].fibre_work['apical_lambda']
                midv_lambda_upper = fibre_work_post[2].fibre_work['midv_lambda']
                basal_lambda_upper = fibre_work_post[2].fibre_work['basal_lambda']

                apical_Ta_lower = fibre_work_post[0].fibre_work['apical_Ta']
                midv_Ta_lower = fibre_work_post[0].fibre_work['midv_Ta']
                basal_Ta_lower = fibre_work_post[0].fibre_work['basal_Ta']
                apical_Ta_baseline = fibre_work_post[1].fibre_work['apical_Ta']
                midv_Ta_baseline = fibre_work_post[1].fibre_work['midv_Ta']
                basal_Ta_baseline = fibre_work_post[1].fibre_work['basal_Ta']
                apical_Ta_upper = fibre_work_post[2].fibre_work['apical_Ta']
                midv_Ta_upper = fibre_work_post[2].fibre_work['midv_Ta']
                basal_Ta_upper = fibre_work_post[2].fibre_work['basal_Ta']

                ax_Ta = fig.add_subplot(gs[1, 1])
                ax_Ta.set_title('Ta(t)')
                ax_Ta.plot(fibrework_t_lower, apical_Ta_lower, 'C0',
                           fibrework_t_lower, midv_Ta_lower, 'C1',
                           fibrework_t_lower, basal_Ta_lower, 'C2', linestyle='dotted')
                ax_Ta.plot(fibrework_t_baseline, apical_Ta_baseline, 'C0',
                           fibrework_t_baseline, midv_Ta_baseline, 'C1',
                           fibrework_t_baseline, basal_Ta_baseline, 'C2', linestyle='dashed')
                ax_Ta.plot(fibrework_t_upper, apical_Ta_upper, 'C0', label='Apical', linestyle='solid')
                ax_Ta.plot(fibrework_t_upper, midv_Ta_upper, 'C1', label='Midventricular', linestyle='solid')
                ax_Ta.plot(fibrework_t_upper, basal_Ta_upper, 'C2', label='Basal', linestyle='solid')
                ax_Ta.fill(np.append(fibrework_t_lower, fibrework_t_upper[::-1]),
                               np.append(apical_Ta_lower, apical_Ta_upper[::-1]), alpha=0.3,
                               edgecolor=None, color='C0')
                ax_Ta.fill(np.append(fibrework_t_lower, fibrework_t_upper[::-1]),
                               np.append(midv_Ta_lower, midv_Ta_upper[::-1]), alpha=0.3,
                               edgecolor=None, color='C1')
                ax_Ta.fill(np.append(fibrework_t_lower, fibrework_t_upper[::-1]),
                               np.append(basal_Ta_lower, basal_Ta_upper[::-1]), alpha=0.3,
                               edgecolor=None, color='C2')
                ax_Ta.legend()

                ax_lambda = fig.add_subplot(gs[1, 2])
                ax_lambda.set_title('Lambda(t)')
                ax_lambda.plot(fibrework_t_lower, apical_lambda_lower, 'C0',
                               fibrework_t_lower, midv_lambda_lower, 'C1',
                            fibrework_t_lower, basal_lambda_lower, 'C2', linestyle='dotted')
                ax_lambda.plot(fibrework_t_baseline, apical_lambda_baseline, 'C0',
                               fibrework_t_baseline, midv_lambda_baseline, 'C1',
                               fibrework_t_baseline, basal_lambda_baseline, 'C2', linestyle='dashed')
                ax_lambda.plot(fibrework_t_upper, apical_lambda_upper, 'C0', label='Apical', linestyle='solid')
                ax_lambda.plot(fibrework_t_upper, midv_lambda_upper, 'C1', label='Midventricular', linestyle='solid')
                ax_lambda.plot(fibrework_t_upper, basal_lambda_upper, 'C2', label='Basal', linestyle='solid')
                ax_lambda.fill(np.append(fibrework_t_lower, fibrework_t_upper[::-1]),
                            np.append(apical_lambda_lower, apical_lambda_upper[::-1]), alpha=0.3,
                            edgecolor=None, color='C0')
                ax_lambda.fill(np.append(fibrework_t_lower, fibrework_t_upper[::-1]),
                            np.append(midv_lambda_lower, midv_lambda_upper[::-1]), alpha=0.3,
                            edgecolor=None, color='C1')
                ax_lambda.fill(np.append(fibrework_t_lower, fibrework_t_upper[::-1]),
                            np.append(basal_lambda_lower, basal_lambda_upper[::-1]), alpha=0.3,
                            edgecolor=None, color='C2')
                ax_lambda.legend()

                ax_fibrework = fig.add_subplot(gs[1, 3])
                ax_fibrework.set_title('Fibre work Ta(lambda)')
                ax_fibrework.plot(apical_lambda_lower, apical_Ta_lower, 'C0',
                                  midv_lambda_lower, midv_Ta_lower, 'C1',
                                  basal_lambda_lower, basal_Ta_lower, 'C2', linestyle='dotted')
                ax_fibrework.plot(apical_lambda_baseline, apical_Ta_baseline, 'C0',
                                  midv_lambda_baseline, midv_Ta_baseline, 'C1',
                                  basal_lambda_baseline, basal_Ta_baseline, 'C2', linestyle='dashed')
                ax_fibrework.plot(apical_lambda_upper, apical_Ta_upper, 'C0',
                                  midv_lambda_upper, midv_Ta_upper, 'C1',
                                  basal_lambda_upper, basal_Ta_upper, 'C2', linestyle='solid')
                ax_fibrework.fill(np.append(apical_lambda_lower, apical_lambda_upper[::-1]),
                            np.append(apical_Ta_lower, apical_Ta_upper[::-1]), alpha=0.3,
                            edgecolor=None, color='C0')
                ax_fibrework.fill(np.append(midv_lambda_lower, midv_lambda_upper[::-1]),
                                  np.append(midv_Ta_lower, midv_Ta_upper[::-1]), alpha=0.3,
                                  edgecolor=None, color='C1')
                ax_fibrework.fill(np.append(basal_lambda_lower, basal_lambda_upper[::-1]),
                                  np.append(basal_Ta_lower, basal_Ta_upper[::-1]), alpha=0.3,
                                  edgecolor=None, color='C2')

            plt.savefig(self.simulation_dir + '/uq_plots_'+parameter_name+'.png')
            plt.show()



    def analyse(self, filename, qois):
        self.qois_db = pd.read_csv(filename, index_col=False)
        names = self.parameter_names
        # selected_qois = self.qois_db.columns.values.tolist() # All QoIs
        selected_qois = qois
        # selected_qois = ['t_pe_mean', 't_peak_mean', 'qt_dur_mean', 't_polarity_mean', 'EDVL', 'LVEF', 'PmaxL', 'SVL']
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
                elif '_lv' in names[param_i]:
                    X[simulation_i,param_i] = dict[names[param_i].split('_lv')[0]][0]
                elif '_rv' in names[param_i]:
                    X[simulation_i, param_i] = dict[names[param_i].split('_rv')[0]][1]
                elif '_myocardium' in names[param_i]:
                    X[simulation_i, param_i] = dict[names[param_i].split('_myocardium')[0]][0]
                elif '_valveplug' in names[param_i]:
                    X[simulation_i, param_i] = dict[names[param_i].split('_valveplug')[0]][1]
                else:
                    X[simulation_i,param_i] = dict[names[param_i]]

        # Scatter plots with correlation coefficients
        if len(Y.shape) == 1:
            num_qois = 1
        else:
            num_qois = Y.shape[1]

        ################################################################################################################
        fig = plt.figure(tight_layout=True, figsize=(18, 10))
        fig.suptitle('N=' + str(Y.shape[0]))
        gs = GridSpec(num_qois, X.shape[1])
        for qoi_i in range(num_qois):
            for param_j in range(X.shape[1]):
                ax = fig.add_subplot(gs[qoi_i, param_j])
                x = X[:,param_j]
                if num_qois == 1:
                    y = Y
                else:
                    y = Y[:,qoi_i]
                sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s':1})
                ax.text(x=np.amin(x), y=np.amax(y), va='top', ha='left',
                                          s='p=%.2f' % (np.corrcoef(x,y)[0,1]))
                if qoi_i == num_qois-1:
                    ax.set_xlabel(names[param_j])
                if param_j == 0:
                    ax.set_ylabel(qoi_names[qoi_i])
        plt.savefig('scatter_qois_vs_parameters.png')
        plt.show()
        quit()
        ###############################################################################################################
        fig = plt.figure(tight_layout=True, figsize=(18, 10))
        fig.suptitle('N=' + str(Y.shape[0]))
        gs = GridSpec(num_qois, num_qois)
        for qoi_i in range(num_qois):
            for qoi_j in range(num_qois):
                if qoi_i >= qoi_j:
                    ax = fig.add_subplot(gs[qoi_i, qoi_j])
                    if num_qois == 1:
                        y = Y
                    else:
                        y = Y[:, qoi_i]
                    x = Y[:, qoi_j]
                    sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s': 1})
                    ax.text(x=np.amin(x), y=np.amax(y), va='top', ha='left',
                            s='p=%.2f' % (np.corrcoef(x, y)[0, 1]))
                if qoi_i == num_qois - 1:
                    ax.set_xlabel(qoi_names[qoi_j])
                if qoi_j == 0:
                    ax.set_ylabel(qoi_names[qoi_i])
        plt.savefig('scatter_qois_vs_qois.png')
        plt.show()
        ####################################################################################################################
        # Tornado plot of sensitivity indices https://seaborn.pydata.org/examples/part_whole_bars.html
        sp.set_results(Y)
        Si = sp.analyze_sobol(print_to_console=False, calc_second_order=False)
        fig = plt.figure(tight_layout=True, figsize=(15,6))
        fig.suptitle('N=' + str(Y.shape[0]))
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

        
    # def interpolate_displacement_to_lower_resolution(self, simulation_dir, output_dir):
    #     with open(simulation_dir + '/all_simulation_dirs.txt', 'r') as f:
    #         all_simulation_dirs = f.readlines()
    #     geometric_desired_data_dir = '/p/project/icei-prace-2022-0003/wang1/Personalisation_projects/meta_data/geometric_data/rodero_05/rodero_05_coarse/'
    #     desired_file_prefix = anatomy_subject_name + '_' + desired_resolution
    #
    #     for simulation_i in range(len(all_simulation_dirs)):

