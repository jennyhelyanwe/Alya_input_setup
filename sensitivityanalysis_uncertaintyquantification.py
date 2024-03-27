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
from healthy_qoi_ranges import HealthyBiomarkerRanges
class SAUQ:
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
        self.healthy_ranges = HealthyBiomarkerRanges().healthy_ranges

    def setup(self, upper_bounds, lower_bounds):
        self.generate_parameter_set( upper_bounds, lower_bounds)
        print('Number of samples/simulations to set up: ', str(self.parameter_set.shape[0]))
        print('Number of parameters included: ', str(self.parameter_set.shape[1]))
        self.generate_alya_simulation_json(self.baseline_json_file)

    def setup_fibre_sa(self, fields, fibre_filenames, sheet_filenames, normal_filenames, map_filename, baseline_json_file):
        print('Set up Alya simulations with difference fibre fields')
        baseline_simulation_dict = json.load(open(baseline_json_file, 'r'))
        map = fields.map_doste_nodes_to_rodero_nodes(fibre_vtk_filename=fibre_filenames[0], map_filename=map_filename)
        for fibre_i in range(len(fibre_filenames)):
            print('Reading in fibres: ', fibre_filenames[fibre_i])
            fields.read_doste_fibre_fields_vtk(
                fibre_vtk_filename=fibre_filenames[fibre_i],
                sheet_vtk_filename=sheet_filenames[fibre_i],
                normal_vtk_filename=normal_filenames[fibre_i], save=False, map=map)
            print('Node fields fibre: ', fields.node_fields.dict['fibre'])
            with open(self.simulation_dir + self.name + '_' + str(fibre_i) + '.json', 'w') as f:
                json.dump(baseline_simulation_dict, f)
            self.alya_format.simulation_dir = self.simulation_dir
            self.alya_format.node_fields = fields.node_fields
            print('Writing out Alya simulation files')
            self.alya_format.do(simulation_json_file=self.simulation_dir + self.name + '_' + str(fibre_i) + '.json',
                                SA_flag=False, baseline_dir=self.baseline_dir)
            self.all_simulation_dirs.append(self.alya_format.output_dir)
        with open(self.simulation_dir + '/all_simulation_dirs.txt', 'w') as f:
            for i in range(len(self.all_simulation_dirs)):
                f.write(self.all_simulation_dirs[i] + '\n')

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
        elif self.sampling_method == 'uniform':
            self.parameter_set = np.linspace(lower_bounds, upper_bounds, self.n)
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
            if qoi_group_name == 'ecg':
                post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                      alya_output_dir=alya_output_dir, protocol='raw',
                                      verbose=self.verbose)
                post.evaluate_ecg_biomarkers(beat=beat)
                post.save_qoi(filename=qoi_save_dir + 'ecg_qoi_' + str(simulation_i) + '.csv')
                postp_objects.append(post)
            elif qoi_group_name == 'pv':
                post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                      alya_output_dir=alya_output_dir, protocol='raw',
                                      verbose=self.verbose)
                post.evaluate_pv_biomarkers(beat=beat)
                post.save_qoi(filename=qoi_save_dir + 'pv_qoi_' + str(simulation_i) + '.csv')
                postp_objects.append(post)
            elif qoi_group_name == 'deformation':
                post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                      alya_output_dir=alya_output_dir, protocol='postprocess',
                                      verbose=self.verbose)
                post.evaluate_deformation_biomarkers(beat=beat)
                post.save_qoi(filename=qoi_save_dir + 'deformation_qoi_' + str(simulation_i) + '.csv')
                postp_objects.append(post)
            elif qoi_group_name == 'strain':
                post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                      alya_output_dir=alya_output_dir, protocol='postprocess',
                                      verbose=self.verbose)
                post.evaluate_strain_biomarkers(beat=beat)
                post.save_qoi(filename=qoi_save_dir + 'strain_qoi_' + str(simulation_i) + '.csv')
                postp_objects.append(post)
            elif qoi_group_name == 'fibre_work':
                post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                      alya_output_dir=alya_output_dir, protocol='postprocess',
                                      verbose=self.verbose)
                post.evaluate_fibre_work_biomarkers(beat=beat)
                post.save_qoi(filename=qoi_save_dir + 'fibrework_qoi_' + str(simulation_i) + '.csv')
                postp_objects.append(post)
            elif qoi_group_name == 'cube_deformation_ta':
                post = PostProcessing(alya=alya, simulation_json_file=json_file,
                                      alya_output_dir=alya_output_dir, protocol='postprocess',
                                      verbose=self.verbose)
                post.evaluate_cube_deformation_ta_biomarkers()
                post.save_qoi(filename=qoi_save_dir + 'cube_deformation_ta_qoi_' + str(simulation_i) + '.csv')
                postp_objects.append(post)


        # Gather all simulated QoIs to a single CSV file for SAUQ or UQ analysis later.
        qois = pd.DataFrame({})
        if qoi_group_name == 'ecg':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'ecg_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = pd.concat([qois, pd.DataFrame([qoi])])
            qois.to_csv(qoi_save_dir + 'ecg_qois.csv')
        elif qoi_group_name == 'pv':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'pv_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = pd.concat([qois, pd.DataFrame([qoi])])
            qois.to_csv(qoi_save_dir + 'pv_qois.csv')
        elif qoi_group_name == 'deformation':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'deformation_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = pd.concat([qois, pd.DataFrame([qoi])])
            qois.to_csv(qoi_save_dir + 'deformation_qois.csv')
        elif qoi_group_name == 'strain':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'strain_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = pd.concat([qois, pd.DataFrame([qoi])])
            qois.to_csv(qoi_save_dir + 'strain_qois.csv')
        elif qoi_group_name == 'fibre_work':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'fibrework_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = pd.concat([qois, pd.DataFrame([qoi])])
            qois.to_csv(qoi_save_dir + 'fibrework_qois.csv')
        elif qoi_group_name == 'cube_deformation':
            for simulation_i in range(len(self.finished_simulation_dirs)):
                filename = qoi_save_dir + 'cube_deformation_qoi_' + str(simulation_i) + '.csv'
                qoi = json.load(open(filename, 'r'))
                qois = pd.concat([qois, pd.DataFrame([qoi])])
            qois.to_csv(qoi_save_dir + 'cube_deformation_qois.csv')
        self.qois_db = qois
        return postp_objects


    def visualise_ed_es_pvr_biomarkers(self, beat, pv_post):
        fig = plt.figure(tight_layout=True, figsize=(18, 10))
        gs = GridSpec(1, 2)
        ax_pv = fig.add_subplot(gs[0, 0])
        ax_edpvr = fig.add_subplot(gs[0, 1])
        lvedp = []
        lvedv = []
        lvesp = []
        lvesv = []
        V0 = pv_post[0].pvs['vls'][beat-1][0] # Volume at EDP = 0 mmHg
        V30 = 0
        # Find V30 - volume when pressure = 30 mmHg = 4 kPa
        for simulation_i in range(len(pv_post)):
            p30_idx = np.where(abs(pv_post[simulation_i].pvs['pls'][beat-1] - 4)<0.01)[0]
            if p30_idx.size:
                V30 = pv_post[simulation_i].pvs['vls'][beat-1][p30_idx[0]]

        for simulation_i in range(len(pv_post)):
            ax_pv.plot(pv_post[simulation_i].pvs['vls'][beat - 1], pv_post[simulation_i].pvs['pls'][beat - 1] / 10000,
                       color='C0', label='LV')
            ax_pv.plot(pv_post[simulation_i].pvs['vrs'][beat - 1], pv_post[simulation_i].pvs['prs'][beat - 1] / 10000,
                       color='C1', label='RV')
            # ax_vt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['vls'][beat - 1],
            #            color='C0', label='LV')
            # ax_vt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['vrs'][beat - 1],
            #            color='C1', label='RV')
            # ax_pt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['pls'][beat - 1] / 10000*7.5,
            #            color='C0', label='LV')
            # ax_pt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['prs'][beat - 1] / 10000*7.5,
            #            color='C1', label='RV')
            lved_idx = np.where(pv_post[simulation_i].pvs['phasels'][beat - 1] == 0)[0][-1]
            lves_idx = np.where(pv_post[simulation_i].pvs['phasels'][beat - 1] == 2)[0][-1]
            lvedv.append(pv_post[simulation_i].pvs['vls'][beat - 1][lved_idx])
            lvedp.append(pv_post[simulation_i].pvs['pls'][beat - 1][lved_idx]/10000)
            lvesv.append(pv_post[simulation_i].pvs['vls'][beat - 1][lves_idx])
            lvesp.append(pv_post[simulation_i].pvs['pls'][beat - 1][lves_idx]/10000)
        ax_pv.set_xlabel('Volume (mL)')
        ax_pv.set_ylabel('Pressure (kPa)')
        # ax_vt.set_xlabel('Time (s)')
        # ax_vt.set_ylabel('Volume (mL)')
        # ax_pt.set_xlabel('Time (s)')
        # ax_pt.set_ylabel('Pressure (mmHg)')
        ed_idx = np.argsort(lvedv)
        es_idx = np.argsort(lvesv)
        lvedv = np.array(lvedv)[ed_idx]
        lvesv = np.array(lvesv)[es_idx]
        lvedp = np.array(lvedp)[ed_idx]
        lvesp = np.array(lvesp)[es_idx]
        es_a, es_b = np.polyfit(lvesv, lvesp, 1)
        # ax_pv.plot(lvedv, lvedp, color='k', linestyle='--')
        # ax_pv.axline((lvesv[0], lvesv[0]), slope=es_a, color='k', linestyle='--')
        # ax_pv.scatter(lvesv, lvesp)
        # # Plot Klotz curve for EDPVR
        # ax_pv.plot(lvedv, self.healthy_ranges['edpvr_a_klotz'] * lvedv ** self.healthy_ranges['edpvr_b_klotz'], color='g', linestyle='--')
        # # Plot healthy ranges for EDVPR and ESPVR
        # es_intercept_0 = lvesp[0] - self.healthy_ranges['espvr'][0] * lvesv[0]
        # es_intercept_1 = lvesp[0] - self.healthy_ranges['espvr'][1] * lvesv[0]
        # ax_pv.fill_between(lvesv, self.healthy_ranges['espvr'][0] * lvesv + es_intercept_0,
        #                    self.healthy_ranges['espvr'][1] * lvesv + es_intercept_1,
        #                    alpha=0.3, facecolor='green')
        plt.show()


    def sort_simulations(self, tag='postprocess'):
        with open(self.simulation_dir + '/all_simulation_dirs.txt', 'r') as f:
            all_simulation_dirs = f.readlines()
        self.finished_simulation_dirs = []
        if tag=='postprocess':
            for simulation_i in range(len(all_simulation_dirs)):
                if os.path.exists(all_simulation_dirs[simulation_i].split()[0] + '/heart-cardiac-cycle.sld.res'):
                    if os.path.exists(all_simulation_dirs[simulation_i].split()[0] + '/results_csv/timeset_1.csv'):
                        self.finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])
        elif tag== 'raw':
            for simulation_i in range(len(all_simulation_dirs)):
                if os.path.exists(all_simulation_dirs[simulation_i].split()[0]+'/heart-cardiac-cycle.sld.res'):
                    self.finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])
        elif tag == 'cube_postprocess':
            for simulation_i in range(len(all_simulation_dirs)):
                if os.path.exists(all_simulation_dirs[simulation_i].split()[0] + '/results_csv/timeset_1.csv'):
                    self.finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])
            # if os.path.exists(all_simulation_dirs[simulation_i].split()[0] + '/heart.post.alyafil'):
            #     with open(all_simulation_dirs[simulation_i].split()[0] + '/heart.post.alyafil', 'r') as f:
            #         if len(f.readlines()) >= 1000:
            #             self.finished_simulation_dirs.append(all_simulation_dirs[simulation_i].split()[0])


    def visualise_sa(self, beat, ecg_post=None, pv_post=None, deformation_post=None,
                     fibre_work_post=None, strain_post=None, cube_deformation_ta_post=None, labels=None, save_filename=None):
        if pv_post:
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
            gs = GridSpec(2, 2)
            ax_pv = fig.add_subplot(gs[0:2, 0])
            ax_vt = fig.add_subplot(gs[0, 1])
            ax_pt = fig.add_subplot(gs[1, 1])
            for simulation_i in range(len(pv_post)):
                if labels:
                    ax_pv.plot(pv_post[simulation_i].pvs['vls'][beat - 1],
                               pv_post[simulation_i].pvs['pls'][beat - 1] / 10000, color='C0', label=labels[simulation_i])
                else:
                    ax_pv.plot(pv_post[simulation_i].pvs['vls'][beat-1], pv_post[simulation_i].pvs['pls'][beat-1]/10000, color='C0', label='LV')
                ax_pv.plot(pv_post[simulation_i].pvs['vrs'][beat - 1], pv_post[simulation_i].pvs['prs'][beat - 1]/10000,
                           color='C1')
                if labels:
                    ax_vt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['vls'][beat - 1],
                               color='C0', label=labels[simulation_i])
                else:
                    ax_vt.plot(pv_post[simulation_i].pvs['ts'][beat-1], pv_post[simulation_i].pvs['vls'][beat-1], color='C0', label='LV')
                ax_vt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['vrs'][beat - 1],
                           color='C1')

                if labels:
                    ax_pt.plot(pv_post[simulation_i].pvs['ts'][beat - 1],
                               pv_post[simulation_i].pvs['pls'][beat - 1] / 10000,
                               color='C0', label=labels[simulation_i])
                else:
                    ax_pt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['pls'][beat - 1]/10000,
                               color='C0', label='LV')
                ax_pt.plot(pv_post[simulation_i].pvs['ts'][beat - 1], pv_post[simulation_i].pvs['prs'][beat - 1]/10000,
                           color='C1')
            ax_pv.set_xlabel('Volume (mL)')
            ax_pv.set_ylabel('Pressure (kPa)')
            ax_vt.set_xlabel('Time (s)')
            ax_vt.set_ylabel('Volume (mL)')
            ax_pt.set_xlabel('Time (s)')
            ax_pt.set_ylabel('Pressure (kPa)')
            if labels:
                ax_pv.legend()
                ax_vt.legend()
                ax_pt.legend()
        if ecg_post:
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
            gs = GridSpec(1, 6)
            ax_V1 = fig.add_subplot(gs[0, 0])
            ax_V2 = fig.add_subplot(gs[0, 1])
            ax_V3 = fig.add_subplot(gs[0, 2])
            ax_V4 = fig.add_subplot(gs[0, 3])
            ax_V5 = fig.add_subplot(gs[0, 4])
            ax_V6 = fig.add_subplot(gs[0, 5])
            for simulation_i in range(len(ecg_post)):
                if labels:
                    ax_V1.plot(ecg_post[simulation_i].ecgs['ts'][beat-1],
                               ecg_post[simulation_i].ecgs['V1s'][beat-1]/ecg_post[simulation_i].ecgs['max_all_leads'],
                               label=labels[simulation_i])
                else:
                    ax_V1.plot(ecg_post[simulation_i].ecgs['ts'][beat - 1],
                               ecg_post[simulation_i].ecgs['V1s'][beat - 1] / ecg_post[simulation_i].ecgs[
                                   'max_all_leads'])
                ax_V2.plot(ecg_post[simulation_i].ecgs['ts'][beat - 1],
                           ecg_post[simulation_i].ecgs['V2s'][beat - 1] / ecg_post[simulation_i].ecgs['max_all_leads'])
                ax_V3.plot(ecg_post[simulation_i].ecgs['ts'][beat - 1],
                           ecg_post[simulation_i].ecgs['V3s'][beat - 1] / ecg_post[simulation_i].ecgs['max_all_leads'])
                ax_V4.plot(ecg_post[simulation_i].ecgs['ts'][beat - 1],
                           ecg_post[simulation_i].ecgs['V4s'][beat - 1] / ecg_post[simulation_i].ecgs['max_all_leads'])
                ax_V5.plot(ecg_post[simulation_i].ecgs['ts'][beat - 1],
                           ecg_post[simulation_i].ecgs['V5s'][beat - 1] / ecg_post[simulation_i].ecgs['max_all_leads'])
                ax_V6.plot(ecg_post[simulation_i].ecgs['ts'][beat - 1],
                           ecg_post[simulation_i].ecgs['V6s'][beat - 1] / ecg_post[simulation_i].ecgs['max_all_leads'])
            ax_V1.set_xlabel('Time (s)')
            ax_V1.set_title('V1')
            ax_V1.set_ylim([-1, 1])
            ax_V2.set_xlabel('Time (s)')
            ax_V2.set_title('V2')
            ax_V2.set_ylim([-1, 1])
            ax_V3.set_xlabel('Time (s)')
            ax_V3.set_title('V3')
            ax_V3.set_ylim([-1, 1])
            ax_V4.set_xlabel('Time (s)')
            ax_V4.set_title('V4')
            ax_V4.set_ylim([-1, 1])
            ax_V5.set_xlabel('Time (s)')
            ax_V5.set_title('V5')
            ax_V5.set_ylim([-1, 1])
            ax_V6.set_xlabel('Time (s)')
            ax_V6.set_title('V6')
            ax_V6.set_ylim([-1, 1])
            if labels:
                ax_V1.legend()
        if deformation_post:
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
            gs = GridSpec(1, 3)
            ax_avpd = fig.add_subplot(gs[0,0])
            ax_apex = fig.add_subplot(gs[0,1])
            ax_wall = fig.add_subplot(gs[0,2])
            for simulation_i in range(len(deformation_post)):
                if labels:
                    ax_avpd.plot(deformation_post[simulation_i].deformation_transients['deformation_t'],
                                 deformation_post[simulation_i].deformation_transients['avpd'],
                                 label=labels[simulation_i])
                    ax_apex.plot(deformation_post[simulation_i].deformation_transients['deformation_t'],
                                 deformation_post[simulation_i].deformation_transients['apical_displacement'],
                                 label=labels[simulation_i])
                    ax_wall.plot(deformation_post[simulation_i].deformation_transients['deformation_t'],
                                 deformation_post[simulation_i].deformation_transients['lv_wall_thickness'],
                                 label=labels[simulation_i])
                else:
                    ax_avpd.plot(deformation_post[simulation_i].deformation_transients['deformation_t'],
                                 deformation_post[simulation_i].deformation_transients['avpd'])
                    ax_apex.plot(deformation_post[simulation_i].deformation_transients['deformation_t'],
                                 deformation_post[simulation_i].deformation_transients['apical_displacement'])
                    ax_apex.plot(deformation_post[simulation_i].deformation_transients['deformation_t'],
                                 deformation_post[simulation_i].deformation_transients['lv_wall_thickness'])
            ax_avpd.set_title('AVPD')
            ax_avpd.set_xlabel('Time (s)')
            ax_avpd.set_ylabel('AVPD (cm)')
            ax_apex.set_title('Apical displacement')
            ax_apex.set_xlabel('Time (s)')
            ax_apex.set_ylabel('(cm)')
            ax_wall.set_title('Wall thickness')
            ax_wall.set_xlabel('Time (s)')
            ax_wall.set_ylabel('Wall thickness (cm)')
            if labels:
                ax_avpd.legend()
                ax_apex.legend()
                ax_wall.legend()
        if fibre_work_post:
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
            gs = GridSpec(1, 2)
            ax_lambda = fig.add_subplot(gs[0,0])
            ax_ta = fig.add_subplot(gs[0,1])
            for simulation_i in range(len(fibre_work_post)):
                if labels:
                    ax_lambda.plot(fibre_work_post[simulation_i].fibre_work_transients['fibrework_t'],
                                 fibre_work_post[simulation_i].fibre_work_transients['mean_lambda'],
                                 label=labels[simulation_i])
                    ax_ta.plot(fibre_work_post[simulation_i].fibre_work_transients['fibrework_t'],
                                 fibre_work_post[simulation_i].fibre_work_transients['mean_Ta'],
                                 label=labels[simulation_i])
                else:
                    ax_lambda.plot(deformation_post[simulation_i].deformation_transients['fibrework_t'],
                                 deformation_post[simulation_i].deformation_transients['mean_lambda'])
                    ax_ta.plot(deformation_post[simulation_i].deformation_transients['fibrework_t'],
                                 deformation_post[simulation_i].deformation_transients['mean_Ta'])
            ax_lambda.set_title('Fibre stretch ratio')
            ax_lambda.set_xlabel('Time (s)')
            ax_lambda.set_ylabel('Lambda')
            ax_ta.set_title('Active tension')
            ax_ta.set_xlabel('Time (s)')
            ax_ta.set_ylabel('Ta (kPa)')
            if labels:
                ax_lambda.legend()
                ax_ta.legend()
        if strain_post:
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
            gs = GridSpec(1, 3)
            ax_Err = fig.add_subplot(gs[0, 0])
            ax_Ecc = fig.add_subplot(gs[0, 1])
            ax_Ell = fig.add_subplot(gs[0, 2])
            for simulation_i in range(len(strain_post)):
                if labels:
                    ax_Ecc.plot(strain_post[simulation_i].strain_transients['strain_t'],
                                strain_post[simulation_i].strain_transients['mean_mid_E_cc'], label=labels[simulation_i])
                    ax_Err.plot(strain_post[simulation_i].strain_transients['strain_t'],
                                strain_post[simulation_i].strain_transients['mean_mid_E_rr'], label=labels[simulation_i])
                    ax_Ell.plot(strain_post[simulation_i].strain_transients['strain_t'],
                                strain_post[simulation_i].strain_transients['mean_four_chamber_E_ll'], label=labels[simulation_i])
                else:
                    ax_Ecc.plot(strain_post[simulation_i].strain_transients['strain_t'], strain_post[simulation_i].strain_transients['mean_mid_E_cc'])
                    ax_Err.plot(strain_post[simulation_i].strain_transients['strain_t'], strain_post[simulation_i].strain_transients['mean_mid_E_rr'])
                    ax_Ell.plot(strain_post[simulation_i].strain_transients['strain_t'], strain_post[simulation_i].strain_transients['mean_four_chamber_E_ll'])
            ax_Err.set_title('Mid-vent Err')
            ax_Err.set_xlabel('Time (s)')
            ax_Ecc.set_title('Mid-vent Ecc')
            ax_Ecc.set_xlabel('Time (s)')
            ax_Ell.set_title('Four chamber Ell')
            ax_Ell.set_xlabel('Time (s)')
            if labels:
                ax_Ecc.legend()
                ax_Err.legend()
                ax_Ell.legend()
        if cube_deformation_ta_post:
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
            gs = GridSpec(1, 2)
            ax_displ = fig.add_subplot(gs[0, 0])
            ax_ta = fig.add_subplot(gs[0, 1])
            for simulation_i in range(len(cube_deformation_ta_post)):
                if labels:
                    ax_displ.plot(cube_deformation_ta_post[simulation_i].deformation_transients['deformation_t'],
                                 cube_deformation_ta_post[simulation_i].deformation_transients['mean_displ_x'], label=labels[simulation_i])
                    ax_ta.plot(cube_deformation_ta_post[simulation_i].deformation_transients['deformation_t'],
                               cube_deformation_ta_post[simulation_i].deformation_transients['mean_ta'],
                               label=labels[simulation_i])
                else:
                    ax_displ.plot(cube_deformation_ta_post[simulation_i].deformation_transients['deformation_t'],
                                  cube_deformation_ta_post[simulation_i].deformation_transients['mean_displ_x'])
                    ax_ta.plot(cube_deformation_ta_post[simulation_i].deformation_transients['deformation_t'],
                                  cube_deformation_ta_post[simulation_i].deformation_transients['mean_ta'])

            ax_displ.set_title('X displacement')
            ax_displ.set_xlabel('Time (s)')
            ax_displ.set_ylabel('(cm)')
            ax_ta.set_title('Ta')
            ax_ta.set_xlabel('Time (s)')
            ax_ta.set_ylabel('(kPa)')
            if labels:
                ax_displ.legend()
                ax_ta.legend()
        print('Showing figure...')
        if save_filename:
            plt.savefig(save_filename)
        # plt.show()


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
            #                       alya_output_dir=aly a_output_dir, verbose=self.verbose)
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

    def analyse(self, filename, qois, show_healthy_ranges=True, save_filename=None):
        self.qois_db = pd.read_csv(filename, index_col=False)
        names = self.parameter_names
        # selected_qois = self.qois_db.columns.values.tolist() # All QoIs
        selected_qois = qois
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
        if X.shape[1] < 2:
            fig = plt.figure(tight_layout=True, figsize=(6,10))
        else:
            fig = plt.figure(tight_layout=True, figsize=(18, 10))
        fig.suptitle('N=' + str(Y.shape[0]))
        gs = GridSpec(num_qois, X.shape[1])
        corrs = np.zeros((num_qois, X.shape[1]))
        ranges = np.zeros((num_qois, X.shape[1]))
        for qoi_i in range(num_qois):
            for param_j in range(X.shape[1]):
                ax = fig.add_subplot(gs[qoi_i, param_j])
                x = X[:,param_j]
                if num_qois == 1:
                    y = Y
                else:
                    y = Y[:,qoi_i]
                x[~np.isfinite(x)] = 0
                y[~np.isfinite(y)] = 0
                sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s':1})
                ax.text(x=np.nanmin(x), y=np.nanmax(y), va='top', ha='left',
                                          s='p=%.2f' % (np.corrcoef(x,y)[0,1]) + 'range=%.2f' % (np.amax(y[np.nonzero(y)]) - np.amin(y[np.nonzero(y)])))
                corrs[qoi_i, param_j] = np.corrcoef(x, y)[0, 1]
                ranges[qoi_i, param_j] = np.amax(y) - np.amin(y)
                if qoi_i == num_qois-1:
                    ax.set_xlabel(names[param_j])
                if param_j == 0:
                    ax.set_ylabel(qoi_names[qoi_i])
                if show_healthy_ranges:
                    if qoi_names[qoi_i] == 'LVEF':
                        ax.axhspan(self.healthy_ranges['LVEF'][0], self.healthy_ranges['LVEF'][1], alpha=0.3, facecolor='green')
                    elif qoi_names[qoi_i] == 'EDVL':
                        ax.axhspan(self.healthy_ranges['LVEDV'][0], self.healthy_ranges['LVEDV'][1], alpha=0.3, facecolor='green')
                    elif qoi_names[qoi_i] == 'ESVL':
                        ax.axhspan(self.healthy_ranges['LVESV'][0], self.healthy_ranges['LVESV'][1], alpha=0.3, facecolor='green')
                    elif qoi_names[qoi_i] == 'PmaxL':
                        ax.axhspan(self.healthy_ranges['LVESP'][0], self.healthy_ranges['LVESP'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 'es_ed_avpd':
                        ax.axhspan(self.healthy_ranges['AVPD'][0], self.healthy_ranges['AVPD'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 'es_ed_apical_displacement':
                        ax.axhspan(self.healthy_ranges['apical_displacement'][0],
                                   self.healthy_ranges['apical_displacement'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 'qt_dur_mean':
                        ax.axhspan(self.healthy_ranges['QTc'][0],
                                   self.healthy_ranges['QTc'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 'qrs_dur_mean':
                        ax.axhspan(self.healthy_ranges['QRS_duration'][0],
                                   self.healthy_ranges['QRS_duration'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 't_pe_mean':
                        ax.axhspan(self.healthy_ranges['Tpe'][0],
                                   self.healthy_ranges['Tpe'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 'dvdt_ejection':
                        ax.axhspan(self.healthy_ranges['dvdt_ejection'][0],
                                   self.healthy_ranges['dvdt_ejection'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 'dvdt_filling':
                        ax.axhspan(self.healthy_ranges['dvdt_filling'][0],
                                   self.healthy_ranges['dvdt_filling'][1], alpha=0.3,
                                   facecolor='green')
                    elif qoi_names[qoi_i] == 'dpdt_max':
                        ax.axhspan(self.healthy_ranges['dpdt_max'][0],
                                   self.healthy_ranges['dpdt_max'][1], alpha=0.3,
                                   facecolor='green')
        if save_filename:
            plt.savefig(save_filename)
        plt.show()
        return corrs, ranges
        # ###############################################################################################################
        # fig = plt.figure(tight_layout=True, figsize=(18, 10))
        # fig.suptitle('N=' + str(Y.shape[0]))
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
        # plt.savefig('scatter_qois_vs_qois.png')
        # plt.show()
        # ####################################################################################################################
        # # Tornado plot of sensitivity indices https://seaborn.pydata.org/examples/part_whole_bars.html
        # sp.set_results(Y)
        # Si = sp.analyze_sobol(print_to_console=False, calc_second_order=False)
        # fig = plt.figure(tight_layout=True, figsize=(15,6))
        # fig.suptitle('N=' + str(Y.shape[0]))
        # gs = GridSpec(1, num_qois)
        # data = Si.to_df()
        # sns.set_theme(style="whitegrid")
        # for qoi_i in range(num_qois):
        #     ax = fig.add_subplot(gs[0, qoi_i])
        #     st_s1_data = pd.concat([data[qoi_i][0], data[qoi_i][1]], axis=1)
        #     sorted_data = st_s1_data.reindex(st_s1_data.abs().sort_values('ST', ascending=False).index)
        #     names = []
        #     for row in sorted_data.index:
        #         names.append(row)
        #     sns.set_color_codes("pastel")
        #     sns.barplot(data=sorted_data, x='ST', y=names, label='ST', color='b')
        #     sns.set_color_codes("muted")
        #     sns.barplot(data=sorted_data, x='S1', y=names, label='S1', color='b')
        #
        #     ax.set( ylabel="",
        #            xlabel=qoi_names[qoi_i])
        #     if qoi_i == num_qois-1:
        #         ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', frameon=True)
        # plt.show()
        # plt.savefig('ST_S1_tornado.png')

    def analyse_qoi_vs_qoi(self, filename_1, filename_2, qois_1, qois_2):
        self.qois_db_1 = pd.read_csv(filename_1, index_col=False)
        self.qois_db_2 = pd.read_csv(filename_2, index_col=False)
        names = self.parameter_names
        sp_1 = ProblemSpec({'num_vars': len(names),
                          'names': names,
                          'bounds': [[0.5, 2]] * len(names),
                          'outputs': qois_1
                          })
        sp_2 = ProblemSpec({'num_vars': len(names),
                            'names': names,
                            'bounds': [[0.5, 2]] * len(names),
                            'outputs': qois_2
                            })
        qoi_names_1 = sp_1.get('outputs')
        qoi_names_2 = sp_2.get('outputs')
        Y_1 = self.qois_db_1[qois_1].values
        Y_2 = self.qois_db_2[qois_2].values

        ################################################################################################################
        # Scatter plots with correlation coefficients
        if len(Y_1.shape) == 1:
            num_qois_1 = 1
        else:
            num_qois_1 = Y_1.shape[1]
        if len(Y_2.shape) == 1:
            num_qois_2 = 1
        else:
            num_qois_2 = Y_2.shape[1]
        fig = plt.figure(tight_layout=True, figsize=(18, 10))
        fig.suptitle('N=' + str(Y_1.shape[0]))
        gs = GridSpec(num_qois_2, num_qois_1)
        for qoi_i_1 in range(num_qois_1):
            for qoi_i_2 in range(num_qois_2):
                ax = fig.add_subplot(gs[qoi_i_2, qoi_i_1])
                y = Y_2[:, qoi_i_2]
                x = Y_1[:, qoi_i_1]
                sns.regplot(x=x, y=y, ax=ax, scatter_kws={'s':1})
                ax.text(x=np.amin(x), y=np.amax(y), va='top', ha='left',
                        s='p=%.2f' % (np.corrcoef(x, y)[0, 1]))
                if qoi_i_2 == num_qois_2-1:
                    ax.set_xlabel(qoi_names_1[qoi_i_1])
                if qoi_i_1 == 0:
                    ax.set_ylabel(qoi_names_2[qoi_i_2])
                if qoi_names_1[qoi_i_1] == 'EDVL':
                    ax.set_xlim(145, 160)
                elif qoi_names_1[qoi_i_1] == 'ESVL':
                    ax.set_xlim(60, 90)
                elif qoi_names_1[qoi_i_1] == 'LVEF':
                    ax.set_xlim(35, 60)
                elif qoi_names_1[qoi_i_1] == 'PmaxL':
                    ax.set_xlim(11, 13)
                if qoi_names_2[qoi_i_2] == 'qt_dur_mean':
                    ax.set_ylim(0.48, 0.50)
                elif qoi_names_2[qoi_i_2] == 't_pe_mean':
                    ax.set_ylim(0.28, 0.29)
        plt.show()


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

