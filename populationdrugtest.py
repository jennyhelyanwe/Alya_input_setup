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
        self.population_sf_filenames = []
        for population_i in range(self.population_size):
            data = pd.read_csv(self.personalised_population_dir + self.name + '_nodefield_personalisation-biomarker_' + str(
                population_i) + '.csv')
            for sf_name in population_sf_names:
                filename = self.personalised_population_dir + self.name + '_populationid_' + str(population_i) + '.' + sf_name
                self.population_sf_filenames.append(filename)
                with open(filename, 'w') as f:
                    for i in range(data[sf_name].values.shape[0]):
                        f.write(str(data[sf_name].index[i] + 1) + '\t' + str(data[sf_name].values[i]) + '\n')
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
            for sf_name in population_sf_names:
                cmd = 'cp ' + self.personalised_population_dir + self.name + '_populationid_' + str(
                    population_i) + '.' + sf_name + ' ' + self.alya_format.output_dir + 'heart.' + sf_name
                os.system(cmd)
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
        linewidth = 1.0
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
        titles = ['I', 'II', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']
        # baseline_dir ='/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/alya_simulations/DTI004_baseline_simulation_ep_DTI004_fine/'
        # clinical_ecg = np.loadtxt(baseline_dir + 'DTI004_clinical_full_ecg.txt', delimiter=',')
        # max_leads = np.amax(clinical_ecg)
        # clinical_ecg = clinical_ecg[:, 100:] / max_leads
        # Plot clinical ECGs
        # clinical_ecg = np.loadtxt(simulation_dir + 'baseline_population/' + 'populationid_' + str(0) + '_baseline' \
        #                       '_' + self.name + '/DTI004_clinical_full_ecg.txt',delimiter=',')
        # max_leads = np.amax(clinical_ecg)
        # dict = {'I': [], 'II': [], 'V1': [], 'V2': [], 'V3': [], 'V4': [], 'V5': [], 'V6': []}
        # dict['I'] = clinical_ecg[0, 100:] / max_leads
        # dict['II'] = clinical_ecg[1, 100:] / max_leads
        # dict['V1'] = clinical_ecg[2, 100:] / max_leads
        # dict['V2'] = clinical_ecg[3, 100:] / max_leads
        # dict['V3'] = clinical_ecg[4, 100:] / max_leads
        # dict['V4'] = clinical_ecg[5, 100:] / max_leads
        # dict['V5'] = clinical_ecg[6, 100:] / max_leads
        # dict['V6'] = clinical_ecg[7, 100:] / max_leads
        # clinical_t = np.arange(0, len(dict['I']), 1)
        # keys = list(dict.keys())
        # for lead_i in range(nb_leads):
        #     axes[lead_i].plot(clinical_t, dict[keys[lead_i]], 'g', linewidth=linewidth*2, label='Clinical')

        # Plot baseline popluation
        for population_i in range(self.population_size):
            current_sim_dir = simulation_dir + 'baseline_population/' + 'populationid_' + str(population_i) + '_baseline' \
                              '_' + self.name + '/'
            if os.path.exists(current_sim_dir + '/heart.exm.vin'):
                ecgs = vis._read_ECG(current_sim_dir + '/heart')
                ecgs['ts'][beat - 1] = ecgs['ts'][beat - 1] * 1000  # Units of milliseconds
                for lead_i in range(nb_leads):
                    if population_i == 0:
                        if lead_i <2:
                            axes[lead_i].plot(ecgs['ts'][beat - 1],
                                              ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_limb_leads'],
                                              'k', linewidth=linewidth, label='Baseline')
                        else:
                            axes[lead_i].plot(ecgs['ts'][beat - 1],
                                              ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_all_leads'],
                                              'k', linewidth=linewidth, label='Baseline')
                    else:
                        if lead_i < 2:
                            axes[lead_i].plot(ecgs['ts'][beat - 1],
                                              ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_limb_leads'],
                                              'k', linewidth=linewidth)
                        else:
                            axes[lead_i].plot(ecgs['ts'][beat - 1],
                                              ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_all_leads'],
                                              'k', linewidth=linewidth)
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
                    ecgs['ts'][beat - 1] = ecgs['ts'][beat - 1] * 1000  # Units of milliseconds
                    for lead_i in range(nb_leads):
                        if population_i == 0:
                            if lead_i < 2:
                                axes[lead_i].plot(ecgs['ts'][beat - 1],
                                                  ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_limb_leads'],
                                                  colours[dose_i], linewidth=linewidth,
                                                  label='Dose :' + str(dose_i + 1))
                            else:
                                axes[lead_i].plot(ecgs['ts'][beat - 1],
                                                  ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_all_leads'],
                                                  colours[dose_i], linewidth=linewidth, label='Dose :'+str(dose_i+1))
                        else:
                            if lead_i < 2:
                                axes[lead_i].plot(ecgs['ts'][beat - 1],
                                                  ecgs[lead_names[lead_i]][beat - 1] / ecgs['max_limb_leads'],
                                                  colours[dose_i], linewidth=linewidth)
                            else:
                                axes[lead_i].plot(ecgs['ts'][beat-1], ecgs[lead_names[lead_i]][beat-1]/ecgs['max_all_leads'],
                                              colours[dose_i], linewidth=linewidth)
                        axes[lead_i].set_title(titles[lead_i], fontsize=20)
                        axes[lead_i].set_ylim([-1.0, 1.0])
                        axes[lead_i].set_xlim([0, 500])
                        for tick in axes[lead_i].xaxis.get_major_ticks():
                            tick.label1.set_fontsize(14)
                        for tick in axes[lead_i].yaxis.get_major_ticks():
                            tick.label1.set_fontsize(14)
                        # if (dose_i==0) & (lead_i == 0) & (population_i == 0):
                        #     lines.append(line)
                        #     legends.append('Dose: '+str(dose_i))
        axes[7].legend()
        plt.show()


    def evaluate_drug_effect_ecgs(self, drug_name, drug_doses, simulation_dir):
        CL = 1.0
        vis = ECGPV_visualisation(CL)
        beat = 1
        # Plot baseline population
        mean_qt_effect = np.zeros((len(drug_doses), self.population_size))
        # std_qt_effect = np.zeros((len(drug_doses)))
        mean_tpe_effect = np.zeros((len(drug_doses), self.population_size))
        # std_tpe_effect = np.zeros((len(drug_doses)))
        mean_tpeak_effect = np.zeros((len(drug_doses), self.population_size))
        concentrations = np.array([0.5, 1, 2, 3, 4, 5, 6])  # nM dofetilide
        population_concentrations = np.zeros((len(drug_doses), self.population_size))
        # if not os.path.exists(simulation_dir + 'ecg_drug_effects.json'):
        for dose_i in range(len(drug_doses)):
            # For each dose of the drug run a population of models
            dose_dir = simulation_dir + 'dose_' + str(dose_i) + '/'
            print('Evaluating T wave biomarkers for dose ' + str(dose_i))
            for population_i in range(self.population_size):
                self.alya_format.simulation_dir = dose_dir
                current_sim_dir = dose_dir + 'populationid_' + str(population_i) + '_' + drug_name + \
                                  '_dose_' + str(dose_i) + '_' + self.name + '/'
                if os.path.exists(current_sim_dir + '/heart.exm.vin'):
                    ecgs = vis._read_ECG(current_sim_dir + '/heart')
                    analysis = vis.analysis_ECG_5leads(ecgs=ecgs, beat=beat, show=False)
                    mean_qt_effect[dose_i, population_i] = analysis['QT'][0] * 1000. # [ms]
                    mean_tpe_effect[dose_i, population_i] = analysis['T_pe_dur'][0] * 1000. # [ms]
                    mean_tpeak_effect[dose_i, population_i] = analysis['T_amp'][0]
                    population_concentrations[dose_i, population_i] = concentrations[dose_i]
                else:
                    print('Simulation results missing for: ', current_sim_dir)
        sim_ecg_drug_effects = {'concentrations':[0.5, 1, 2, 3, 4, 5, 6],
                     'QT': mean_qt_effect,
                     'QT_mean': list(np.mean(mean_qt_effect, axis=1)),
                     'QT_std': list(np.std(mean_qt_effect, axis=1)),
                     'Tpe': mean_tpe_effect,
                     'Tpe_mean': list(np.mean(mean_tpe_effect, axis=1)),
                     'Tpe_std': list(np.std(mean_tpe_effect, axis=1)),
                     'Tpeak': mean_tpeak_effect,
                     'Tpeak_mean': list(np.mean(mean_tpeak_effect, axis=1)),
                     'Tpeak_std': list(np.std(mean_tpeak_effect, axis=1))}
        #     # with open(simulation_dir + 'ecg_drug_effects.json', 'w') as f:
        #     #     json.dump(sim_ecg_drug_effects, f)
        # else:
        #     sim_ecg_drug_effects = json.load(open(simulation_dir + 'ecg_drug_effects.json', 'r'))

        # Compare effects against average cohort biomarkers
        concentrations = np.array([0.5, 1, 2, 3, 4, 5, 6]) # nM dofetilide
        cohort_qtc = np.array([390, 388, 388, 410, 440, 463, 442]) # ms
        cohort_qtc_std = np.array([27, 22, 11, 28, 35, 45, 31]) # ms
        cohort_tpe = np.array([76, 75, 74, 89, 104, 121, 245]) # ms
        cohort_tpe_std = np.array([5, 7, 6, 24, 10, 41, 32]) # ms
        cohort_tpeak = np.array([680, 703, 673, 672, 548, 475, 605]) # mV
        cohort_tpeak_std = np.array([154, 158, 185, 133, 178, 229, 100]) # mV
        colours = ['#fde725', '#90d743', '#35b779', '#21918c', '#31688e', '#443983', '#440154']
        fig = plt.figure()
        gs = GridSpec(1,1)
        axis1 = fig.add_subplot(gs[0,0])
        # axis2 = fig.add_subplot(gs[0,1])
        axis1.plot(concentrations, cohort_qtc, colours[0])
        axis1.fill_between(concentrations, cohort_qtc - cohort_qtc_std, cohort_qtc + cohort_qtc_std, alpha=0.5,
                           edgecolor=None, facecolor=colours[0], label='Clinical QT ranges')
        axis1.plot(np.reshape(population_concentrations, population_concentrations.shape[0]*population_concentrations.shape[1]),
                   np.reshape(sim_ecg_drug_effects['QT'], sim_ecg_drug_effects['QT'].shape[0]*sim_ecg_drug_effects['QT'].shape[1]),
                   'b.', label='Simulated QT')
        # top = np.array(sim_ecg_drug_effects['QT_mean']) - np.array(sim_ecg_drug_effects['QT_std'])
        # bottom = np.array(sim_ecg_drug_effects['QT_mean']) + np.array(sim_ecg_drug_effects['QT_std'])
        # axis1.fill_between(concentrations, top, bottom, alpha=0.5,
        #                    edgecolor=None, facecolor=colours[2], label='Simulated QT')

        axis1.plot(concentrations, cohort_tpe, 'cyan')
        axis1.fill_between(concentrations, cohort_tpe - cohort_tpe_std, cohort_tpe + cohort_tpe_std, alpha=0.5,
                           edgecolor=None, facecolor='cyan', label='Clinical Tpe ranges')
        axis1.plot(np.reshape(population_concentrations, population_concentrations.shape[0]*population_concentrations.shape[1]),
                   np.reshape(sim_ecg_drug_effects['Tpe'], sim_ecg_drug_effects['Tpe'].shape[0]*sim_ecg_drug_effects['Tpe'].shape[1]),
                   'k.', label='Simulated Tpe')
        # top = np.array(sim_ecg_drug_effects['Tpe']) - np.array(sim_ecg_drug_effects['Tpe_std'])
        # bottom = np.array(sim_ecg_drug_effects['Tpe']) + np.array(sim_ecg_drug_effects['Tpe_std'])
        # axis1.fill_between(concentrations, top, bottom, alpha=0.5,
        #                    edgecolor=None, facecolor=colours[5], label='Simulated Tpe')
        axis1.legend()
        plt.show()

    def extract_vms_for_download(self, simulation_dir, dir_for_download, drug_name, drug_doses):
        # Sort baseline CSVs
        print('Saving baseline CSVs to ' + dir_for_download)
        if not os.path.exists(dir_for_download):
            os.mkdir(dir_for_download)
        destination_root = dir_for_download + '/baseline_population/'
        if not os.path.exists(destination_root):
            os.mkdir(destination_root)
        for population_i in range(self.population_size):
            print('Processing population id: ' + str(population_i))
            csv_dir = simulation_dir + 'baseline_population/' + 'populationid_' + str(population_i) + '_baseline' \
                          '_' + self.name + '/results_csv/'
            destination_dir = destination_root + str(population_i) + '/'
            if not os.path.exists(destination_dir):
                os.mkdir(destination_dir)
            if os.path.exists(csv_dir):
                os.system('cp ' + csv_dir + '* '+destination_dir)
                filenames = os.listdir(csv_dir)
                filenames_shared = pymp.shared.list(filenames)
                threadsNum = multiprocessing.cpu_count()
                with pymp.Parallel(min(threadsNum, len(filenames))) as p1:
                    for conf_i in p1.range(len(filenames)):
                        # for file in filenames:
                        file = filenames_shared[conf_i]
                        os.system('cp ' + csv_dir + file + ' ' + destination_dir + '/' + file.replace('heart', self.name))
        # Sort dosage CSVs
        print('Saving drug simulation CSVs to ' + dir_for_download)
        if not os.path.exists(dir_for_download + drug_name):
            os.mkdir(dir_for_download + drug_name)
        for dose_i in range(len(drug_doses)):
            destination_root = dir_for_download + drug_name + '/dose_' + str(dose_i) + '/'
            if not os.path.exists(destination_root):
                os.mkdir(destination_root)
            # For each dose of the drug run a population of models
            dose_dir = simulation_dir + 'dose_' + str(dose_i) + '/'
            for population_i in range(self.population_size):
                print('Processing population id: ' + str(population_i))
                csv_dir = dose_dir + 'populationid_' + str(population_i) + '_' + drug_name + '_dose_'+str(dose_i)+'_' + self.name + '/results_csv/'
                destination_dir = destination_root + str(population_i) + '/'
                if not os.path.exists(destination_dir):
                    os.mkdir(destination_dir)
                if os.path.exists(csv_dir):
                    os.system('cp ' + csv_dir + '* ' + destination_dir)
                    # filenames = os.listdir(csv_dir)
                    # filenames_shared = pymp.shared.array((len(filenames)), dtype=str)
                    # filenames_shared[:] = filenames
                    # threadsNum = multiprocessing.cpu_count()
                    # with pymp.Parallel(min(threadsNum, len(filenames))) as p1:
                    #     for conf_i in p1.range(len(filenames)):
                    # # for file in filenames:
                    #         file = filenames_shared[conf_i]
                    #         os.system('cp ' + csv_dir + file + ' ' + destination_dir + '/' + file.replace('heart', self.name))





