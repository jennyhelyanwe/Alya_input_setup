from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib import ProblemSpec
from SALib.plotting.bar import plot as barplot
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.gridspec import GridSpec
import os
import pandas as pd
import scipy

protocol = 'UQ'
if protocol == 'SA':
    flag = 'sample'
    parameter_names = ['INa', 'INaL', 'Ito', 'IKr', 'IKs', 'IK1', 'INCX', 'INaK', 'ICaL', 'Jrel', 'Jup', 'ca50', 'kuw', 'kws', 'Istim']
    sp = ProblemSpec({'num_vars': len(parameter_names),
                   'names': parameter_names,
                   'bounds': [[0.5, 2]]*len(parameter_names),
                   'outputs': ['APD40', 'APD50', 'APD90', 'CTD50', 'CTD90', 'CaTmax', 'CaTmin', 'Tamax', 'Tamin', 'TaD50', 'TaD90', 'dTadtmax'] #, 'APD90', 'CTD50', 'CTD90', 'CaTmin', 'CaTmax', 'Tamax', 'Tamin'
                   })
    qoi_names = sp.get('outputs')
    # if flag == 'sample':
    sample_size = 2 ** 7
    sp.sample_sobol(sample_size, calc_second_order=True)
    param_values= sp.samples
    print('Number of samples: ', param_values.shape[0])
    output_dir = 'sa_cell_em_n_' + str(param_values.shape[0]) + '/'
    # np.save('sa_param_values.txt', param_values, delimiter=',')
    # with open('sa_output_dir.txt', 'w') as f:
    #     f.write(output_dir)
    #
    ########################################################################################################################
    # Run simulations in MATLAB
    # os.system('matlab -nodisplay -r "evaluateSA"')

    ########################################################################################################################
    # elif flag == 'analyse':
    qoi_names_ta = ['APD90', 'Tamax', 'TaD50', 'dTadtmax']  # qoi_names[7:]
    qoi_names_with_units_ta = ['APD90 (ms)', 'Tamax (kPa)', 'TaD50 (ms)', 'dTadtmax (kPa/ms)']
    sp['outputs'] = qoi_names_ta
    celltype = ['endo', 'epi', 'mid']
    for cell in celltype:
        # Y = np.loadtxt(output_dir + cell + '_output.txt', float, delimiter=',', skiprows=1, header=1)
        print('Reading ' + output_dir + cell + '_output.txt')
        Y = pd.read_csv(output_dir + cell + '_output.txt')
        Y_qoi = Y[qoi_names_ta].values
        X = np.loadtxt('sa_param_values.txt', delimiter=',')
        ####################################################################################################################
        # Scatter plots with correlation coefficients
        fig = plt.figure(tight_layout=True, figsize=(18,10))
        gs = GridSpec(len(qoi_names_with_units_ta), 8)
        axes = [[0] * 8] * len(qoi_names_with_units_ta)
        for qoi_i in range(len(qoi_names_with_units_ta)):
            for param_j in range(8):
                axes[qoi_i][param_j] = fig.add_subplot(gs[qoi_i, param_j])
                x = X[:,param_j]
                y = Y_qoi[:,qoi_i]
                sns.regplot(x=x, y=y, ax=axes[qoi_i][param_j], scatter_kws={'s':1})
                axes[qoi_i][param_j].text(x=np.amin(x), y=np.amax(y), va='top', ha='left',
                                          s='p=%.2f' % (np.corrcoef(x,y)[0,1]))
                if qoi_i == len(qoi_names_with_units_ta)-1:
                    axes[len(qoi_names_with_units_ta)-1][param_j].set_xlabel(parameter_names[param_j])
            axes[qoi_i][0].set_ylabel(qoi_names_with_units_ta[qoi_i])
        plt.savefig(output_dir + cell + '_scatter_plot_part1.png')
        plt.show()
        fig = plt.figure(tight_layout=True, figsize=(18,10))
        gs = GridSpec(len(qoi_names_with_units_ta), 7)
        axes = [[0] * 7] * len(qoi_names_with_units_ta)
        for qoi_i in range(len(qoi_names_with_units_ta)):
            for param_j in range(7):
                axes[qoi_i][param_j] = fig.add_subplot(gs[qoi_i, param_j])
                x = X[:,param_j+8]
                y = Y_qoi[:,qoi_i]
                sns.regplot(x=x, y=y, ax=axes[qoi_i][param_j], scatter_kws={'s':1})
                axes[qoi_i][param_j].text(x=np.amin(x), y=np.amax(y), va='top', ha='left', s='p=%.2f' % (np.corrcoef(x,y)[0,1]))
                if qoi_i == len(qoi_names_with_units_ta)-1:
                    axes[len(qoi_names_with_units_ta)-1][param_j].set_xlabel(parameter_names[param_j + 8])
            axes[qoi_i][0].set_ylabel(qoi_names_with_units_ta[qoi_i])
        plt.savefig(output_dir + cell + '_scatter_plot_part2.png')
        plt.show()
        ###################################################################################################################
        # Tornado plot of sensitivity indices https://seaborn.pydata.org/examples/part_whole_bars.html
        sp.set_results(Y_qoi)
        Si = sp.analyze_sobol(print_to_console=False, calc_second_order=True)
        analysis = sp._analysis
        QOIs = sp.get('outputs')
        fig = plt.figure(tight_layout=True, figsize=(15,6))
        gs = GridSpec(1, len(qoi_names_with_units_ta))
        data = Si.to_df()
        sns.set_theme(style="whitegrid")
        for qoi_i in range(len(qoi_names_with_units_ta)):
            ax = fig.add_subplot(gs[0, qoi_i])
            st_s1_data = pd.concat([data[qoi_i][0], data[qoi_i][1]], axis=1)
            sorted_data = st_s1_data.reindex(st_s1_data.abs().sort_values('S1', ascending=False).index)
            sorted_parameter_names = []
            for row in sorted_data.index:
                sorted_parameter_names.append(row)
            sns.set_color_codes("pastel")
            sns.barplot(data=sorted_data, x='ST', y=sorted_parameter_names, label='ST', color='b')
            sns.set_color_codes("muted")
            sns.barplot(data=sorted_data, x='S1', y=sorted_parameter_names, label='S1', color='b')

            ax.set( ylabel="",
                   xlabel=qoi_names_with_units_ta[qoi_i])
            if qoi_i == 3:
                ax.legend(ncol=2, loc="lower right", frameon=True)
        plt.savefig(output_dir + cell+'_ST_S1_tornado.png')
        plt.show()
        ####################################################################################################################
        # Second order effects heat map
        # sp.set_results(Y)
        # Si = sp.analyze_sobol(print_to_console=False, calc_second_order=True)
        # analysis = sp._analysis
        # QOIs = sp.get('outputs')
        # fig = plt.figure()
        # names = []
        # data = Si.to_df()
        # for row in data.index:
        #     names.append(row)
        # for qoi_i in range(4):
        #     print(cell +'_' +QOIs[qoi_i+7])
        #     ax = sns.heatmap(analysis[QOIs[qoi_i+7]]['S2'])
        #     ax.set_yticklabels(names, rotation=45, fontsize=14)
        #     ax.set_xticklabels(names, rotation=45,
        #                        horizontalalignment='right', fontsize=14)
        #     ax.set_title(cell +'_' +QOIs[qoi_i+7])
        #     plt.savefig(output_dir + cell +'_' + QOIs[qoi_i+7]+'_S2_heatmap.png')
        # plt.show()

    #print(sp[0][2]['S2'])
    # ax = sns.heatmap(analysis[0]['S2']) #, vmin=-1, vmax=1, center=0,
                     # square=True, annot_kws={"size": 20, "color": "k"}) # cmap=sns.diverging_palette(220, 20, n=200),
    # for (j, i), label in np.ndenumerate(sp.get('names')):
    #     ax.text(i + 0.5, j + 0.5, np.round(label, 2),
    #             fontdict=dict(ha='center', va='center', color='k', fontsize=10))
    # sns.pairplot(dataframe_ecg_biomarkers)

    #sp.plot()
    # sp.plot()
    # plt.savefig('Overall_plot.png')
run = False
visualise = True
style = '8_traces' # '8_traces' or 'min_max'
if protocol == 'UQ':
    ########################################################################################################################
    # Uncertainty quantification - using 0.5 and 2 fold as ranges
    print('Propagating prescribed uncertainty of the following parameters and quantifying effect on QoIs')
    uncertain_parameters = ['ca50', 'ICaL', 'kws', 'Jup', 'INaL']
    parameter_names = ['INa', 'INaL', 'Ito', 'IKr', 'IKs', 'IK1', 'INCX', 'INaK', 'ICaL', 'Jrel', 'Jup', 'ca50', 'kuw', 'kws',
             'Istim']
    for uq_param in uncertain_parameters:
        if style == 'min_max':
            print('Uncertainty propagation for ', uq_param)
            param_values = np.ones((3, len(parameter_names)))
            param_values[0, parameter_names.index(uq_param)] = 0.5
            param_values[1, parameter_names.index(uq_param)] = 1.0
            param_values[2, parameter_names.index(uq_param)] = 2.0
        elif style == '8_traces':
            print('Uncertainty propagation for ', uq_param)
            param_values = np.ones((8, len(parameter_names)))
            for i in range(8):
                param_values[i, parameter_names.index(uq_param)] = np.linspace(0.5, 2, 8)[i]
        output_dir = 'uq_cell_em_' + style + '_' + uq_param + '/'
        if run:
            np.savetxt('uq_param_values.txt', param_values, delimiter=',')
            with open('uq_output_dir.txt', 'w') as f:
                f.write(output_dir)
            # Run simulations in MATLAB
            os.system('matlab -nodisplay -r "evaluateUQ"')
        if visualise:
            print('Visualising...')
            if style == 'min_max':
                # Load outputs
                V = scipy.io.loadmat(output_dir + 'V.mat')['V'][0]
                cai = scipy.io.loadmat(output_dir + 'cai.mat')['cai'][0]
                Ta = scipy.io.loadmat(output_dir + 'Ta.mat')['Ta'][0]
                time = scipy.io.loadmat(output_dir + 'time.mat')['time'][0]
                fig = plt.figure(tight_layout=True, figsize=(15, 6))
                gs = GridSpec(1, 3)
                ax_V = fig.add_subplot(gs[0, 0])
                ax_V.plot(time[1][0], V[1][0], 'k')
                ax_V.plot(time[0][0], V[0][0], 'b')
                ax_V.plot(time[2][0], V[2][0], 'r')
                ax_V.fill(np.append(time[0][0], time[2][0][::-1]), np.append(V[0][0], V[2][0][::-1]), alpha=0.3, edgecolor=None, color='k')

                ax_cai = fig.add_subplot(gs[0, 1])
                ax_cai.plot(time[1][0], cai[1][0], 'k')
                ax_cai.plot(time[0][0], cai[0][0], 'b')
                ax_cai.plot(time[2][0], cai[2][0], 'r')
                ax_cai.fill(np.append(time[0][0], time[2][0][::-1]),
                          np.append(cai[0][0], cai[2][0][::-1]), alpha=0.3, edgecolor=None, color='k')

                ax_ta = fig.add_subplot(gs[0, 2])
                ax_ta.plot(time[1][0], Ta[1][0], 'k')
                ax_ta.plot(time[0][0], Ta[0][0], 'b')
                ax_ta.plot(time[2][0], Ta[2][0], 'r')
                ax_ta.fill(np.append(time[0][0], time[2][0][::-1]),
                          np.append(Ta[0][0], Ta[2][0][::-1]), alpha=0.3, edgecolor=None, color='k')
                plt.savefig(output_dir + 'uq_' + uq_param + '_V_cai_Ta_ranges_fill.png')
            elif style == '8_traces':
                celltypes = ['endo', 'epi', 'mid']
                for celltype in celltypes:
                    # Load outputs
                    V = scipy.io.loadmat(output_dir + 'V_'+ celltype+'.mat')['V'][0]
                    cai = scipy.io.loadmat(output_dir + 'cai_'+ celltype+'.mat')['cai'][0]
                    Ta = scipy.io.loadmat(output_dir + 'Ta_'+ celltype+'.mat')['Ta'][0]
                    time = scipy.io.loadmat(output_dir + 'time_'+ celltype+'.mat')['time'][0]
                    fig = plt.figure(tight_layout=True, figsize=(15, 6))
                    gs = GridSpec(1, 3)
                    ax_V = fig.add_subplot(gs[0, 0])
                    for i in range(8):
                        ax_V.plot(time[i][0], V[i][0])
                    ax_cai = fig.add_subplot(gs[0, 1])
                    for i in range(8):
                        ax_cai.plot(time[i][0], cai[i][0])
                    ax_Ta = fig.add_subplot(gs[0, 2])
                    for i in range(8):
                        ax_Ta.plot(time[i][0], Ta[i][0])
                    print('Saving figure to: '+ output_dir + 'uq_' + uq_param + '_V_cai_Ta_8_traces_'+ celltype+'.png' )
                    plt.savefig(output_dir + 'uq_' + uq_param + '_V_cai_Ta_8_traces_'+ celltype+'.png')
                    plt.show()
