from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib import ProblemSpec
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns

flag = 'analyse'
sp = ProblemSpec({'num_vars': 16,
               'names': ['INa', 'INaL', 'Ito', 'IKr', 'IKs', 'IK1', 'INCX', 'INaK', 'ICaL', 'Jrel', 'Jup', 'ca50', 'kuw', 'kws', 'ksu', 'Istim'],
               'bounds': [[0.5, 2]]*16,
               'outputs': ['APD40', 'APD50', 'APD90', 'CTD50', 'CTD90', 'CaTmax', 'CaTmin', 'Tamax', 'Tamin', 'TaD50', 'TaD90'] #, 'APD90', 'CTD50', 'CTD90', 'CaTmin', 'CaTmax', 'Tamax', 'Tamin'
               })
if flag == 'sample':
    sample_size = 4
    sp.sample_sobol(sample_size, calc_second_order=True)
    param_values= sp.samples
    print ('Number of samples: ', param_values.shape[0])
    np.savetxt('param_values.txt', param_values, delimiter=',')
elif flag == 'analyse':
    Y = np.loadtxt('endo_output.txt', float, delimiter=',', skiprows=1)
    print(Y.shape)
    sp.set_results(Y)
    sp.analyze_sobol(print_to_console=False, calc_second_order=True)
    analysis = sp._analysis
    for str in sp.get('outputs'):
        print(str)
        print(analysis[str]['S2'])
        ax = sns.heatmap(analysis[str]['S2'])
        ax.set_yticklabels(sp.get('names'), rotation=45, fontsize=14)
        ax.set_xticklabels(sp.get('names'), rotation=45,
                           horizontalalignment='right', fontsize=14)
        ax.set_title(str)
        plt.show()
        plt.savefig(str+'_S2.png')
        sp.plot()
        plt.savefig(str+'_S1.png')
        plt.show()
    #print(sp[0][2]['S2'])
    # ax = sns.heatmap(analysis[0]['S2']) #, vmin=-1, vmax=1, center=0,
                     # square=True, annot_kws={"size": 20, "color": "k"}) # cmap=sns.diverging_palette(220, 20, n=200),
    # for (j, i), label in np.ndenumerate(sp.get('names')):
    #     ax.text(i + 0.5, j + 0.5, np.round(label, 2),
    #             fontdict=dict(ha='center', va='center', color='k', fontsize=10))
    # sns.pairplot(dataframe_ecg_biomarkers)

    #sp.plot()

