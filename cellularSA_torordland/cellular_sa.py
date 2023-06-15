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
    sp.analyze_sobol(print_to_console=True, calc_second_order=True)
    sp.heatmap()

    sp.plot()
    plt.show()

