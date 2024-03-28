import matplotlib
from matplotlib.ticker import MaxNLocator
matplotlib.use('tkagg')
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from healthy_qoi_ranges import HealthyBiomarkerRanges

########################################################################################################################
data = pd.read_csv('SA_summary_OAT_corrs.csv').values
correlation_matrix = data[:, 1:]
data = pd.read_csv('SA_summary_OAT_ranges.csv').values
ranges_matrix = data[:, 1:]
ranges_matrix_normalised = ranges_matrix
maxima = np.nanmax(ranges_matrix, axis=0)
ranges_matrix_normalised[:, np.nonzero(maxima)] = ranges_matrix[:, np.nonzero(maxima)] / maxima[np.nonzero(maxima)]
qoi_names = pd.read_csv('SA_summary_OAT_corrs.csv').columns.values[1:]
qoi = {}
qoi_ticks = []
for i in range(len(qoi_names)):
    qoi[qoi_names[i]] = len(qoi_names) - i
    qoi_ticks.append(len(qoi_names)-i)

param_names = data[:,0]
param = {}
param_ticks = []
for i in range(len(param_names)):
    param[param_names[i]] = len(param_names) - i
    param_ticks.append(len(param_names)-i)

# #######################################################################################################################
# # All together plot
# ax = plt.figure(figsize=[8,10]).gca()
# for param_i in range(len(param_names)):
#     for qoi_i in range(len(qoi_names)):
#         if correlation_matrix[param_i, qoi_i] < 0.0:
#             c = 'b'
#         elif correlation_matrix[param_i, qoi_i] > 0.0:
#             c = 'r'
#         else:
#             continue
#         # if ranges_matrix_normalised[param_i, qoi_i] > 0.1:
#         w = ranges_matrix_normalised[param_i, qoi_i]*1.5
#         # else:
#         #     w = 0
#         if abs(correlation_matrix[param_i, qoi_i]) > 0.2:
#             plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
#                      alpha=abs(correlation_matrix[param_i, qoi_i]), color=c, linewidth=w)
#
# ax.set_xticklabels([])
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# ax.set_yticks(param_ticks)
# ax.set_yticklabels(param_names)
# ax.set_xlim([0, 1])
# ax.set_ylim([0, len(qoi_names)])
# secax = ax.twinx()
# secax.set_ylim(0, len(qoi_names))
# secax.set_yticks(qoi_ticks)
# secax.set_yticklabels(qoi_names)
# plt.tight_layout()
# plt.show()
########################################################################################################################
# # Highlight separate groups of parameters
# # group = ['INaL', 'IKr', 'INaK', 'ICaL', 'ISERCA', 'Ca50', 'K_ws']
# group = ['Kct', 'a', 'af', 'as', 'afs', 'Kepi', 'Tref_scaling']
# # group = ['R_lv', 'C_lv', 'Gain_IVR_lv', 'dGain_IVR_lv', 'P_ejection_lv']
# ax = plt.figure(figsize=[8,10]).gca()
# for param_i in range(len(param_names)):
#     for qoi_i in range(len(qoi_names)):
#         if correlation_matrix[param_i, qoi_i] < 0.0:
#             c = 'b'
#         elif correlation_matrix[param_i, qoi_i] > 0.0:
#             c = 'r'
#         else:
#             continue
#         # if ranges_matrix_normalised[param_i, qoi_i] > 0.1:
#         w = ranges_matrix_normalised[param_i, qoi_i]*1.5
#         # else:
#         #     w = 0
#         if abs(correlation_matrix[param_i, qoi_i]) > 0.2:
#             if param_i in group:
#                 plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
#                          alpha=abs(correlation_matrix[param_i, qoi_i]), color=c, linewidth=w)
#             else:
#                 plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
#                          alpha=0.1, color='k', linewidth=w)
# ax.set_xticklabels([])
# ax.xaxis.set_major_locator(MaxNLocator(integer=True))
# ax.set_yticks(param_ticks)
# ax.set_yticklabels(param_names)
# ax.set_xlim([0, 1])
# ax.set_ylim([0, len(qoi_names)])
# secax = ax.twinx()
# secax.set_ylim(0, len(qoi_names))
# secax.set_yticks(qoi_ticks)
# secax.set_yticklabels(qoi_names)
# plt.tight_layout()
# plt.show()

########################################################################################################################
# Highlight separate groups of QoIs
# group = ['QT', 'Tpe', 'Tamp']
group = ['LVEF']
# group = ['LVEDV', 'LVESV','SVL','LVESP','LVEF', 'dPdt', 'PER', 'PFR']
# group = ['LVEDV', 'LVESV', 'SVL', 'LVESP']
# group = ['RVEDV', 'RVESV', 'SVR', 'RVESP']
# group = ['LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVSV', 'RVSV', 'LVESP', 'RVESP', 'LVEF', 'dPdt', 'PER', 'PFR']
# group = ['dThickness', 'Peak AVPD', 'Apex displacement']
# group = ['Peak Eff', 'Peak Ell', 'Peak Err']
ax = plt.figure(figsize=[8,10]).gca()
for param_i in range(len(param_names)):
    for qoi_i in range(len(qoi_names)):
        if correlation_matrix[param_i, qoi_i] < 0.0:
            c = 'b'
        elif correlation_matrix[param_i, qoi_i] > 0.0:
            c = 'r'
        else:
            continue
        # if ranges_matrix_normalised[param_i, qoi_i] > 0.1:
        w = ranges_matrix_normalised[param_i, qoi_i]*1.5
        # else:
        #     w = 0
        if abs(correlation_matrix[param_i, qoi_i]) > 0.2:
            if qoi_names[qoi_i] in group:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=abs(correlation_matrix[param_i, qoi_i]), color=c, linewidth=w)
            else:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=0.1, color='k', linewidth=w)
ax.set_xticklabels([])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_yticks(param_ticks)
ax.set_yticklabels(param_names)
ax.set_xlim([0, 1])
ax.set_ylim([0, len(qoi_names)])
secax = ax.twinx()
secax.set_ylim(0, len(qoi_names))
secax.set_yticks(qoi_ticks)
secax.set_yticklabels(qoi_names)
plt.tight_layout()
plt.show()




