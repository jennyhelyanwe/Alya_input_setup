import matplotlib
from matplotlib.ticker import MaxNLocator
matplotlib.use('tkagg')
from matplotlib import pyplot as plt
import numpy as np


param_names = ['INaL', 'IKr', 'INaK', 'ICaL', 'ISERCA', 'Ca50', 'K_ws', 'Kct',
                'a', 'af', 'as', 'afs', 'Kepi', 'Tref_scaling', 'R_lv', 'C_lv',
                'Gain_IVR_lv', 'dGain_IVR_lv', 'P_ejection_lv']
param = {}
param_ticks = []
for i in range(len(param_names)):
    param[param_names[i]] = len(param_names) - i
    param_ticks.append(len(param_names)-i)

qoi_names = ['Tpe', 'Tamp', 'QT', 'LVEDV', 'LVESV', 'RVEDV', 'RVESV',
             'LVSV', 'RVSV', 'LVESP', 'RVESP', 'LVEF', 'dPdt', 'PER', 'PFR',
             'dThickness', 'Peak AVPD', 'Apex d', 'Peak Eff',
             'Peak Ell', 'Peak Err']
qoi = {}
qoi_ticks = []
for i in range(len(qoi_names)):
    qoi[qoi_names[i]] = len(qoi_names) - i
    qoi_ticks.append(len(qoi_names)-i)
correlations = [['Kct','Peak AVPD',  -0.42],
                ['afs','Peak AVPD', -0.55],
                ['Ca50', 'Peak AVPD', 0.43],
                ['Tref_scaling', 'Peak AVPD', 0.25],
                ['dGain_IVR_lv', 'Peak AVPD', 0.31],
                ['P_ejection_lv', 'Peak AVPD', -0.62],
                ['Kct', 'Apex d', -0.45],
                ['afs', 'Apex d', -0.51],
                ['Ca50', 'Apex d', 0.43],
                ['Tref_scaling', 'Apex d', 0.3],
                ['P_ejection_lv', 'LVEF', -0.59],
                ['Kct', 'LVEF', -0.41],
                ['afs', 'LVEF', -0.49],
                ['Ca50', 'LVEF', 0.47],
                ['Tref_scaling', 'LVEF', 0.29],
                ['ICaL', 'LVEF', 0.27],
                ['ISERCA', 'LVEF', -0.48],
                ['Kct', 'LVESP', -0.32],
                ['afs', 'LVESP', -0.32],
                ['Ca50', 'LVESP', 0.38],
                ['P_ejection_lv', 'LVESP', 0.71],
                ['dGain_IVR_lv', 'LVESP', -0.29],
                ['Kepi', 'LVESP', -0.21],
                ['ISERCA', 'LVESP', -0.42],
                ['INaK', 'RVESP', -0.27],
                ['ICaL', 'RVESP', 0.44],
                ['Kepi', 'PER', 0.24],
                ['Kct', 'PER', 0.42],
                ['afs', 'PER', 0.34],
                ['Ca50', 'PER', 0.34],
                ['Tref_scaling', 'PER', -0.25],
                ['ICaL', 'PER', 0.25],
                ['ISERCA', 'PER', -0.26],
                ['Kepi', 'PFR', 0.17],
                ['Kct', 'PFR', -0.19],
                ['afs', 'PFR', -0.33],
                ['Ca50', 'PFR', 0.34],
                ['Tref_scaling', 'PFR', 0.2],
                ['ICaL', 'PFR', 0.31],
                ['ISERCA', 'PFR', -0.47],
                ['R_lv', 'LVESV', 0.14],
                ['dGain_IVR_lv', 'LVESV', -0.3],
                ['Gain_IVR_lv', 'LVESV', 0.19],
                ['P_ejection_lv', 'LVESV', 0.59],
                ['afs', 'LVESV', 0.4],
                ['Kct', 'LVESV', 0.35],
                ['Ca50', 'LVESV', -0.37],
                ['Tref_scaling', 'LVESV', -0.2],
                ['ICaL', 'LVESV', -0.36],
                ['ISERCA', 'LVESV', 0.39],
                ['Kepi', 'LVEDV', -0.14],
                ['Kct', 'LVEDV', 0.49],
                ['a', 'LVEDV', -0.34],
                ['afs', 'LVEDV', -0.66],
                ['Ca50', 'LVEDV', 0.6],
                ['Tref_scaling', 'LVEDV', 0.5],
                ['ICaL', 'LVEDV', -0.43],
                ['ISERCA', 'LVEDV', -0.41],
                ['INaK', 'RVEDV', 0.25],
                ['ICaL', 'RVEDV', -0.49],
                ['ICaL', 'RVESV', -0.35],
                ['ISERCA', 'RVESV', 0.4],
                ['ICaL', 'LVSV', 0.23],
                ['ISERCA', 'LVSV', -0.51],
                ['ICaL', 'RVSV', 0.33],
                ['Kct', 'Tpe', -0.23],
                ['af', 'Tpe', -0.13],
                ['afs', 'Tpe', -0.33],
                ['Tref_scaling', 'Tpe', 0.16],
                ['INaL', 'Tpe', 0.21],
                ['IKr', 'Tpe', -0.29],
                ['ICaL', 'Tpe', 0.3],
                ['Kct', 'QT', -0.22],
                ['af', 'QT', -0.13],
                ['afs', 'QT', -0.33],
                ['Tref_scaling', 'QT', 0.15],
                ['INaL', 'QT', 0.21],
                ['IKr', 'QT', -0.29],
                ['ICaL', 'QT', 0.3],
                ['Kepi', 'Tamp',0.21],
                ['a', 'Tamp', -0.19],
                ['af', 'Tamp', 0.12],
                ['Tref_scaling', 'Tamp',0.14],
                ['INaL', 'Tamp',0.22],
                ['ICaL', 'Tamp',-0.22],
                ['ICaL', 'dPdt',0.35],
                ['ISERCA', 'dPdt',0.3],]

# # All together plot
# ax = plt.figure().gca()
# for i in range(len(correlations)):
#     if correlations[i][2] < 0.0:
#         c = 'r' # Negative correlation in red
#     else:
#         c = 'b' # Positive correlation in blue
#     plt.plot([0, 1], [param[correlations[i][0]], qoi[correlations[i][1]]], alpha=abs(correlations[i][2]), color=c)
#
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
#
#
# # Highlight separate groups of parameters
# group = ['INaL', 'IKr', 'INaK', 'ICaL', 'ISERCA', 'Ca50', 'K_ws']
# group = ['Kct', 'a', 'af', 'as', 'afs', 'Kepi', 'Tref_scaling']
# # group = ['R_lv', 'C_lv', 'Gain_IVR_lv', 'dGain_IVR_lv', 'P_ejection_lv']
# ax = plt.figure().gca()
# for i in range(len(correlations)):
#     if correlations[i][2] < 0.0:
#         c = 'r' # Negative correlation in red
#     else:
#         c = 'b' # Positive correlation in blue
#     if correlations[i][0] in group:
#         plt.plot([0, 1], [param[correlations[i][0]], qoi[correlations[i][1]]], alpha=abs(correlations[i][2]), color=c)
#     else:
#         plt.plot([0, 1], [param[correlations[i][0]], qoi[correlations[i][1]]], alpha=0.1, color='k')
#
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

# Highlight separate groups of qois
qoi_names = ['Tpe', 'Tamp', 'QT', 'LVEDV', 'LVESV', 'RVEDV', 'RVESV',
             'LVSV', 'RVSV', 'LVESP', 'RVESP', 'LVEF', 'dPdt', 'PER', 'PFR',
             'dThickness', 'Peak AVPD', 'Apex d', 'Peak Eff',
             'Peak Ell', 'Peak Err']
group = ['Tpe', 'Tamp', 'QT',]
# group = ['LVEDV', 'LVESV', 'RVEDV', 'RVESV', 'LVSV', 'RVSV', 'LVESP', 'RVESP', 'LVEF', 'dPdt', 'PER', 'PFR']
# group = ['dThickness', 'Peak AVPD', 'Apex d', 'Peak Eff', 'Peak Ell', 'Peak Err']
ax = plt.figure().gca()
for i in range(len(correlations)):
    if correlations[i][2] < 0.0:
        c = 'r' # Negative correlation in red
    else:
        c = 'b' # Positive correlation in blue
    if correlations[i][1] in group:
        plt.plot([0, 1], [param[correlations[i][0]], qoi[correlations[i][1]]], alpha=abs(correlations[i][2]), color=c)
    else:
        plt.plot([0, 1], [param[correlations[i][0]], qoi[correlations[i][1]]], alpha=0.1, color='k')


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
