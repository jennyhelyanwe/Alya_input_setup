import matplotlib
from matplotlib.ticker import MaxNLocator
matplotlib.use('tkagg')
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from healthy_qoi_ranges import HealthyBiomarkerRanges

########################################################################################################################
slopes = pd.read_csv('ALL_SA_OAT_slopes.csv').values[:, 1:]
p_values = pd.read_csv('ALL_SA_OAT_p_values.csv').values[:, 1:]
r_values = pd.read_csv('ALL_SA_OAT_r_values.csv').values[:, 1:]
data = pd.read_csv('ALL_SA_OAT_ranges.csv').values
ranges_matrix = data[:, 1:]
ranges_matrix_normalised = ranges_matrix
maxima = np.nanmax(abs(ranges_matrix), axis=0)
ranges_matrix_normalised[:, np.nonzero(maxima)] = ranges_matrix[:, np.nonzero(maxima)] / maxima[np.nonzero(maxima)]
qoi_names = pd.read_csv('ALL_SA_OAT_ranges.csv').columns.values[1:]
qoi_names = np.delete(qoi_names, np.argwhere(qoi_names=='j_point_mean')) # I don't trust the evaulation of the J point...
# qoi_names = qoi_names[0:4]
qoi = {}
qoi_ticks = []
for i in range(len(qoi_names)):
    qoi[qoi_names[i]] = len(qoi_names) - i
    qoi_ticks.append(len(qoi_names)-i)
param_names = data[:, 0]
# param_names = data[5:-3, 0] # Only cross bridge cycline, mechanics, and haemo parameters.
param = {}
param_ticks = []
for i in range(len(param_names)):
    param[param_names[i]] = len(param_names) - i
    param_ticks.append(len(param_names)-i)

# #######################################################################################################################
# All together plot
ax = plt.figure(figsize=[8, 10]).gca()
for param_i in range(len(param_names)):
    for qoi_i in range(len(qoi_names)):
        if (slopes[param_i, qoi_i] < 0.0) & (p_values[param_i, qoi_i] < 0.05): # Showing statisticall significant regression only
            c = 'b'
        elif (slopes[param_i, qoi_i] > 0.0) & (p_values[param_i, qoi_i] < 0.05):
            c = 'r'
        else:
            continue
        w = ranges_matrix_normalised[param_i, qoi_i]*1.5
        if abs(r_values[param_i, qoi_i]) > 0.6: # Showing linear correlations only
            plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                     alpha=abs(r_values[param_i, qoi_i]), color=c, linewidth=w)

ax.set_xticklabels([])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_yticks(param_ticks)
ax.set_yticklabels(param_names)
ax.set_xlim([0, 1])
ax.set_ylim([0, max(len(qoi_names), len(param_names))])
secax = ax.twinx()
secax.set_ylim(0, max(len(qoi_names), len(param_names)))
secax.set_yticks(qoi_ticks)
secax.set_yticklabels(qoi_names)
plt.tight_layout()
plt.savefig('ALL_OAT_SA_QOIS.png')
plt.show()

######################################################################################################################
# Highlight separate groups of QoIs
group = ['qrs_dur_mean', 'qrs_peak_v3', 't_dur_mean', 'qt_dur_mean', 't_pe_mean', 'jt_dur_mean', 't_peak_mean']
# group = ['qrs_dur_mean', 't_dur_mean', 'qt_dur_mean', 't_pe_mean', 'jt_dur_mean']
ax = plt.figure(figsize=[8,10]).gca()
for param_i in range(len(param_names)):
    for qoi_i in range(len(qoi_names)):
        if (slopes[param_i, qoi_i] < 0.0) & (
                p_values[param_i, qoi_i] < 0.05):  # Showing statisticall significant regression only
            c = 'b'
        elif (slopes[param_i, qoi_i] > 0.0) & (p_values[param_i, qoi_i] < 0.05):
            c = 'r'
        else:
            continue
        # if ranges_matrix_normalised[param_i, qoi_i] > 0.1:
        w = ranges_matrix_normalised[param_i, qoi_i]*1.5
        # else:
        #     w = 0
        if abs(r_values[param_i, qoi_i]) > 0.6:
            if qoi_names[qoi_i] in group:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=abs(r_values[param_i, qoi_i]), color=c, linewidth=w)
            else:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=0.1, color='k', linewidth=w)
ax.set_xticklabels([])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_yticks(param_ticks)
ax.set_yticklabels(param_names)
ax.set_xlim([0, 1])
ax.set_ylim([0, max(len(qoi_names), len(param_names))])
secax = ax.twinx()
secax.set_ylim(0, max(len(qoi_names), len(param_names)))
secax.set_yticks(qoi_ticks)
secax.set_yticklabels(qoi_names)
plt.tight_layout()
plt.savefig('OAT_SA_ECG_QOIS.png')
plt.show()


# Highlight separate groups of QoIs
group = ['EDVL', 'ESVL', 'PmaxL', 'LVEF', 'SVL', 'dvdt_ejection', 'dvdt_filling', 'dpdt_max']
ax = plt.figure(figsize=[8,10]).gca()
for param_i in range(len(param_names)):
    for qoi_i in range(len(qoi_names)):
        if (slopes[param_i, qoi_i] < 0.0) & (
                p_values[param_i, qoi_i] < 0.05):  # Showing statisticall significant regression only
            c = 'b'
        elif (slopes[param_i, qoi_i] > 0.0) & (p_values[param_i, qoi_i] < 0.05):
            c = 'r'
        else:
            continue
        # if ranges_matrix_normalised[param_i, qoi_i] > 0.1:
        w = ranges_matrix_normalised[param_i, qoi_i]*1.5
        # else:
        #     w = 0
        if abs(r_values[param_i, qoi_i]) > 0.6:
            if qoi_names[qoi_i] in group:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=abs(r_values[param_i, qoi_i]), color=c, linewidth=w)
            else:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=0.1, color='k', linewidth=w)
ax.set_xticklabels([])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_yticks(param_ticks)
ax.set_yticklabels(param_names)
ax.set_xlim([0, 1])
ax.set_ylim([0, max(len(qoi_names), len(param_names))])
secax = ax.twinx()
secax.set_ylim(0, max(len(qoi_names), len(param_names)))
secax.set_yticks(qoi_ticks)
secax.set_yticklabels(qoi_names)
plt.tight_layout()
plt.savefig('OAT_SA_PV_QOIS.png')
plt.show()

# Highlight separate groups of QoIs
group = ['es_ed_avpd', 'es_ed_apical_displacement', 'diff_lv_wall_thickness', 'percentage_volume_change']
ax = plt.figure(figsize=[8,10]).gca()
for param_i in range(len(param_names)):
    for qoi_i in range(len(qoi_names)):
        if (slopes[param_i, qoi_i] < 0.0) & (
                p_values[param_i, qoi_i] < 0.05):  # Showing statisticall significant regression only
            c = 'b'
        elif (slopes[param_i, qoi_i] > 0.0) & (p_values[param_i, qoi_i] < 0.05):
            c = 'r'
        else:
            continue
        # if ranges_matrix_normalised[param_i, qoi_i] > 0.1:
        w = ranges_matrix_normalised[param_i, qoi_i]*1.5
        # else:
        #     w = 0
        if abs(r_values[param_i, qoi_i]) > 0.6:
            if qoi_names[qoi_i] in group:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=abs(r_values[param_i, qoi_i]), color=c, linewidth=w)
            else:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=0.1, color='k', linewidth=w)
ax.set_xticklabels([])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_yticks(param_ticks)
ax.set_yticklabels(param_names)
ax.set_xlim([0, 1])
ax.set_ylim([0, max(len(qoi_names), len(param_names))])
secax = ax.twinx()
secax.set_ylim(0, max(len(qoi_names), len(param_names)))
secax.set_yticks(qoi_ticks)
secax.set_yticklabels(qoi_names)
plt.tight_layout()
plt.savefig('OAT_SA_DEFORMATION_QOIS.png')
plt.show()


# Highlight separate groups of QoIs
group = ['max_mid_Ecc', 'min_mid_Ecc', 'max_mid_Err', 'min_mid_Err', 'max_four_chamber_Ell', 'min_four_chamber_Ell']
ax = plt.figure(figsize=[8,10]).gca()
for param_i in range(len(param_names)):
    for qoi_i in range(len(qoi_names)):
        if (slopes[param_i, qoi_i] < 0.0) & (
                p_values[param_i, qoi_i] < 0.05):  # Showing statisticall significant regression only
            c = 'b'
        elif (slopes[param_i, qoi_i] > 0.0) & (p_values[param_i, qoi_i] < 0.05):
            c = 'r'
        else:
            continue
        # if ranges_matrix_normalised[param_i, qoi_i] > 0.1:
        w = ranges_matrix_normalised[param_i, qoi_i]*1.5
        # else:
        #     w = 0
        if abs(r_values[param_i, qoi_i]) > 0.6:
            if qoi_names[qoi_i] in group:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=abs(r_values[param_i, qoi_i]), color=c, linewidth=w)
            else:
                plt.plot([0, 1], [param[param_names[param_i]], qoi[qoi_names[qoi_i]]],
                         alpha=0.1, color='k', linewidth=w)
ax.set_xticklabels([])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.set_yticks(param_ticks)
ax.set_yticklabels(param_names)
ax.set_xlim([0, 1])
ax.set_ylim([0, max(len(qoi_names), len(param_names))])
secax = ax.twinx()
secax.set_ylim(0, max(len(qoi_names), len(param_names)))
secax.set_yticks(qoi_ticks)
secax.set_yticklabels(qoi_names)
plt.tight_layout()
plt.savefig('OAT_SA_STRAIN_QOIS.png')
plt.show()


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