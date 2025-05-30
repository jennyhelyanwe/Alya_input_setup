import os, sys
import numpy
import matplotlib
import numpy as np
# matplotlib.use('tkagg')  # For use on CSCS daint
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.animation as animation
import time
import argparse
import pickle
from scipy.signal import lfilter, resample
from scipy import io

class ECGPV_visualisation:
    def __init__(self, CL):
        self.CL = CL
        self.beat_fig_size = [5, 5]
        self.pv_fig_size = [6, 7]

    def read_ecg_pv(self, name, dir):
        meshname = dir+name
        print(meshname)
        if os.path.exists(dir+'ecgs.pl'):
            ecgs = pickle.load(open(dir+'ecgs.pl', 'rb'))
        else:
            ecgs = self._read_ECG(meshname)
            pickle.dump(ecgs, open(dir+'ecgs.pl', 'wb'))
        if os.path.exists(dir+'pvs.pl'):
            pvs = pickle.load(open(dir+'pvs.pl', 'rb'))
        else:
            if os.path.exists(meshname + '-cardiac-cycle.sld.res'):
                pvs = self._read_PV(meshname)
                pickle.dump(pvs, open(dir+'pvs.pl', 'wb'))
        # Save as .mat file for delineation using Karlsruhe code:
        raw_leads = numpy.vstack([ecgs['I']/ecgs['max_all_leads'],ecgs['II']/ecgs['max_all_leads'],ecgs['III']/ecgs['max_all_leads'],ecgs['aVL']/ecgs['max_all_leads'],ecgs['aVR']/ecgs['max_all_leads'],ecgs['aVF']/ecgs['max_all_leads'],ecgs['V1']/ecgs['max_all_leads'],ecgs['V2']/ecgs['max_all_leads'],ecgs['V3']/ecgs['max_all_leads'],ecgs['V4']/ecgs['max_all_leads'],ecgs['V5']/ecgs['max_all_leads'],ecgs['V6']/ecgs['max_all_leads']]).T

        raw_leads_resampled = numpy.zeros((int(ecgs['t'][-1]/0.002), 12))
        # Down sample to 500 Hz
        for i in range(0, 12):
            nsample = int(ecgs['t'][-1]/0.002)
            raw_leads_resampled[:,i] = resample(raw_leads[:,i], nsample)
        io.savemat('raw_ecg_leads.mat', {'signal':raw_leads_resampled})
        io.savemat('ecgs.mat', {'ecgs': ecgs})
        if os.path.exists(meshname + '-cardiac-cycle.sld.res'):
            io.savemat('pvs.mat', {'pvs': pvs})
        else:
            pvs = []
        return ecgs, pvs

    def _read_ECG(self, meshname):
        filename = meshname + '.exm.vin'
        data = numpy.loadtxt(filename, skiprows=8)
        t = data[:,1]
        LA = data[:,3]
        RA = data[:,4]
        LL = data[:,5]
        RL = data[:,6]
        V1 = data[:,7]
        V2 = data[:,8]
        V3 = data[:,9]
        V4 = data[:,10]
        V5 = data[:,11]
        V6 = data[:,12]

        # Ealuate Wilson's central terminal
        VW = 1.0/3.0*(RA + LA + LL)

        # Evaluate simulated ECG lead traces
        V1 = V1 - VW
        V2 = V2 - VW
        V3 = V3 - VW
        V4 = V4 - VW
        V5 = V5 - VW
        V6 = V6 - VW
        I = LA - LL
        II = LL - RA
        III = LL - LA
        aVL = LA - (RA + LL)/2.0
        aVF = LL - (LA + RA)/2.0
        aVR = RA - (LA + LL)/2.0

        all_leads = numpy.concatenate((V1,V2,V3,V4,V5,V6,I,II,III))
        max_all_leads = max(abs(all_leads))
        max_limb_leads = max(abs(numpy.concatenate((I,II,III,aVR,aVL,aVF))))

        # Divide into beats
        ts, V1s = self._divide_signal(t, V1, self.CL)
        ts, V2s = self._divide_signal(t, V2, self.CL)
        ts, V3s = self._divide_signal(t, V3, self.CL)
        ts, V4s = self._divide_signal(t, V4, self.CL)
        ts, V5s = self._divide_signal(t, V5, self.CL)
        ts, V6s = self._divide_signal(t, V6, self.CL)
        ts, aVRs = self._divide_signal(t, aVR, self.CL)
        ts, aVLs = self._divide_signal(t, aVL, self.CL)
        ts, aVFs = self._divide_signal(t, aVF, self.CL)
        ts, Is = self._divide_signal(t, I, self.CL)
        ts, IIs = self._divide_signal(t, II, self.CL)
        ts, IIIs = self._divide_signal(t, III, self.CL)

        output_dict = {'t':t,'ts':ts,'V1':V1,'V2':V2,'V3':V3,'V4':V4,'V5':V5,'V6':V6,
                    'aVR':aVR,'aVL':aVL,'aVF':aVF,'I':I,'II':II,'III':III,
                    'V1s':V1s,'V2s':V2s,'V3s':V3s,'V4s':V4s,'V5s':V5s,'V6s':V6s,
                    'aVRs':aVRs,'aVLs':aVLs,'aVFs':aVFs,'Is':Is,'IIs':IIs,'IIIs':IIIs,
                    'max_all_leads':max_all_leads, 'max_limb_leads':max_limb_leads}
        return output_dict

    def _read_PV(self,meshname):
        filename = meshname + '-cardiac-cycle.sld.res'
        data = numpy.loadtxt(filename, skiprows=18)
        curtime = data[:,1]
        phasel = data[:,4]
        vl = data[:,5]
        pl = data[:,6]
        al = data[:, 7]
        phaser = data[:, 10]
        vr = data[:, 11]
        pr = data[:, 12]
        ar = data[:, 13]
        # if (len(data)>18):
        # 	pl = numpy.zeros(len(data)-17)
        # 	pr = numpy.zeros(len(data)-17)
        # 	vl= numpy.zeros(len(data)-17)
        # 	vr = numpy.zeros(len(data)-17)
        # 	curtime = numpy.zeros(len(data)-17)
        # 	phasel = numpy.zeros(len(data)-17)
        # 	phaser = numpy.zeros(len(data)-17)
        # 	vl[0] = float(data[18].split()[5])
        # 	vr[0] = float(data[18].split()[-3])
        #
        # 	for i in range(18, len(data)):
        # 		pl[i-17] = float(data[i].split()[6])/10000
        # 		pr[i-17] = float(data[i].split()[-2])/10000
        # 		vl[i-17] = float(data[i].split()[5])
        # 		vr[i-17] = float(data[i].split()[-3])
        # 		curtime[i-17] = float(data[i].split()[1])
        # 		phasel[i-17] = float(data[i].split()[4])
        # 		phaser[i-17] = float(data[i].split()[10])

        # Divide into beats
        ts, pls = self._divide_signal(curtime, pl, self.CL)
        ts, vls = self._divide_signal(curtime, vl, self.CL)
        ts, prs = self._divide_signal(curtime, pr, self.CL)
        ts, vrs = self._divide_signal(curtime, vr, self.CL)
        ts, phasels = self._divide_signal(curtime, phasel, self.CL)
        ts, phasers = self._divide_signal(curtime, phaser, self.CL)
        output_dict = {'t':curtime, 'ts':ts, 'pl':pl, 'vl':vl, 'pr':pr, 'vr':vr, 'phasel':phasel, 'phaser':phaser,
                        'pls':pls, 'vls':vls, 'prs':prs, 'vrs':vrs, 'phasels':phasels, 'phasers':phasers,
                        'plabel':'Pressure (kPa)', 'vlabel':'Volume (mL)'}
        return output_dict

    def _divide_signal(self, curtime, signal, CL):
        if (curtime.max() > CL):
            t_offsets = []
            signals = []
            idx_start = []
            idx_end = []
            for i in range(0, int(curtime.max()/CL)+1):
                idx_start.append(numpy.where(curtime>=i*CL)[0][0])
                if curtime.max() > (i+1)*CL:
                    idx_end.append(numpy.where(curtime>=(i+1)*CL)[0][0])
                else:
                    idx_end.append(len(curtime)-1)
                t_offsets.append(curtime[idx_start[i]:idx_end[i]]-curtime[idx_start[i]])
                signals.append(signal[idx_start[i]:idx_end[i]])
        else:
            t_offsets = [curtime]
            signals = [signal]
        return t_offsets, signals

    def plot_ecgpv_live(self, ecgs, pvs, title, show, ecgs2=[], pvs2=[]):
        print('Plotting ECG PV live')
        matplotlib.rcParams.update({'font.size':'11'})
        matplotlib.rcParams.update({'text.color':'black'})
        matplotlib.rcParams.update({'lines.linewidth':'1'})
        fig = plt.figure(tight_layout=True, figsize=[15,7])
        fig.suptitle(title)
        gs = GridSpec(3,6)
        axs = []
        axs.append(fig.add_subplot(gs[:,0]))
        axs.append(fig.add_subplot(gs[0,1]))
        axs.append(fig.add_subplot(gs[1,1]))
        axs.append(fig.add_subplot(gs[2,1]))
        axs.append(fig.add_subplot(gs[0,2]))
        axs.append(fig.add_subplot(gs[1,2]))
        axs.append(fig.add_subplot(gs[2,2]))
        axs.append(fig.add_subplot(gs[0,3]))
        axs.append(fig.add_subplot(gs[1,3]))
        axs.append(fig.add_subplot(gs[2,3]))
        axs.append(fig.add_subplot(gs[0,4]))
        axs.append(fig.add_subplot(gs[1,4]))
        axs.append(fig.add_subplot(gs[2,4]))
        axs.append(fig.add_subplot(gs[0,5]))
        axs.append(fig.add_subplot(gs[1,5]))
        axs.append(fig.add_subplot(gs[2,5]))

        def animate_single(i):
            # Plot PV
            axs[0].clear()
            axs[0].plot(pvs['vl'], pvs['pl'], pvs['vr'], pvs['pr'])
            axs[1].clear()
            axs[1].plot(pvs['t'], pvs['pl'], pvs['t'], pvs['pr'])
            axs[2].clear()
            axs[2].plot(pvs['t'], pvs['vl'], pvs['t'], pvs['vr'])
            axs[4].clear()
            axs[3].plot(pvs['t'], pvs['phasel'], pvs['t'], pvs['phaser'])

            # Plot ECGs:
            axs[4].clear()
            axs[4].plot(ecgs['t'], ecgs['I']/ecgs['max_all_leads'])
            axs[5].clear()
            axs[5].plot(ecgs['t'], ecgs['II']/ecgs['max_all_leads'])
            axs[6].clear()
            axs[6].plot(ecgs['t'], ecgs['III']/ecgs['max_all_leads'])
            axs[7].clear()
            axs[7].plot(ecgs['t'], ecgs['aVR']/ecgs['max_all_leads'])
            axs[8].clear()
            axs[8].plot(ecgs['t'], ecgs['aVL']/ecgs['max_all_leads'])
            axs[9].clear()
            axs[9].plot(ecgs['t'], ecgs['aVF']/ecgs['max_all_leads'])
            axs[10].clear()
            axs[10].plot(ecgs['t'], ecgs['V1']/ecgs['max_all_leads'])
            axs[11].clear()
            axs[11].plot(ecgs['t'], ecgs['V2']/ecgs['max_all_leads'])
            axs[12].clear()
            axs[12].plot(ecgs['t'], ecgs['V3']/ecgs['max_all_leads'])
            axs[13].clear()
            axs[13].plot(ecgs['t'], ecgs['V4']/ecgs['max_all_leads'])
            axs[14].clear()
            axs[14].plot(ecgs['t'], ecgs['V5']/ecgs['max_all_leads'])
            axs[15].clear()
            axs[15].plot(ecgs['t'], ecgs['V6']/ecgs['max_all_leads'])

            axs[0].set_xlabel('Volume (mL)')
            axs[1].set_xlabel('Time (ms)')
            axs[2].set_xlabel('Time (ms)')
            axs[3].set_xlabel('Time (ms)')
            axs[0].set_ylabel('Pressure (kPa)')
            axs[1].set_ylabel('Pressure (kPa)')
            axs[2].set_ylabel('Volume (mL)')
            axs[3].set_ylabel('Phase')
            axs[3].set_ylim([0,5])
            axs[1].set_title('Pressure transient')
            axs[2].set_title('Volume transient')
            axs[3].set_title('Cardiac cycle')
            axs[4].set_xlabel('Time (s)')
            axs[5].set_xlabel('Time (s)')
            axs[6].set_xlabel('Time (s)')
            axs[7].set_xlabel('Time (s)')
            axs[8].set_xlabel('Time (s)')
            axs[9].set_xlabel('Time (s)')
            axs[10].set_xlabel('Time (s)')
            axs[11].set_xlabel('Time (s)')
            axs[12].set_xlabel('Time (s)')
            axs[13].set_xlabel('Time (s)')
            axs[14].set_xlabel('Time (s)')
            axs[15].set_xlabel('Time (s)')
            axs[4].set_ylabel('Normalised ECG')
            axs[5].set_ylabel('Normalised ECG')
            axs[6].set_ylabel('Normalised ECG')
            axs[7].set_ylabel('Normalised ECG')
            axs[8].set_ylabel('Normalised ECG')
            axs[9].set_ylabel('Normalised ECG')
            axs[10].set_ylabel('Normalised ECG')
            axs[11].set_ylabel('Normalised ECG')
            axs[12].set_ylabel('Normalised ECG')
            axs[13].set_ylabel('Normalised ECG')
            axs[14].set_ylabel('Normalised ECG')
            axs[15].set_ylabel('Normalised ECG')
            axs[4].set_ylim(-1,1)
            axs[5].set_ylim(-1,1)
            axs[6].set_ylim(-1,1)
            axs[7].set_ylim(-1,1)
            axs[8].set_ylim(-1,1)
            axs[9].set_ylim(-1,1)
            axs[10].set_ylim(-1,1)
            axs[11].set_ylim(-1,1)
            axs[12].set_ylim(-1,1)
            axs[13].set_ylim(-1,1)
            axs[14].set_ylim(-1,1)
            axs[15].set_ylim(-1,1)
            axs[4].set_title('I')
            axs[5].set_title('II')
            axs[6].set_title('III')
            axs[7].set_title('aVR')
            axs[8].set_title('aVL')
            axs[9].set_title('aVF')
            axs[10].set_title('V1')
            axs[11].set_title('V2')
            axs[12].set_title('V3')
            axs[13].set_title('V4')
            axs[14].set_title('V5')
            axs[15].set_title('V6')

        def animate_double(i):
            # Plot PV
            axs[0].clear()
            axs[0].plot(pvs['vl'], pvs['pl'], pvs['vr'], pvs['pr'], pvs2['vl'],pvs2['pl'],'--', pvs2['vr'], pvs2['pr'], '--')
            axs[1].clear()
            axs[1].plot(pvs['t'], pvs['pl'], pvs['t'], pvs['pr'],pvs2['t'], pvs2['pl'],'--', pvs2['t'], pvs2['pr'], '--')
            axs[2].clear()
            axs[2].plot(pvs['t'], pvs['vl'], pvs['t'], pvs['vr'],pvs2['t'], pvs2['vl'], '--', pvs2['t'],  pvs2['vr'], '--')
            axs[4].clear()
            axs[3].plot(pvs['t'], pvs['phasel'], pvs['t'], pvs['phaser'],pvs2['t'], pvs2['phasel'], '--', pvs2['t'], pvs2['phaser'], '--')
            axs[0].set_xlabel('Volume (mL)')
            axs[1].set_xlabel('Time (ms)')
            axs[2].set_xlabel('Time (ms)')
            axs[3].set_xlabel('Time (ms)')
            axs[0].set_ylabel('Pressure (kPa)')
            axs[1].set_ylabel('Pressure (kPa)')
            axs[2].set_ylabel('Volume (mL)')
            axs[3].set_ylabel('Phase')
            axs[3].set_ylim([0,5])

            axs[1].set_title('Pressure transient')
            axs[2].set_title('Volume transient')
            axs[3].set_title('Cardiac cycle')

            # Plot ECGs:
            max_all_leads = max([ecgs['max_all_leads'], ecgs2['max_all_leads']])
            axs[4].clear()
            axs[4].plot(ecgs['t'], ecgs['I']/max_all_leads, ecgs2['t'], ecgs2['I']/max_all_leads, '--')
            axs[5].clear()
            axs[5].plot(ecgs['t'], ecgs['II']/max_all_leads, ecgs2['t'],ecgs2['II']/max_all_leads, '--')
            axs[6].clear()
            axs[6].plot(ecgs['t'], ecgs['III']/max_all_leads, ecgs2['t'], ecgs2['III']/max_all_leads, '--')
            axs[7].clear()
            axs[7].plot(ecgs['t'], ecgs['aVR']/max_all_leads, ecgs2['t'], ecgs2['aVR']/max_all_leads, '--')
            axs[8].clear()
            axs[8].plot(ecgs['t'], ecgs['aVL']/max_all_leads, ecgs2['t'], ecgs2['aVL']/max_all_leads, '--')
            axs[9].clear()
            axs[9].plot(ecgs['t'], ecgs['aVF']/max_all_leads, ecgs2['t'], ecgs2['aVF']/max_all_leads, '--')
            axs[10].clear()
            axs[10].plot(ecgs['t'], ecgs['V1']/max_all_leads, ecgs2['t'], ecgs2['V1']/max_all_leads, '--')
            axs[11].clear()
            axs[11].plot(ecgs['t'], ecgs['V2']/max_all_leads, ecgs2['t'], ecgs2['V2']/max_all_leads, '--')
            axs[12].clear()
            axs[12].plot(ecgs['t'], ecgs['V3']/max_all_leads, ecgs2['t'], ecgs2['V3']/max_all_leads, '--')
            axs[13].clear()
            axs[13].plot(ecgs['t'], ecgs['V4']/max_all_leads, ecgs2['t'], ecgs2['V4']/max_all_leads, '--')
            axs[14].clear()
            axs[14].plot(ecgs['t'], ecgs['V5']/max_all_leads, ecgs2['t'], ecgs2['V5']/max_all_leads, '--')
            axs[15].clear()
            axs[15].plot(ecgs['t'], ecgs['V6']/max_all_leads, ecgs2['t'], ecgs2['V6']/max_all_leads, '--')
            axs[4].set_xlabel('Time (s)')
            axs[5].set_xlabel('Time (s)')
            axs[6].set_xlabel('Time (s)')
            axs[7].set_xlabel('Time (s)')
            axs[8].set_xlabel('Time (s)')
            axs[9].set_xlabel('Time (s)')
            axs[10].set_xlabel('Time (s)')
            axs[11].set_xlabel('Time (s)')
            axs[12].set_xlabel('Time (s)')
            axs[13].set_xlabel('Time (s)')
            axs[14].set_xlabel('Time (s)')
            axs[15].set_xlabel('Time (s)')
            axs[4].set_ylabel('Normalised ECG')
            axs[5].set_ylabel('Normalised ECG')
            axs[6].set_ylabel('Normalised ECG')
            axs[7].set_ylabel('Normalised ECG')
            axs[8].set_ylabel('Normalised ECG')
            axs[9].set_ylabel('Normalised ECG')
            axs[10].set_ylabel('Normalised ECG')
            axs[11].set_ylabel('Normalised ECG')
            axs[12].set_ylabel('Normalised ECG')
            axs[13].set_ylabel('Normalised ECG')
            axs[14].set_ylabel('Normalised ECG')
            axs[15].set_ylabel('Normalised ECG')
            axs[4].set_ylim(-1,1)
            axs[5].set_ylim(-1,1)
            axs[6].set_ylim(-1,1)
            axs[7].set_ylim(-1,1)
            axs[8].set_ylim(-1,1)
            axs[9].set_ylim(-1,1)
            axs[10].set_ylim(-1,1)
            axs[11].set_ylim(-1,1)
            axs[12].set_ylim(-1,1)
            axs[13].set_ylim(-1,1)
            axs[14].set_ylim(-1,1)
            axs[15].set_ylim(-1,1)
            axs[4].set_title('I')
            axs[5].set_title('II')
            axs[6].set_title('III')
            axs[7].set_title('aVR')
            axs[8].set_title('aVL')
            axs[9].set_title('aVF')
            axs[10].set_title('V1')
            axs[11].set_title('V2')
            axs[12].set_title('V3')
            axs[13].set_title('V4')
            axs[14].set_title('V5')
            axs[15].set_title('V6')

        # if ecgs2:
        #     ani = animation.FuncAnimation(fig, animate_double, fargs=( axs, ecgs, pvs, ecgs2, pvs2), interval=1000)
        # else:
        ani = animation.FuncAnimation(fig, animate_single, interval=1000)
        if show:
            plt.show(block=True)

    def _set_ecg_ticks(self,ax, t_end, CL):
        t_start = 0
        t_end = numpy.ceil(t_end/CL) * CL
        t_end = int(numpy.ceil(t_end*10.0/2.0)*2)/10.0
        #t_end = numpy.round(t_end*10.0)/10.0
        #t_end = numpy.round(t_end)
        #t_end = self.CL
        minor_ticks = numpy.arange(0, t_end+0.04, 0.04)
        major_ticks = numpy.arange(0, t_end+0.2, 0.2)
        #minor_ticks = numpy.linspace(0, t_end, int((t_end-t_start)/0.04) + 1)
        #major_ticks = numpy.linspace(0, t_end, int((t_end-t_start)/0.2) + 1)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_xticks(major_ticks)
        minor_ticks = numpy.linspace(-1, 1, 21)
        major_ticks = numpy.linspace(-1, 1, 5)
        ax.set_yticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.grid(which="minor", color='r', linestyle='-', linewidth=1, alpha=0.5)
        ax.grid(which="major", color='r', linestyle='-', linewidth=2, alpha=0.5)

    def _set_full_ecg_ticks(self,ax, t_end, v_max, v_min, CL):
        t_start = 0
        t_end = numpy.ceil(t_end/CL) * CL
        t_end = int(numpy.ceil(t_end*10.0/2.0)*2)/10.0
        minor_ticks = numpy.arange(0, t_end+0.04, 0.04)
        major_ticks = numpy.arange(0, t_end+0.2, 0.2)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_xticks(major_ticks)
        v_max = numpy.ceil(v_max)
        minor_ticks = numpy.arange(v_min, v_max+0.1, 0.1)
        major_ticks = numpy.arange(v_min, v_max+0.5, 0.5)
        ax.set_yticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.grid(which="minor", color='r', linestyle='-', linewidth=1, alpha=0.3)
        ax.grid(which="major", color='r', linestyle='-', linewidth=2, alpha=0.3)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1line.set_visible(False)
            tic.tick2line.set_visible(False)
            tic.label1.set_visible(False)
            tic.label2.set_visible(False)
        for tic in ax.xaxis.get_minor_ticks():
            tic.tick1line.set_visible(False)
            tic.tick2line.set_visible(False)
            tic.label1.set_visible(False)
            tic.label2.set_visible(False)
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1line.set_visible(False)
            tic.tick2line.set_visible(False)
            tic.label1.set_visible(False)
            tic.label2.set_visible(False)
        for tic in ax.yaxis.get_minor_ticks():
            tic.tick1line.set_visible(False)
            tic.tick2line.set_visible(False)
            tic.label1.set_visible(False)
            tic.label2.set_visible(False)

    def _set_pv_ticks(self,ax):
        x_ticks = numpy.arange(50, 250, 50)
        y_ticks = numpy.arange(0, 140, 20)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        #ax.grid(which="major", color='k', linestyle='--', linewidth=1, alpha=0.5)

    def plot_ecg_lead(self, ecgs, lead_name, filename,show, beat=0, ecg2=[]):
        print('Plotting ECG lead '+str(lead_name))
        matplotlib.rcParams.update({'font.size':'24'})
        matplotlib.rcParams.update({'text.color':'black'})
        matplotlib.rcParams.update({'lines.linewidth':'3'})
        if ecg2:
            if beat > 0:
                max_all_leads = max([ecgs['max_all_leads'], ecgs2['max_all_leads']])
                self._plot_double(self.beat_fig_size,
                            ecgs['ts'][beat-1],ecgs[lead_name+'s'][beat-1]/max_all_leads,'-',
                            ecgs2['ts'][beat-1], ecgs[lead_name+'s'][beat-1]/max_all_leads,'--',
                            'Time (s)', lead_name, filename+'_beat'+str(beat)+'_comparison.png',
                            show,ecg_grid=True)
            else:
                figsize=[self.beat_fig_size[0]*(int(ecgs['t'][-1]/1.0)+1), self.beat_fig_size[1]]
                self._plot_double(figsize, ecgs['t'], ecgs[lead_name]/max_all_leads,'-',
                                ecgs2['t'], ecgs2[lead_name]/max_all_leads,'--'
                                'Time (s)', lead_name, filename+'_full_comparison.png',show, ecg_grid=True)
        else:
            if beat > 0:
                self._plot_single(self.beat_fig_size, ecgs['ts'][beat-1],
                            ecgs[lead_name+'s'][beat-1]/ecgs['max_all_leads'],
                            'Time (s)', lead_name, filename+'_beat'+str(beat)+'.png',
                            show, ecg_grid=True)
            else:
                figsize=[self.beat_fig_size[0]*(int(ecgs['t'][-1]/1.0)+1), self.beat_fig_size[1]]
                self._plot_single(figsize, ecgs['t'], ecgs[lead_name]/ecgs['max_all_leads'],
                                'Time (s)', lead_name, filename+'_full.png', show,ecg_grid=True)

    def plot_ecg_all_leads(self, ecgs, show, beat=0, ecgs2=[]):
        print('Plotting all ECG leads')
        matplotlib.rcParams.update({'font.size':'24'})
        matplotlib.rcParams.update({'text.color':'black'})
        matplotlib.rcParams.update({'lines.linewidth':'1'})

        if ecgs2:
            print('Full ECG comparison not yet implemented...')
        else:
            if beat > 0:
                # t_padding
                t_padding = 0.2
                v_padding = 0.5

                # Generate calibrating step function
                t_res = ecgs['t'][1] - ecgs['t'][0]
                t_calib = numpy.arange(0, 0.2+t_res, t_res)
                v_calib = numpy.zeros(numpy.shape(t_calib))
                for i in range(0, len(t_calib)):
                    if t_calib[i] < 0.04:
                        v_calib[i] = 0.0
                    elif t_calib[i] > t_calib[-1]-0.04:
                        v_calib[i] = 0.0
                    else:
                        v_calib[i] = 1.0

                # Concatenate the leads together to plot in a single figure
                t_end = ecgs['ts'][beat-1][-1]
                full_t = numpy.concatenate([t_calib, ecgs['ts'][beat-1]+t_calib[-1], ecgs['ts'][beat-1]+t_end+t_calib[-1], ecgs['ts'][beat-1]+2*t_end+t_calib[-1], ecgs['ts'][beat-1]+3*t_end+t_calib[-1]]) + t_padding
                scale = 1
                fig_size = [numpy.ceil(full_t[-1]+t_padding)*0.5/0.2*scale,6*scale]
                fig = plt.figure(tight_layout=True, figsize=fig_size)
                gs = GridSpec(1,1)
                ax = fig.add_subplot(gs[0,0])

                # Bottom row: III, aVF, V3, V6
                bottom_V = numpy.concatenate([v_calib, ecgs['IIIs'][beat-1]/ecgs['max_all_leads'], ecgs['aVFs'][beat-1]/ecgs['max_all_leads'], ecgs['V3s'][beat-1]/ecgs['max_all_leads'], ecgs['V6s'][beat-1]/ecgs['max_all_leads']])
                ax.plot(full_t, bottom_V, 'k')

                # Middle row: II, aVL, V2, V5
                offset = 2
                midrow_V = numpy.concatenate([v_calib+offset, ecgs['IIs'][beat-1]/ecgs['max_all_leads']+offset, ecgs['aVLs'][beat-1]/ecgs['max_all_leads']+offset, ecgs['V2s'][beat-1]/ecgs['max_all_leads']+offset, ecgs['V5s'][beat-1]/ecgs['max_all_leads']+offset])
                ax.plot(full_t, midrow_V, 'k')

                # Top row: I, aVR, V1, V4
                toprow_V = numpy.concatenate([v_calib+offset*2, ecgs['Is'][beat-1]/ecgs['max_all_leads']+offset*2, ecgs['aVRs'][beat-1]/ecgs['max_all_leads']+offset*2, ecgs['V1s'][beat-1]/ecgs['max_all_leads']+offset*2, ecgs['V4s'][beat-1]/ecgs['max_all_leads']+offset*2])
                ax.plot(full_t, toprow_V, 'k')

                # ax.set_xlim([0, t_end*4])
                # ax.set_ylim([-1, offset*2+1])
                #
                # self._set_full_ecg_ticks(ax, full_t[-1]+t_padding, offset*2+1, self.CL)

                # Add ECG red grid
                self._set_full_ecg_ticks(ax, full_t[-1]+t_padding*2, offset*2 + 1 + v_padding, -1-v_padding, self.CL)
                ax.set_xlim([0, full_t[-1]+t_padding*2])
                ax.set_ylim([-1-v_padding, offset*2 + 1 +v_padding])

                # Label the leads
                label_fontsize = 'large'
                label_t_offset = t_calib[-1] + t_padding
                label_v_offset = 0.1
                plt.text(0+label_t_offset, -1+label_v_offset, 'III', fontsize=label_fontsize)
                plt.text(t_end+label_t_offset, -1+label_v_offset, 'aVF', fontsize=label_fontsize)
                plt.text(2*t_end+label_t_offset, -1+label_v_offset, 'V3', fontsize=label_fontsize)
                plt.text(3*t_end+label_t_offset, -1+label_v_offset, 'V6', fontsize=label_fontsize)
                plt.text(0+label_t_offset, -1+offset+label_v_offset, 'II', fontsize=label_fontsize)
                plt.text(t_end+label_t_offset, -1+offset+label_v_offset, 'aVL', fontsize=label_fontsize)
                plt.text(2*t_end+label_t_offset, -1+offset+label_v_offset, 'V2', fontsize=label_fontsize)
                plt.text(3*t_end+label_t_offset, -1+offset+label_v_offset, 'V5', fontsize=label_fontsize)
                plt.text(0+label_t_offset, -1+offset*2+label_v_offset, 'I', fontsize=label_fontsize)
                plt.text(t_end+label_t_offset, -1+offset*2+label_v_offset, 'aVR', fontsize=label_fontsize)
                plt.text(2*t_end+label_t_offset, -1+offset*2+label_v_offset, 'V1', fontsize=label_fontsize)
                plt.text(3*t_end+label_t_offset, -1+offset*2+label_v_offset, 'V4', fontsize=label_fontsize)
                plt.savefig('full_ecg_beat_'+str(beat)+'.png')
                if show:
                    plt.show()
            else:
                # Padding
                t_padding = 0.2
                v_padding = 0.5

                # Generate calibrating step function
                t_res = ecgs['t'][1] - ecgs['t'][0]
                t_calib = numpy.arange(0, 0.2+t_res, t_res)
                v_calib = numpy.zeros(numpy.shape(t_calib))
                for i in range(0, len(t_calib)):
                    if t_calib[i] < 0.04:
                        v_calib[i] = 0.0
                    elif t_calib[i] > t_calib[-1]-0.04:
                        v_calib[i] = 0.0
                    else:
                        v_calib[i] = 1.0

                # Concatenate the leads together to plot in a single figure
                t_end = ecgs['t'][-1]
                full_t = numpy.concatenate([t_calib, ecgs['t']+t_calib[-1], ecgs['t']+t_end+t_calib[-1], ecgs['t']+2*t_end+t_calib[-1], ecgs['t']+3*t_end+t_calib[-1]]) + t_padding
                scale = 0.8
                fig_size = [numpy.ceil(full_t[-1]+t_padding)*0.5/0.2*scale,6*scale]
                fig = plt.figure(tight_layout=True, figsize=fig_size)
                gs = GridSpec(1,1)
                ax = fig.add_subplot(gs[0,0])
                # Bottom row: III, aVF, V3, V6
                bottom_V = numpy.concatenate([v_calib, ecgs['III']/ecgs['max_all_leads'], ecgs['aVF']/ecgs['max_all_leads'], ecgs['V3']/ecgs['max_all_leads'], ecgs['V6']/ecgs['max_all_leads']])
                ax.plot(full_t, bottom_V, 'k')

                # Middle row: II, aVL, V2, V5
                offset = 2
                midrow_V = numpy.concatenate([v_calib+offset, ecgs['II']/ecgs['max_all_leads']+offset, ecgs['aVL']/ecgs['max_all_leads']+offset, ecgs['V2']/ecgs['max_all_leads']+offset, ecgs['V5']/ecgs['max_all_leads']+offset])
                ax.plot(full_t, midrow_V, 'k')

                # Top row: I, aVR, V1, V4
                toprow_V = numpy.concatenate([v_calib+offset*2, ecgs['I']/ecgs['max_all_leads']+offset*2, ecgs['aVR']/ecgs['max_all_leads']+offset*2, ecgs['V1']/ecgs['max_all_leads']+offset*2, ecgs['V4']/ecgs['max_all_leads']+offset*2])
                ax.plot(full_t, toprow_V, 'k')

                # Add ECG red grid
                self._set_full_ecg_ticks(ax, full_t[-1]+t_padding*2, offset*2 + 1 + v_padding, -1-v_padding, self.CL)
                ax.set_xlim([0, full_t[-1]+t_padding*2])
                ax.set_ylim([-1-v_padding, offset*2 + 1 +v_padding])

                # Label the leads
                label_fontsize = 'xx-small'
                label_t_offset = t_calib[-1] + t_padding
                label_v_offset = 0.1
                plt.text(0+label_t_offset, -1+label_v_offset, 'III', fontsize=label_fontsize)
                plt.text(t_end+label_t_offset, -1+label_v_offset, 'aVF', fontsize=label_fontsize)
                plt.text(2*t_end+label_t_offset, -1+label_v_offset, 'V3', fontsize=label_fontsize)
                plt.text(3*t_end+label_t_offset, -1+label_v_offset, 'V6', fontsize=label_fontsize)
                plt.text(0+label_t_offset, -1+offset+label_v_offset, 'II', fontsize=label_fontsize)
                plt.text(t_end+label_t_offset, -1+offset+label_v_offset, 'aVL', fontsize=label_fontsize)
                plt.text(2*t_end+label_t_offset, -1+offset+label_v_offset, 'V2', fontsize=label_fontsize)
                plt.text(3*t_end+label_t_offset, -1+offset+label_v_offset, 'V5', fontsize=label_fontsize)
                plt.text(0+label_t_offset, -1+offset*2+label_v_offset, 'I', fontsize=label_fontsize)
                plt.text(t_end+label_t_offset, -1+offset*2+label_v_offset, 'aVR', fontsize=label_fontsize)
                plt.text(2*t_end+label_t_offset, -1+offset*2+label_v_offset, 'V1', fontsize=label_fontsize)
                plt.text(3*t_end+label_t_offset, -1+offset*2+label_v_offset, 'V4', fontsize=label_fontsize)
                plt.savefig('full_ecg_all_beats.png')
                if show:
                    plt.show()



    def plot_pv_signal(self, pvs, signal_name, filename, show, beat=0, pvs2=[]):
        matplotlib.rcParams.update({'font.size':'28'})
        matplotlib.rcParams.update({'text.color':'black'})
        matplotlib.rcParams.update({'lines.linewidth':'3'})
        if pvs2:
            if beat > 0:
                if (signal_name == 'p') | (signal_name == 'v'):
                    self._plot_quadruple(self.pv_fig_size,
                                pvs['ts'][beat-1], pvs[signal_name+'ls'][beat-1],'b--',
                                pvs['ts'][beat-1], pvs[signal_name+'rs'][beat-1],'g--',
                                pvs2['ts'][beat-1], pvs2[signal_name+'ls'][beat-1],'b-',
                                pvs2['ts'][beat-1], pvs2[signal_name+'rs'][beat-1],'g-',
                                'Time (s)',pvs[signal_name+'label'], filename+'_beat'+str(beat)+'_compare.png', show)
                elif (signal_name == 'pv'):
                    self._plot_quadruple(self.pv_fig_size,
                                pvs['vls'][beat-1],pvs['pls'][beat-1]*7.5,'b--',
                                pvs['vrs'][beat-1],pvs['prs'][beat-1]*7.5,'g--',
                                pvs2['vls'][beat-1],pvs2['pls'][beat-1]*7.5,'b-',
                                pvs2['vrs'][beat-1],pvs2['prs'][beat-1]*7.5,'g-',
                                'Volume (mL)', 'Pressure (mmHg)',
                                filename+'_beat'+str(beat)+'_compare.png', show)
            else:
                if (signal_name == 'p') | (signal_name == 'v'):
                    figsize=[self.pv_fig_size[0]*(int(pvs['t'][-1]/1.0)+1), self.pv_fig_size[1]]
                    self._plot_quadruple(figsize,
                                    pvs['t'], pvs[signal_name+'l'],'b--',
                                    pvs['t'], pvs[signal_name+'r'],'g--',
                                    pvs2['t'], pvs2[signal_name+'l'],'b-',
                                    pvs2['t'], pvs2[signal_name+'r'],'g-','Time (s)',
                                    pvs[signal_name+'label'], filename+'_full_compare.png', show)
                elif (signal_name == 'pv'):
                    self._plot_quadruple(self.pv_fig_size,
                                    pvs['vl'], pvs['pl']*7.5,'b--',
                                    pvs['vr'], pvs['pr']*7.5, 'g--',
                                    pvs2['vl'], pvs2['pl']*7.5,'b-',
                                    pvs2['vr'], pvs2['pr']*7.5, 'g-',
                                    'Volume (mL)','Pressure (mmHg)', filename+'_full_compare.png', show)
        else:
            if beat > 0:
                if (signal_name == 'p') | (signal_name == 'v'):
                    self._plot_double(self.pv_fig_size,
                                pvs['ts'][beat-1],pvs[signal_name+'ls'][beat-1],'b-',
                                pvs['ts'][beat-1],pvs[signal_name+'rs'][beat-1],'g-',
                                'Time (s)',pvs[signal_name+'label'], filename+'_beat'+str(beat)+'.png',show)
                elif (signal_name == 'pv'):
                    self._plot_double(self.pv_fig_size,
                                pvs['vls'][beat-1],pvs['pls'][beat-1], 'b-',
                                pvs['vrs'][beat-1],pvs['prs'][beat-1], 'g-',
                                'Volume (mL)', 'Pressure (kPa)',
                                filename+'_beat'+str(beat)+'.png',show)
            else:
                if (signal_name == 'p') | (signal_name == 'v'):
                    figsize=[self.pv_fig_size[0]*(int(pvs['t'][-1]/1.0)+1), self.pv_fig_size[1]]
                    self._plot_double(figsize, pvs['t'], pvs[signal_name+'l'],'b-',
                                    pvs['t'], pvs[signal_name+'r'],'g-',
                                    'Time (s)',pvs[signal_name+'label'], filename+'_full.png',show)
                elif (signal_name == 'pv'):
                    self._plot_double(self.pv_fig_size, pvs['vl'], pvs['pl'],'b-',
                                    pvs['vr'], pvs['pr'], 'r-',
                                    'Volume (mL)','Pressure (kPa)', filename+'_full.png',show)

    def _plot_single(self, fig_size, x, y, xlabel, ylabel, filename, show, ecg_grid=False):
        fig = plt.figure(tight_layout=True, figsize=fig_size)
        gs = GridSpec(1,1)
        ax = fig.add_subplot(gs[0,0])
        ax.plot(x, y)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if ecg_grid:
            print('setting ecg grids')
            self._set_ecg_ticks(ax, x[-1], self.CL)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(left=False, right=False, top=False, bottom=False, labelleft=False, labelbottom=False)
            ax.set_xlim([0, numpy.ceil(x[-1]/self.CL)*CL])
        if show:
            plt.show(block=True)
        plt.savefig(filename)

    def _plot_double(self, fig_size, x, y, linestyle, x1, y1, linestyle1, xlabel, ylabel, filename, show, ecg_grid=False):
        fig = plt.figure(tight_layout=True, figsize=fig_size)
        gs = GridSpec(1,1)
        ax = fig.add_subplot(gs[0,0])
        ax.plot(x, y, linestyle, x1, y1, linestyle1)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if ecg_grid:
            self._set_ecg_ticks(ax, x[-1],self.CL)
        if show:
            plt.show(block=True)
        plt.savefig(filename)

    def _plot_quadruple(self, fig_size, x, y, linestyle,  x1, y1,linestyle1, x2, y2, linestyle2, x3, y3,linestyle3, xlabel, ylabel, filename, show, ecg_grid=False):
        fig = plt.figure(tight_layout=True, figsize=fig_size)
        gs = GridSpec(1,1)
        ax = fig.add_subplot(gs[0,0])
        ax.plot(x, y, linestyle,x1, y1,linestyle1, x2, y2, linestyle2, x3, y3,linestyle3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_ylim([0,120])
        self._set_pv_ticks(ax)
        if ecg_grid:
            self._set_ecg_ticks(ax, x[-1], self.CL)
        if show:
            plt.show(block=True)
        plt.savefig(filename)

    def read_single_cell(self, dir, material):
        # Read membrane potentials
        filename = dir+'torord_ode.m'+str(material)+'c3.csv'
        with open(filename, 'r') as f:
            data = f.readlines()
        data = data[1:]
        epi = numpy.zeros((len(data),5))
        for i in range(0, len(epi)): epi[i,1] = float(data[i].split()[2])*1000000.0 # Calcium transient
        for i in range(0, len(epi)): epi[i,2] = float(data[i].split()[3]) # Membrane potential
        for i in range(0, len(epi)): epi[i,0] = float(data[i].split()[0])/1000.0 # real time
        for i in range(0, len(epi)): epi[i,3] = float(data[i].split()[1])/1000.0 # each beat time
        for i in range(0, len(epi)): epi[i,4] = float(data[i].split()[4]) # Ta

        filename = dir+'torord_ode.m'+str(material)+'c2.csv'
        with open(filename, 'r') as f:
            data = f.readlines()
        data = data[1:]
        mid = numpy.zeros((len(data),5))
        for i in range(0, len(mid)): mid[i,1] = float(data[i].split()[2])*1000000.0 # Calcium transient
        for i in range(0, len(mid)): mid[i,2] = float(data[i].split()[3]) # Membrane potential
        for i in range(0, len(mid)): mid[i,0] = float(data[i].split()[0])/1000.0 # real time
        for i in range(0, len(mid)): mid[i,3] = float(data[i].split()[1])/1000.0 # each beat time
        for i in range(0, len(mid)): mid[i,4] = float(data[i].split()[4]) # Ta

        filename = dir+'torord_ode.m'+str(material)+'c1.csv'
        with open(filename, 'r') as f:
            data = f.readlines()
        data = data[1:]
        endo = numpy.zeros((len(data),5))
        for i in range(0, len(endo)): endo[i,1] = float(data[i].split()[2])*1000000.0 # Calcium transient
        for i in range(0, len(endo)): endo[i,2] = float(data[i].split()[3]) # Membrane potential
        for i in range(0, len(endo)): endo[i,0] = float(data[i].split()[0])/1000.0 # real time
        for i in range(0, len(endo)): endo[i,3] = float(data[i].split()[1])/1000.0 # each beat time
        for i in range(0, len(endo)): endo[i,4] = float(data[i].split()[4]) # Ta
        return epi, mid, endo

    def plot_single_cell(self, epi, mid, endo, material,analysis, show, epi2=[], mid2=[], endo2=[]):
        # Set up subplots
        fig = plt.figure(tight_layout=True, figsize=[15,5])
        gs = GridSpec(3,7)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0])
        ax3 = fig.add_subplot(gs[2,0])
        ax4 = fig.add_subplot(gs[:,1:3])
        ax5 = fig.add_subplot(gs[:,3:5])
        ax6 = fig.add_subplot(gs[:,5:7])

        if len(epi2)>0:
            ax1.plot(epi[:,0], epi[:,2], mid[:,0], mid[:,2], endo[:,0], endo[:,2],
                    epi2[:,0], epi2[:,2],'--', mid2[:,0], mid2[:,2],'--', endo2[:,0], endo2[:,2], '--')
            ax2.plot(epi[:,0], epi[:,1],mid[:,0],mid[:,1],endo[:,0],endo[:,1],
                    epi2[:,0], epi2[:,1],'--',mid2[:,0],mid2[:,1],'--',endo2[:,0],endo2[:,1],'--')
            ax3.plot(epi[:,0], epi[:,4],mid[:,0],mid[:,4],endo[:,0],endo[:,4],
                    epi2[:,0], epi2[:,4],'--',mid2[:,0],mid2[:,4],'--',endo2[:,0],endo2[:,4],'--')
            dt = epi[1,0] - epi[0,0]
            idx = int(self.CL / dt)
            t_s_epi = epi[-idx,0]
            t_s_mid = mid[-idx,0]
            t_s_endo = endo[-idx,0]
            t_s_epi2 = epi2[-idx,0]
            t_s_mid2 = mid2[-idx,0]
            t_s_endo2 = endo2[-idx,0]
            ax4.plot(epi[-idx:,0]-t_s_epi, epi[-idx:,2], mid[-idx:,0]-t_s_mid, mid[-idx:,2], endo[-idx:,0]-t_s_endo, endo[-idx:,2],
                    epi2[-idx:,0]-t_s_epi2, epi2[-idx:,2],'--', mid2[-idx:,0]-t_s_mid2, mid2[-idx:,2],'--', endo2[-idx:,0]-t_s_endo2, endo2[-idx:,2], '--')
            ax5.plot(epi[-idx:,0]-t_s_epi, epi[-idx:,1], mid[-idx:,0]-t_s_mid, mid[-idx:,1], endo[-idx:,0]-t_s_endo, endo[-idx:,1],
                    epi2[-idx:,0]-t_s_epi2, epi2[-idx:,1],'--',mid2[-idx:,0]-t_s_mid2, mid2[-idx:,1],'--', endo2[-idx:,0]-t_s_endo2, endo2[-idx:,1],'--')
            ax6.plot(epi[-idx:,0]-t_s_epi, epi[-idx:,4], mid[-idx:,0]-t_s_mid, mid[-idx:,4], endo[-idx:,0]-t_s_endo, endo[-idx:,4],
                    epi2[-idx:,0]-t_s_epi2, epi2[-idx:,4],'--',mid2[-idx:,0]-t_s_mid2, mid2[-idx:,4],'--', endo2[-idx:,0]-t_s_endo2, endo2[-idx:,4],'--')
            if analysis:
                print ('1: Diastolic calcium (nM) epi: '+str(epi[-1,1])+', mid: '+str(mid[-1,1])+', endo: '+str(endo[-1,1]))
                print ('1: Calcium amplitude (nM) epi: '+str(max(epi[-idx:,1]))+', mid: '+str(max(mid[-idx:,1]))+', endo: '+str(max(endo[-idx:,1])))
                print ('1: Diastolic active tension (kPa) epi: '+str(epi[-1,4])+', mid: '+str(mid[-1,4])+', endo: '+str(endo[-1,4]))
                print ('1: Active tension amplitude (kPa) epi: '+str(max(epi[-idx:,4]))+', mid: '+str(max(mid[-idx:,4]))+', endo: '+str(max(endo[-idx:,4])))
                print ('2: Diastolic calcium (nM) epi: '+str(epi2[-1,1])+', mid: '+str(mid2[-1,1])+', endo: '+str(endo2[-1,1]))
                print ('2: Calcium amplitude (nM) epi: '+str(max(epi2[-idx:,1]))+', mid: '+str(max(mid2[-idx:,1]))+', endo: '+str(max(endo2[-idx:,1])))
                print ('2: Diastolic active tension (kPa) epi: '+str(epi2[-1,4])+', mid: '+str(mid2[-1,4])+', endo: '+str(endo2[-1,4]))
                print ('2: Active tension amplitude (kPa) epi: '+str(max(epi2[-idx:,4]))+', mid: '+str(max(mid2[-idx:,4]))+', endo: '+str(max(endo2[-idx:,4])))
        else:
            ax1.plot(epi[:,0], epi[:,2], mid[:,0], mid[:,2], endo[:,0], endo[:,2])
            ax2.plot(epi[:,0],epi[:,1],mid[:,0],mid[:,1],endo[:,0],endo[:,1])
            ax3.plot(epi[:,0],epi[:,4],mid[:,0],mid[:,4],endo[:,0],endo[:,4])
            dt = epi[1,0] - epi[0,0]
            idx = int(self.CL / dt)
            t_s_epi = epi[-idx,0]
            t_s_mid = mid[-idx,0]
            t_s_endo = endo[-idx,0]
            ax4.plot(epi[-idx:,0]-t_s_epi, epi[-idx:,2], mid[-idx:,0]-t_s_mid, mid[-idx:,2], endo[-idx:,0]-t_s_endo, endo[-idx:,2])
            ax5.plot(epi[-idx:,0]-t_s_epi, epi[-idx:,1], mid[-idx:,0]-t_s_mid, mid[-idx:,1], endo[-idx:,0]-t_s_endo, endo[-idx:,1])
            ax6.plot(epi[-idx:,0]-t_s_epi, epi[-idx:,4], mid[-idx:,0]-t_s_mid, mid[-idx:,4], endo[-idx:,0]-t_s_endo, endo[-idx:,4])
            if analysis:
                print ('Diastolic calcium (nM) epi: {:.2f}, mid: {:.2f}, endo: {:.2f}'.format(epi[-1,1], mid[-1,1], endo[-1,1]))
                print ('Calcium amplitude (nM) epi: {:.2f}, mid: {:.2f}, endo: {:.2f}'.format(max(epi[-idx:,1]), max(mid[-idx:,1]), max(endo[-idx:,1])))
                print ('Diastolic active tension (kPa) epi: {:.2f}, mid: {:.2f}, endo: {:.2f}'.format(epi[-1,4], mid[-1,4], endo[-1,4]))
                print ('Active tension amplitude (kPa) epi: {:.2f}, mid: {:.2f}, endo: {:.2f}'.format(max(epi[-idx:,4]), max(mid[-idx:,4]), max(endo[-idx:,4])))

        ax1.set_ylabel('(mV)')
        ax1.set_xlabel('Time (s)')
        #ax1.set_title('Action potentials')
        ax1.legend(['Epi', 'Mid', 'Endo'])

        #ax2.set_ylabel('Intracellular Calcium (nM)')
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('CaT (nM)')

        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Active tension (kPa)')

        dt = epi[1,0] - epi[0,0]
        idx = int(0.8 / dt)

        ax4.set_ylabel('Vm (mV)')
        ax4.set_xlabel('Time (s)')
        ax4.set_ylim([-100, 40])
        ax4.grid()
        ax4.set_xticks(numpy.arange(0, self.CL+0.2, 0.2))

        ax5.set_ylabel('CaT (nM)')
        ax5.set_xlabel('Time (s)')
        ax5.set_ylim([0, 1000])
        ax5.grid()
        ax5.set_xticks(numpy.arange(0, self.CL+0.2, 0.2))

        ax6.set_ylabel('Active tension (kPa)')
        ax6.set_xlabel('Time (s)')
        ax6.set_ylim([0, 70])
        ax6.grid()
        ax6.set_xticks(numpy.arange(0, self.CL+0.2, 0.2))

        if show:
            plt.show(block=True)
        plt.savefig('single_cell_'+str(material)+'.png')


    def analysis_PV(self, pvs, beat=0):
        if beat > 0:
            # More sophisticated method of evaluating EDV and ESV making use of the phase information.
            ejection_vls = pvs['vls'][beat - 1][numpy.where(pvs['phasels'][beat-1] == 2)]
            ejection_vrs = pvs['vrs'][beat - 1][numpy.where(pvs['phasers'][beat - 1] == 2)]
            if len(ejection_vls) > 0:
                ESVL = min(ejection_vls)
                EDVL = max(ejection_vls)
            else:
                ESVL = np.nan
                EDVL = max(pvs['vls'][beat-1])
                plt.plot(pvs['ts'][beat-1], pvs['phasels'][beat-1], pvs['ts'][beat-1], pvs['vls'][beat-1])
                plt.show()
            if len(ejection_vrs) > 0 :
                ESVR = min(ejection_vrs)
                EDVR = max(ejection_vrs)
            else:
                ESVR = np.nan
                EDVR = max(pvs['vrs'][beat - 1])
            # EDVL = max(pvs['vls'][beat-1])
            # EDVR = max(pvs['vrs'][beat-1])
            # ESVL = min(pvs['vls'][beat-1])
            # ESVR = min(pvs['vrs'][beat-1])
            PmaxL = max(pvs['pls'][beat-1])
            PmaxR = max(pvs['prs'][beat-1])
            LVEF = (EDVL-ESVL)/EDVL*100
            RVEF = (EDVR-ESVR)/EDVR*100
            SVL = EDVL - ESVL
            SVR = EDVR - ESVR
        else:
            LVEF = []
            RVEF = []
            EDVL = []
            EDVR = []
            ESVL = []
            ESVR = []
            PmaxL = []
            PmaxR = []
            SVL = []
            SVR = []
            for i in range(0, len(pvs['vls'])):
                EDVL.append(max(pvs['vls'][i]))
                ESVL.append(min(pvs['vls'][i]))
                LVEF.append(int((EDVL[i]-ESVL[i])/EDVL[i] * 100))
                PmaxL.append(int(max(pvs['pls'][i])))
                SVL.append(int(EDVL[i]-ESVL[i]))
                EDVR.append(max(pvs['vrs'][i]))
                ESVR.append(min(pvs['vrs'][i]))
                RVEF.append(int((EDVR[i]-ESVR[i])/EDVR[i] * 100))
                PmaxR.append(int(max(pvs['prs'][i])))
                SVR.append(int(EDVR[i]-ESVR[i]))
        analysis = {'EDVL':EDVL,
            'EDVR':EDVR,
            'ESVL':ESVL,
            'ESVR':ESVR,
            'LVEF':LVEF,
            'RVEF':RVEF,
            'PmaxL':PmaxL,
            'PmaxR':PmaxR,
            'SVL':SVL,
            'SVR':SVR}
        # print('EDVL:'+str(int(EDVL))+' mL,\tEDVR: '+str(int(EDVR))+' mL')
        # print('ESVL: '+str(int(ESVL))+' mL,\tESVR: '+str(int(ESVR))+' mL')
        # print('LVEF: '+str(LVEF)+' %,\tRVEF: '+str(RVEF)+' %')
        # print('PmaxL: '+str(PmaxL)+' kPa,\tPmaxR: '+str(PmaxR)+' kPa')
        # print('SVL: '+str(SVL)+' mL,\tSVR: '+str(SVR)+' mL')
        return analysis

    def analysis_ECG(self,ecgs,beat,show):
        if show:
            fig = plt.figure(tight_layout=True, figsize=[11,9])
            gs = GridSpec(3,4)
            ax1 = fig.add_subplot(gs[0,0])
            ax2 = fig.add_subplot(gs[1,0])
            ax3 = fig.add_subplot(gs[2,0])
            ax4 = fig.add_subplot(gs[0,1])
            ax5 = fig.add_subplot(gs[1,1])
            ax6 = fig.add_subplot(gs[2,1])
            ax7 = fig.add_subplot(gs[0,2])
            ax8 = fig.add_subplot(gs[1,2])
            ax9 = fig.add_subplot(gs[2,2])
            ax10 = fig.add_subplot(gs[0,3])
            ax11 = fig.add_subplot(gs[1,3])
            ax12 = fig.add_subplot(gs[2,3])

            ecg_biomarkers = numpy.zeros((12,6))
            CL = self.CL
            ecg_biomarkers[0, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['Is'][beat-1], 'I', CL, ax=ax1)
            ecg_biomarkers[1, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['IIs'][beat-1], 'II', CL, ax=ax2)
            ecg_biomarkers[2, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['IIIs'][beat-1], 'III', CL, ax=ax3)
            ecg_biomarkers[3, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['aVRs'][beat-1], 'aVR', CL, ax=ax4)
            ecg_biomarkers[4, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['aVLs'][beat-1], 'aVL', CL, ax=ax5)
            ecg_biomarkers[5, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['aVFs'][beat-1], 'aVF', CL, ax=ax6)
            ecg_biomarkers[6, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V1s'][beat-1], 'V1', CL, ax=ax7)
            ecg_biomarkers[7, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V2s'][beat-1], 'V2', CL, ax=ax8)
            ecg_biomarkers[8, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V3s'][beat-1], 'V3', CL, ax=ax9)
            ecg_biomarkers[9, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V4s'][beat-1], 'V4', CL, ax=ax10)
            ecg_biomarkers[10, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V5s'][beat-1], 'V5', CL, ax=ax11)
            ecg_biomarkers[11, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V6s'][beat-1], 'V6', CL, ax=ax12)

            plt.show()
            plt.savefig('ecg_analysis.png')
        else:
            ecg_biomarkers = numpy.zeros((12, 6))
            CL = self.CL
            ecg_biomarkers[0, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['Is'][beat - 1], 'I', CL)
            ecg_biomarkers[1, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['IIs'][beat - 1], 'II', CL)
            ecg_biomarkers[2, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['IIIs'][beat - 1], 'III', CL)
            ecg_biomarkers[3, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['aVRs'][beat - 1], 'aVR', CL)
            ecg_biomarkers[4, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['aVLs'][beat - 1], 'aVL', CL)
            ecg_biomarkers[5, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['aVFs'][beat - 1], 'aVF', CL)
            ecg_biomarkers[6, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['V1s'][beat - 1], 'V1', CL)
            ecg_biomarkers[7, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['V2s'][beat - 1], 'V2', CL)
            ecg_biomarkers[8, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['V3s'][beat - 1], 'V3', CL)
            ecg_biomarkers[9, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                             ecgs['V4s'][beat - 1], 'V4', CL)
            ecg_biomarkers[10, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                              ecgs['V5s'][beat - 1], 'V5', CL)
            ecg_biomarkers[11, :] = self._plot_with_landmarks(ecgs['ts'][beat - 1], ecgs['max_all_leads'],
                                                              ecgs['V6s'][beat - 1], 'V6', CL)

        QRS_mean = numpy.average(ecg_biomarkers[:,0])
        QRS_std = numpy.std(ecg_biomarkers[:,0])
        T_dur_mean = numpy.average(ecg_biomarkers[:,1])
        T_dur_std = numpy.std(ecg_biomarkers[:,1])
        T_pe_mean = numpy.average(ecg_biomarkers[:,2])
        T_pe_std = numpy.std(ecg_biomarkers[:,2])
        T_ep_mean = numpy.average(ecg_biomarkers[:,3])
        T_ep_std = numpy.std(ecg_biomarkers[:,3])
        T_op_mean = numpy.average(ecg_biomarkers[:,4])
        T_op_std = numpy.std(ecg_biomarkers[:,4])
        QT_dur_mean = numpy.average(ecg_biomarkers[:,5])
        QT_dur_std = numpy.std(ecg_biomarkers[:,5])
        QT_dispersion_6 = ecg_biomarkers[6:12,5].max() - ecg_biomarkers[6:12,5].min()
        QT_dispersion_12 = ecg_biomarkers[:,5].max() - ecg_biomarkers[:,5].min()
        analysis = {'QRS':[QRS_mean, QRS_std, ecg_biomarkers[:,0]],
            'T_dur':[T_dur_mean, T_dur_std,ecg_biomarkers[:,1]],
            'T_peak':[T_pe_mean, T_pe_std, ecg_biomarkers[:,2]],
            'T_ep_dur':[T_ep_mean, T_ep_std, ecg_biomarkers[:,3]],
            'T_op_dur':[T_op_mean, T_op_std, ecg_biomarkers[:,4]],
            'QT':[QT_dur_mean, QT_dur_std, ecg_biomarkers[:,5]],
            'QT_dispersion_6':QT_dispersion_6,
            'QT_dispersion_12':QT_dispersion_12}
        # print('QRS: {:.3f} +- {:.3f}'.format(QRS_mean, QRS_std))
        # print('T duration: {:.3f} +- {:.3f}'.format(T_dur_mean, T_dur_std))
        # print('T pe: {:.3f} +- {:.3f}'.format(T_pe_mean, T_pe_std))
        # print('T op: {:.3f} +- {:.3f}'.format(T_op_mean, T_op_std))
        # print('QT duration: {:.3f} +- {:.3f}'.format(QT_dur_mean, QT_dur_std))
        # print('QT dispersion (precordial): {:.3f}'.format(QT_dispersion_6))
        # print('QT dispersion (12 leads): {:.3f}'.format(QT_dispersion_12))
        return analysis

    def analysis_ECG_6leads(self,ecgs,beat,show):
        # fig = plt.figure(tight_layout=True, figsize=[11,9])
        # gs = GridSpec(3,4)
        # ax5 = fig.add_subplot(gs[0,0])
        # ax6 = fig.add_subplot(gs[1,0])
        # ax7 = fig.add_subplot(gs[2,0])
        # ax8 = fig.add_subplot(gs[0,1])
        # ax9 = fig.add_subplot(gs[1,1])
        # ax10 = fig.add_subplot(gs[2,1])
        # ax11 = fig.add_subplot(gs[0,2])
        # ax12 = fig.add_subplot(gs[1,2])
        # ax13 = fig.add_subplot(gs[2,2])
        # ax14 = fig.add_subplot(gs[0,3])
        # ax15 = fig.add_subplot(gs[1,3])
        # ax16 = fig.add_subplot(gs[2,3])

        # ecg_biomarkers = numpy.zeros((8,6))
        # CL = self.CL
        # ecg_biomarkers[0, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_limb_leads'], ecgs['Is'][beat-1], 'I', CL)
        # ecg_biomarkers[1, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_limb_leads'], ecgs['IIs'][beat-1], 'II', CL)
        # ecg_biomarkers[2, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V1s'][beat-1], 'V1', CL)
        # ecg_biomarkers[3, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V2s'][beat-1], 'V2', CL)
        # ecg_biomarkers[4, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V3s'][beat-1], 'V3', CL)
        # ecg_biomarkers[5, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V4s'][beat-1], 'V4', CL)
        # ecg_biomarkers[6, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V5s'][beat-1], 'V5', CL)
        # ecg_biomarkers[7, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V6s'][beat-1], 'V6', CL)

        ecg_biomarkers = numpy.zeros((6, 6))
        CL = self.CL
        ecg_biomarkers[0, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V1s'][beat-1], 'V1', CL)
        ecg_biomarkers[1, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V2s'][beat-1], 'V2', CL)
        ecg_biomarkers[2, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V3s'][beat-1], 'V3', CL)
        ecg_biomarkers[3, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V4s'][beat-1], 'V4', CL)
        ecg_biomarkers[4, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V5s'][beat-1], 'V5', CL)
        ecg_biomarkers[5, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V6s'][beat-1], 'V6', CL)

        if show:
            plt.show()
        plt.savefig('ecg_analysis.png')
        #qrs_dur, t_dur, t_pe, t_op, t_amplitude, qt_dur, landmarks
        QRS_mean = numpy.average(ecg_biomarkers[:,0])
        QRS_std = numpy.std(ecg_biomarkers[:,0])
        T_dur_mean = numpy.average(ecg_biomarkers[:,1])
        T_dur_std = numpy.std(ecg_biomarkers[:,1])
        T_pe_mean = numpy.average(ecg_biomarkers[:,2])
        T_pe_std = numpy.std(ecg_biomarkers[:,2])
        T_op_mean = numpy.average(ecg_biomarkers[:,3])
        T_op_std = numpy.std(ecg_biomarkers[:,3])
        T_amp_mean = numpy.average(ecg_biomarkers[:,4])
        T_amp_std = numpy.std(ecg_biomarkers[:,4])
        QT_dur_mean = numpy.average(ecg_biomarkers[:,5])
        QT_dur_std = numpy.std(ecg_biomarkers[:,5])
        # QT_dispersion_6 = ecg_biomarkers[6:12,5].max() - ecg_biomarkers[6:12,5].min()
        # QT_dispersion_12 = ecg_biomarkers[:,5].max() - ecg_biomarkers[:,5].min()
        analysis = {'QRS':[QRS_mean, QRS_std, ecg_biomarkers[:,0]],
            'T_dur':[T_dur_mean, T_dur_std,ecg_biomarkers[:,1]],
            'T_pe_dur':[T_pe_mean, T_pe_std, ecg_biomarkers[:,2]],
            'T_op_dur':[T_op_mean, T_op_std, ecg_biomarkers[:,3]],
            'T_amp': [T_amp_mean, T_amp_std, ecg_biomarkers[:, 3]],
            'QT':[QT_dur_mean, QT_dur_std, ecg_biomarkers[:,5]]}
            # 'QT_dispersion_6':QT_dispersion_6,
            # 'QT_dispersion_12':QT_dispersion_12}
        # print('QRS: {:.3f} +- {:.3f}'.format(QRS_mean, QRS_std))
        # print('T duration: {:.3f} +- {:.3f}'.format(T_dur_mean, T_dur_std))
        # print('T pe: {:.3f} +- {:.3f}'.format(T_pe_mean, T_pe_std))
        # print('T op: {:.3f} +- {:.3f}'.format(T_op_mean, T_op_std))
        # print('QT duration: {:.3f} +- {:.3f}'.format(QT_dur_mean, QT_dur_std))
        # print('QT dispersion (precordial): {:.3f}'.format(QT_dispersion_6))
        # print('QT dispersion (12 leads): {:.3f}'.format(QT_dispersion_12))
        return analysis

    def analysis_ECG_5leads(self,ecgs,beat,show):
        ecg_biomarkers = numpy.zeros((5, 6))
        CL = self.CL
        ecg_biomarkers[0, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V1s'][beat-1], 'V1', CL)
        ecg_biomarkers[1, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V2s'][beat-1], 'V2', CL)
        ecg_biomarkers[2, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V3s'][beat-1], 'V3', CL)
        ecg_biomarkers[3, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V4s'][beat-1], 'V4', CL)
        ecg_biomarkers[4, :] = self._plot_with_landmarks(ecgs['ts'][beat-1], ecgs['max_all_leads'], ecgs['V5s'][beat-1], 'V5', CL)
        if show:
            plt.show()
        plt.savefig('ecg_analysis.png')
        #qrs_dur, t_dur, t_pe, t_op, t_amplitude, qt_dur, landmarks
        QRS_mean = numpy.average(ecg_biomarkers[:,0])
        QRS_std = numpy.std(ecg_biomarkers[:,0])
        T_dur_mean = numpy.average(ecg_biomarkers[:,1])
        T_dur_std = numpy.std(ecg_biomarkers[:,1])
        T_pe_mean = numpy.average(ecg_biomarkers[:,2])
        T_pe_std = numpy.std(ecg_biomarkers[:,2])
        T_op_mean = numpy.average(ecg_biomarkers[:,3])
        T_op_std = numpy.std(ecg_biomarkers[:,3])
        T_amp_mean = numpy.average(ecg_biomarkers[:,4])
        T_amp_std = numpy.std(ecg_biomarkers[:,4])
        QT_dur_mean = numpy.average(ecg_biomarkers[:,5])
        QT_dur_std = numpy.std(ecg_biomarkers[:,5])
        analysis = {'QRS':[QRS_mean, QRS_std, ecg_biomarkers[:,0]],
            'T_dur':[T_dur_mean, T_dur_std,ecg_biomarkers[:,1]],
            'T_pe_dur':[T_pe_mean, T_pe_std, ecg_biomarkers[:,2]],
            'T_op_dur':[T_op_mean, T_op_std, ecg_biomarkers[:,3]],
            'T_amp': [T_amp_mean, T_amp_std, ecg_biomarkers[:, 3]],
            'QT':[QT_dur_mean, QT_dur_std, ecg_biomarkers[:,5]]}
        return analysis

    def _plot_with_landmarks(self, t, max_all_leads, V, name, CL, ax=None):
    	# qrs_dur, t_dur, t_pe, t_op, t_amplitude, qt_dur, landmarks =self._measurements_with_qrs_dur(V, 0.0, CL, t, 2e-5, 0.01)
        qrs_dur, t_dur, t_pe, t_op, t_amplitude, qt_dur, landmarks = self._measurements_Holmes_Smith(V, t, width=3)
        if ax:
            ax.clear()
            ax.plot(t, V/max_all_leads, landmarks[:,0], landmarks[:,1]/max_all_leads, '*')
            ax.set_title(name + ' '+str(int(qrs_dur*1000))+' '+str(int(qt_dur*1000))+'\n'+str(int(t_dur*1000))+' '+str(int(t_amp/max_all_leads)))
            ax.set_xlim(0.0,CL)
            ax.set_ylim(-1,1)
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Normalised ECG')
        return qrs_dur, t_dur, t_pe, t_op, t_amplitude, qt_dur

    def _measurements_Holmes_Smith(self, V, T, width=3):
        # Resample V and T to 1000 Hz
        V = resample(V, 1000)
        T = resample(T, 1000)
        def get_window(signal, i, width):
            window = signal[(i - width):(i + width + 1)]
            return numpy.mean(window)

        dV = numpy.gradient(V)
        ddV = numpy.gradient(dV)

        T_ex = numpy.concatenate([[-3, -2, -1, 0, 1], T])
        V_ex = numpy.concatenate([[0, 0, 0], V])
        dV_ex = numpy.gradient(V_ex)
        ddV_ex = numpy.gradient(dV_ex)
        dV_windowed = numpy.zeros(300)
        for i in range(width, 301):
            dV_windowed[i - width] = abs(get_window(dV_ex, i, width))

        QRS_start_tol = 0.01 * max(abs(V)) / 30
        QRS_end_tol_ddV = 0.1 * max(abs(V)) / (30 * width + 2)
        QRS_end_tol_dV = 0.07 * max(abs(dV_windowed))
        T_start_tol = 0.12 * max(abs(V)) / 30
        T_end_tol = 0.01 * max(abs(V)) / 30

        # Determine QRS start time
        for i in range(width, 101):
            QRS_window = abs(get_window(ddV_ex, i, width))
            if (QRS_window > QRS_start_tol):
                QRS_start_idx = i - width
                break

        # # Determining QRS end time, and QRS duration
        QRS_window2 = numpy.zeros(501 - 30)
        for i in range(30, 501):
            QRS_window2[i - 30] = get_window(abs(ddV_ex), i, 2)
        QRS_end_idx_lst = 30 + numpy.where(QRS_window2 < QRS_end_tol_ddV)[0] - 1

        for idx in QRS_end_idx_lst:
            dV_window = get_window(abs(dV), idx, width + 2)
            if (abs(dV_window) < QRS_end_tol_dV):
                QRS_end_idx = idx - (width - 2)
                break

        QRS_start_time = T[QRS_start_idx]
        QRS_end_time = T[QRS_end_idx]
        QRS_duration = QRS_end_time - QRS_start_time

        segment = V[QRS_end_idx + 100:1000]  # Assuming ~100 of ST segment
        t_magnitude = max(abs(segment))
        peak_idx = numpy.where(abs(segment) == t_magnitude)[-1][0]
        t_peak_idx = QRS_end_idx + 100 + peak_idx
        t_sign = numpy.sign(V[t_peak_idx])
        t_magnitude_true = t_sign * t_magnitude

        t_wave_peak_time = T[t_peak_idx]

        # Find T-wave end point
        t_wave_end_idx = 0
        for i in range(len(V) - width - 1, t_peak_idx, -1):
            tend_window = abs(get_window(dV, i, width))
            if tend_window > T_end_tol:
                t_wave_end_idx = i
                break
        t_wave_end_time = T[t_wave_end_idx]

        # Find T-wave start point and calculate QT
        try:
            for i in range(t_peak_idx - 20, QRS_end_idx, -1):
                tstart_window = abs(get_window(dV, i, width))
                if tstart_window < T_start_tol:
                    if (max(abs(V[i - 30:i])) < abs(V[i])) & (abs(V[i]) < (0.5 * t_magnitude)):
                        t_wave_start_idx = i
                        break
            t_wave_start_time = T[t_wave_start_idx]
        except:
            t_wave_start_idx = numpy.nan
            t_wave_start_time = numpy.nan

        t_wave_duration = t_wave_end_time - t_wave_start_time
        QT_duration = t_wave_end_time - QRS_start_time
        t_peak_end = t_wave_end_time - t_wave_peak_time
        t_start_peak = t_wave_peak_time - t_wave_start_time

        segment = V[t_wave_end_idx - 10:t_wave_end_idx]
        if max(abs(segment)) > t_magnitude * 0.1:
            if max(segment) > V(t_wave_end_idx):
                inverse = True
            else:
                inverse = False
        else:
            inverse = False
        landmarks = numpy.array(
            [[T[QRS_start_idx], V[QRS_start_idx]], [T[QRS_end_idx], V[QRS_end_idx]], [T[t_peak_idx], V[t_peak_idx]],
             [T[t_wave_end_idx], V[t_wave_end_idx]]])
        return QRS_duration, QRS_duration, t_peak_end, t_start_peak, t_magnitude_true, t_magnitude_true, landmarks

    def _measurements_with_qrs_dur(self, V, start_t, end_t, t, t_tol, v_tol):
        dV = abs(numpy.gradient(V))
        ddV = abs(numpy.gradient(dV))
        dV[0:2] = 0.0  # remove gradient artefacts
        ddV[0:2] = 0.0
        dVTOL_end_of_Twave = 0.0002  # mV/ms
        q_start_idx = 0
        qrs_end_t = 0.1
        qrs_end_idx = numpy.where(abs(t-0.09)<1e-6)[0][0]
        qrs_dur = qrs_end_idx - q_start_idx

        # Find T peak and amplitude
        segment = V[qrs_end_idx:]
        t_amplitude = abs(segment).max()
        t_peak_idx = numpy.where(abs(segment) == t_amplitude)[0][0] + qrs_end_idx
        t_sign = numpy.sign(segment[t_peak_idx - qrs_end_idx])
        t_peak = t_sign * t_amplitude
        t_min = numpy.amin(segment)
        t_max = abs(numpy.amax(segment))
        # t_polarity = t_max/t_min * 1/(max(abs(t_max),abs(t_min))) # Value close to 1 is positive monophasic, close to 0 is negative monophasic, around 0.5 is biphasic.
        t_polarity = (t_max + t_min) / (max(abs(t_max), abs(t_min)))
        # Find T-wave end
        for i in range(len(V) - 1, t_peak_idx, -1):
            if (dV[i] > dVTOL_end_of_Twave):
                break
        t_end_idx = i
        t_dur = 0
        qt_dur = t[t_end_idx] - t[q_start_idx]
        t_pe = t[t_end_idx] - t[t_peak_idx]
        qtpeak_dur = t[t_peak_idx] - t[q_start_idx]
        t_op = 0
        landmarks = numpy.array([[t[q_start_idx], V[q_start_idx]], [t[qrs_end_idx], V[qrs_end_idx]], [t[t_peak_idx], V[t_peak_idx]],
                                 [t[t_end_idx], V[t_end_idx]]])
        return qrs_dur, t_dur, t_pe, t_op, t_amplitude, qt_dur, landmarks

    def _measurements(self, V, start_t, end_t, t, t_tol, v_tol):
    	idx_start = numpy.where(abs(numpy.array(t-start_t)) < t_tol)[0][0]
    	try:
    		idx_end = numpy.where(abs(numpy.array(t-end_t)) < t_tol)[0][0]
    	except:
    		idx_end = len(V)-1
    	# Offset voltage using t_end
    	n = 100
    	b = [1.0 / n] * n
    	V = V - V[idx_end]
    	V = lfilter(b,1, V)
    	dV = numpy.gradient(V)
    	dV[0:2] = 0.0 # Remove gradient artefacts

    	n = 500
    	b = [1.0 / n] * n
    	dV = lfilter(b, 1, dV)
    	dV = abs(dV)/abs(dV).max()
    	ddV = numpy.gradient(dV)
    	ddV[0:2] = 0.0 # Remove gradient artefacts
    	ddV = lfilter(b, 1, ddV)
    	ddV = abs(ddV)/abs(ddV).max()

    	TOL = v_tol

    	# Find Q start
    	for i in range(idx_start, len(V)):
    		if (dV[i] > TOL) & (i > 10):
    			break
    	q_start_idx = i
    	q_start_t = t[i]

    	# Find T end
    	for i in range(idx_end-1, idx_start, -1):
    		if (dV[i] > TOL):
    			break
    	t_end_idx = i
    	t_end_t = t[i]

    	# Find QRS end
    	# window = 5000
    	# dV = abs(dV)
    	# max_window = numpy.zeros(t_end_idx-q_start_idx+1)
    	# for i in range(q_start_idx, t_end_idx):
    	# 	max_window[i-q_start_idx] = dV[i:(i+window)].max()
    	# # print('evaluated max window')
    	# min = max_window[0]
    	# for i in range(0, len(max_window)):
    	# 	if (min >= max_window[i]):
    	# 		min = max_window[i]
    	# 	elif (min < (max_window.min()+0.05)):
    	# 		break
    	qrs_end_idx = i + q_start_idx
    	qrs_end_t = t[qrs_end_idx]

    	qrs_dur = qrs_end_t - q_start_t

    	# Find T peak timing and amplitude
    	segment = V[qrs_end_idx:(t_end_idx+1)]
    	t_amplitude = abs(segment).max()
    	t_peak_idx = numpy.where(abs(segment) == t_amplitude)[0][0] + qrs_end_idx
    	t_sign = numpy.sign(segment[t_peak_idx-qrs_end_idx])
    	t_amplitude = t_sign * t_amplitude
    	t_peak_t = t[t_peak_idx]

    	# # Find T start
    	# segment = ddV[qrs_end_idx:t_peak_idx]
    	# min_dd_idx = numpy.where(segment == segment.min())[0][0] + qrs_end_idx
    	# for i in range(min_dd_idx, t_peak_idx):
    	# 	if (abs(ddV[i]) > 0.0005):
    	# 		break
    	# t_start_idx = i
    	# t_start_t = t[t_start_idx]

    	t_dur = 0 #t_end_t - t_start_t

    	qt_dur = t_end_t - q_start_t
    	t_pe = t_end_t - t_peak_t
    	t_op = 0 #t_peak_t - t_start_t

    	landmarks = numpy.array([[q_start_t, V[q_start_idx]], [qrs_end_t, V[qrs_end_idx]], [t_peak_t, V[t_peak_idx]], [t_end_t, V[t_end_idx]]])
    	return qrs_dur, t_dur, t_pe, t_op, t_amplitude, qt_dur, landmarks


if __name__=='__main__':
    false = ['False', 'false', 'F', 'f']
    true = ['True', 'true', 'T', 't']
    parser = argparse.ArgumentParser()
    parser.add_argument("--type", help='Either ecg or pv, or single cell')
    parser.add_argument("--name", help='Name of simulation', default='heart_remeshed_3D')
    parser.add_argument("--refresh", help='Delete previous reading of ECG and PV', default=False)
    parser.add_argument("--lead", help='Name of the signal lead to plot')
    parser.add_argument("--CL", help='Cycle length of simulation, (s)',default=0.8)
    parser.add_argument("--beat", help='Beat number for single beat plots', default=0)
    parser.add_argument("--figure_title", help='Title of the figure', default = '')
    parser.add_argument("--analysis", help='Title of the figure', default =False)
    parser.add_argument("--compare", help='Directory to compare with', default='')
    parser.add_argument("--material", help='Directory to compare with', default=1)
    parser.add_argument("--show", help='Toggle whether to show plot', default=True)
    args = parser.parse_args()
    plot_type = args.type
    name = args.name
    refresh = bool(args.refresh)
    if refresh in true:
        os.system('rm ecgs.pl pvs.pl')
    lead_name = args.lead
    CL = float(args.CL)
    beat = int(args.beat)
    figure_title = args.figure_title
    analysis = args.analysis
    compare = args.compare
    material = int(args.material)
    show = args.show

    if show in false:
        show = False
    if analysis in true:
        analysis = True

    # Set up visualisation
    a = ECGPV_visualisation(CL)
    if not (plot_type=='cell'):
        ecgs, pvs = a.read_ecg_pv(name, './')
    else:
        epi, mid, endo = a.read_single_cell('./', material)
    if compare:
        ecgs2, pvs2 = a.read_ecg_pv(name, compare)
        epi2, mid2, endo2 = a.read_single_cell(compare, material)
        if analysis:
            if (plot_type=='live')|(plot_type=='p')|(plot_type=='v')|(plot_type=='pv'):
                a.analysis_PV(pvs, beat)
                print('second one...')
                a.analysis_PV(pvs2, beat)
            elif (plot_type=='ecg'):
                a.analysis_ECG(ecgs,beat,show)
                a.analysis_ECG(ecgs2,beat,show)
        if plot_type == 'live':
            a.plot_ecgpv_live(ecgs, pvs, figure_title, show, ecgs2, pvs2)
        elif plot_type == 'ecg':
            if lead_name:
                a.plot_ecg_lead(ecgs, lead_name, lead_name, show, beat, ecgs2)
            else:
                a.plot_ecg_all_leads(ecgs, show, beat, ecgs2)
        elif plot_type == 'p':
            a.plot_pv_signal(pvs, 'p', 'Pt', show, beat, pvs2)
        elif plot_type == 'v':
            a.plot_pv_signal(pvs, 'v', 'Vt', show, beat, pvs2)
        elif plot_type == 'pv':
            a.plot_pv_signal(pvs, 'pv', 'PV_loop', show, beat, pvs2)
        elif plot_type == 'cell':
            a.plot_single_cell(epi, mid, endo, material,analysis, show, epi2, mid2, endo2)
    else:
        if analysis:
            if (plot_type=='live')|(plot_type=='p')|(plot_type=='v')|(plot_type=='pv'):
                a.analysis_PV(pvs, beat)
            elif (plot_type=='ecg'):
                a.analysis_ECG(ecgs,beat,show)
        if plot_type == 'live':
            a.plot_ecgpv_live(ecgs, pvs, figure_title, show)
        elif plot_type == 'ecg':
            if lead_name:
                assert lead_name, 'Please enter lead name using option: lead='
                a.plot_ecg_lead(ecgs, lead_name, lead_name, show, beat)
            else:
                a.plot_ecg_all_leads(ecgs, show, beat)
        elif plot_type == 'p':
            a.plot_pv_signal(pvs, 'p', 'Pt', show, beat)
        elif plot_type == 'v':
            a.plot_pv_signal(pvs, 'v', 'Vt', show, beat)
        elif plot_type == 'pv':
            a.plot_pv_signal(pvs,'pv', 'PV_loop', show, beat)
        elif plot_type == 'cell':
            a.plot_single_cell(epi, mid, endo, material, analysis, show)
