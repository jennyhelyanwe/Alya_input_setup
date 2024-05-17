import matplotlib
# matplotlib.use('TKagg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import os
import json
# from mpi4py import MPI #needs install
# from multiprocessing import Queue
import pymp, multiprocessing
from ECGPV_visualisation import ECGPV_visualisation
from meshstructure import MeshStructure
from myformat import Fields
import pandas as pd
import healthy_qoi_ranges
from scipy.signal import resample

class PostProcessing(MeshStructure):
    def __init__(self, alya, simulation_json_file, alya_output_dir, protocol, verbose):
        # super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        if hasattr(alya, 'geometry'):
            self.geometry = alya.geometry
        self.name = alya.name
        if hasattr(alya, 'node_fields'):
            self.node_fields = alya.node_fields
        if hasattr(alya, 'element_fields'):
            self.element_fields = alya.element_fields
        self.alya_output_dir = alya_output_dir
        self.alyacsv_dir = self.alya_output_dir+'/results_csv/'
        self.results_dir = self.alya_output_dir + '/results_analysis/'
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)
        self.post_nodefield = Fields(self.name, field_type='postnodefield', verbose=verbose)
        self.post_elementfield = Fields(self.name, field_type='postelementfield', verbose=verbose)
        if protocol == 'postprocess':
            print ('Reading in post_nodefield CSV from: ', self.results_dir)
            self.post_nodefield.read_csv_to_attributes(input_dir=self.results_dir, field_type='postnodefield')

            print('Reading in post_elementfield CSV from: ', self.results_dir)
            self.post_elementfield.read_csv_to_attributes(input_dir=self.results_dir, field_type='postelementfield')
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.verbose = verbose
        self.qoi = {}
        self.healthy_ranges = healthy_qoi_ranges.HealthyBiomarkerRanges().healthy_ranges
        # print('Read Alya log outputs')
        # self.read_log_files()
        # super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)

        # if 'INTRA' not in self.post_nodefield.dict.keys():
        #     self.read_csv_fields('INTRA', 'scalar')
        # if 'lat' not in self.post_nodefield.dict.keys():
        #     self.evaluate_ep_maps()
        # # self.save_postprocessing()
        # self.evaluate_ecg_biomarkers()
        # quit()
        # print('Reading postprocessing csv fields')
        # self.post_elementfield = Fields(name, field_type='postelementfield', verbose=verbose)
        # # self.post_elementfield.read_csv_to_attributes(input_dir=results_dir, field_type='postelementfield')
        # # Evaluate various biomarkers
        # if 'INTRA' not in self.post_nodefield.dict.keys():
        #     self.read_csv_fields()
        #     self.save_postprocessing()
        # if 'lat' not in self.post_nodefield.dict.keys():
        #     self.evaluate_ep_maps()
        #     self.save_postprocessing()
        # self.qoi = {}
        # self.evaluate_ventricular_cvs()
        # self.evaluate_ep_maps()
        #
        # self.evaluate_ventricular_cvs()
        # if 'SOLIDZ' in self.simulation_dict['physics']:
        #     self.evaluate_basal_displacement()
        # self.save_postprocessing()
        # print('Writing out simulation biomarkers to '+self.results_dir+'qoi.json')
        # print(self.qoi)
        # json.dump(self.qoi, open(self.results_dir+'qoi.json', 'w'))

    def save_qoi(self, filename):
        # results = pd.DataFrame(self.qoi.items(), columns=self.qoi.keys())
        print('Saving QoIs to ', filename)
        json.dump(self.qoi, open(filename, 'w'))

    def save_postprocessing_fields(self):
        # Save to CSV
        print('Saving postprocessing fields to CSV at '+self.results_dir)
        # self.geometry.save_to_csv(self.results_dir)
        self.post_nodefield.save_to_csv(self.results_dir)
        self.post_elementfield.save_to_csv(self.results_dir)
        # Save to ensight
        # print('Saving geometry and fields to Ensight at ' + self.results_dir + 'ensight/')
        # if not os.path.exists(self.results_dir+'ensight/'):
        #     os.mkdir(self.results_dir+'ensight/')
        # # self.geometry.save_to_ensight(self.results_dir + 'ensight/')
        # self.post_nodefield.save_to_ensight(output_dir=self.results_dir + 'ensight/',
        #                                     casename=self.name + '_postnodefield',
        #                                     geometry=self.geometry)
        # self.post_elementfield.save_to_ensight(output_dir=self.results_dir + 'ensight/',
        #                                        casename=self.name + '_postelementfield',
        #                                        geometry=self.geometry)

    def read_ecg_pv(self):
        ecgpv = ECGPV_visualisation(CL=1.0)
        ecg_filename = self.alya_output_dir + self.simulation_dict['name']
        print('Reading ECG from: ', ecg_filename)
        self.ecgs = ecgpv._read_ECG(ecg_filename)
        if 'SOLIDZ' in self.simulation_dict['physics']:
            pv_filename = self.alya_output_dir + self.simulation_dict['name']
            print('Reading PV from: ', pv_filename)
            self.pvs = ecgpv._read_PV(pv_filename)

    def evaluate_ecg_biomarkers(self, beat, show_landmarks=True):
        # if 'INTRA' not in self.post_nodefield.dict.keys():
        #     self.read_csv_fields(read_field_name='INTRA', read_field_type='scalar')
        #     self.evaluate_ep_maps()
        #     # self.save_postprocessing_fields()
        # if 'lat' not in self.post_nodefield.dict.keys():
        #     self.evaluate_ep_maps()
        #     self.save_postprocessing_fields()
        ecgpv = ECGPV_visualisation(CL=1.0)
        ecg_filename = self.alya_output_dir + self.simulation_dict['name']
        print('Reading ECG from: ', ecg_filename)
        self.ecgs = ecgpv._read_ECG(ecg_filename)
        lead_names = ['Is', 'IIs', 'IIIs', 'aVLs', 'aVRs', 'aVFs', 'V1s', 'V2s', 'V3s', 'V4s', 'V5s', 'V6s']
        qrs_dur = np.zeros((len(lead_names)))
        qt_dur = np.zeros((len(lead_names)))
        t_pe = np.zeros((len(lead_names)))
        t_op = np.zeros((len(lead_names)))
        t_peak = np.zeros((len(lead_names)))
        qtpeak_dur = np.zeros((len(lead_names)))
        t_dur = np.zeros((len(lead_names)))
        for i, lead_i in enumerate(lead_names):
            # qrs_dur[i], qt_dur[i], t_pe[i], t_peak[i], qtpeak_dur[i], t_polarity[i], landmarks_temp = \
            #     self.calculate_ecg_biomarkers(time=self.ecgs['ts'][beat-1], V=self.ecgs[lead_i][beat-1],
            #                                   LAT=self.post_nodefield.dict['lat'])
            # qrs_dur[i], qt_dur[i], t_pe[i], t_peak[i], qtpeak_dur[i], t_polarity[i], landmarks_temp = \
            #         self.calculate_ecg_biomarkers(time=self.ecgs['ts'][beat-1], V=self.ecgs[lead_i][beat-1],
            #                                       qrs_end_t=0.1)
            # QRS_duration, QT_duration, QTpeak_duration, t_wave_duration, t_peak_end, t_start_peak, t_magnitude_true,  landmarks
            qrs_dur[i], qt_dur[i],qtpeak_dur[i], t_dur[i], t_pe[i], t_op[i], t_peak[i], landmarks_temp = \
                        self.calculate_ecg_biomarkers_HolmesSmith(T=self.ecgs['ts'][beat-1], V=self.ecgs[lead_i][beat-1], show=show_landmarks)

        qoi = {}
        qoi['qrs_dur_mean'] = np.mean(qrs_dur)
        qoi['qt_dur_mean'] = np.mean(qt_dur)
        qoi['qtpeak_dur_mean'] = np.mean(qtpeak_dur)
        qoi['t_dur_mean'] = np.mean(t_dur)
        qoi['t_pe_mean'] = np.mean(t_pe)
        qoi['t_op_mean'] = np.mean(t_op)
        qoi['t_peak_mean'] = np.mean(t_peak)

        qoi['qrs_dur_std'] = np.std(qrs_dur)
        qoi['qt_dur_std'] = np.std(qt_dur)
        qoi['qtpeak_dur_std'] = np.std(qtpeak_dur)
        qoi['t_dur_std'] = np.std(t_dur)
        qoi['t_pe_std'] = np.std(t_pe)
        qoi['t_op_std'] = np.std(t_op)
        qoi['t_peak_std'] = np.std(t_peak)
        self.qoi.update(qoi)

    def evaluate_pv_biomarkers(self, beat):
        pv_filename = self.alya_output_dir + self.simulation_dict['name']
        ecgpv = ECGPV_visualisation(CL=1.0)
        print('Reading PV from: ', pv_filename)
        self.pvs = ecgpv._read_PV(pv_filename)
        pv_analysis = ecgpv.analysis_PV(pvs=self.pvs, beat=beat)

        def get_first_derivative(self, t, y):
            dydt = [0]
            for i in range(len(y) - 1):
                i += 1
                dydt.append((y[i] - y[i - 1]) / (t[i] - t[i - 1]))
            assert len(dydt) == len(y)
            return np.array(dydt)
        # Get volume flow rates
        t = self.pvs['ts'][beat-1]
        y = self.pvs['vls'][beat-1]
        end_systole_idx = np.argmin(y)
        dvdt = self.get_first_derivative(t=t, y=y)
        if end_systole_idx > 10:
            dvdt_ejection = abs(np.amin(dvdt[10:end_systole_idx]))
            dvdt_filling = abs(np.amax(dvdt[end_systole_idx:]))
        else:
            dvdt_ejection = np.nan
            dvdt_filling = np.nan

        t = self.pvs['ts'][beat - 1]
        y = self.pvs['pls'][beat - 1]
        end_systole_idx = np.argmax(y)
        dpdt = self.get_first_derivative(t=t, y=y)
        if end_systole_idx > 10:
            dpdt_max = abs(np.amax(dvdt[10:end_systole_idx]))
        else:
            dpdt_max = np.nan

        qoi = {}
        qoi['EDVL'] = pv_analysis['EDVL']
        qoi['EDVR'] = pv_analysis['EDVR']
        qoi['ESVL'] = pv_analysis['ESVL']
        qoi['ESVR'] = pv_analysis['ESVR']
        qoi['LVEF'] = pv_analysis['LVEF']
        qoi['RVEF'] = pv_analysis['RVEF']
        qoi['PmaxL'] = pv_analysis['PmaxL']/10000
        qoi['PmaxR'] = pv_analysis['PmaxR']/10000
        qoi['SVL'] = pv_analysis['SVL']
        qoi['SVR'] = pv_analysis['SVR']
        qoi['dvdt_ejection'] = dvdt_ejection
        qoi['dvdt_filling'] = dvdt_filling
        qoi['dpdt_max'] = dpdt_max
        self.qoi.update(qoi)

    def evaluate_deformation_biomarkers(self, beat):
        self.deformation_transients = {}
        # self.read_csv_fields(read_field_name='DISPL', read_field_type='vector')
        self.read_binary_outputs(read_field_name='DISPL', read_field_type='vector')
        # Select time segment according to specified beat # TODO assuming CL is always 1.0 s
        CL = 1.0
        time_idx = []
        deformation_t = []
        for time_i in range(self.post_nodefield.dict['time'].shape[0]):
            if (self.post_nodefield.dict['time'][time_i] > (beat - 1) * CL) & (
                    self.post_nodefield.dict['time'][time_i] < beat * CL):
                deformation_t.append(self.post_nodefield.dict['time'][time_i] - (beat - 1) * CL)
                time_idx.append(time_i)
        deformation_t = np.array(deformation_t)
        time_idx = np.array(time_idx, dtype=int)

        # AVPD atrioventricular plane displacement
        print('Evaluating AVPD and apical displacement')
        base_cutoff = 0.8
        truncated_mesh_nodes = np.nonzero((self.node_fields.dict['ab'] < base_cutoff) & (self.node_fields.dict['tv'] == self.geometry.lv))[0]
        basal_mesh_nodes = np.nonzero((self.node_fields.dict['ab'] >= base_cutoff)& (self.node_fields.dict['tv'] == self.geometry.lv))[0]
        mean_ab_vector = np.mean(self.node_fields.dict['longitudinal-vector'][truncated_mesh_nodes, :], axis=0)
        apical_cutoff = 0.2
        apical_mesh_nodes = np.nonzero((self.node_fields.dict['ab'] <= apical_cutoff)& (self.node_fields.dict['tv'] == self.geometry.lv))[0]
        avpd = pymp.shared.array(time_idx.shape[0], dtype=float)
        apical_displacement = pymp.shared.array(time_idx.shape[0], dtype=float)
        displacement_shared = pymp.shared.array(self.post_nodefield.dict['DISPL'].shape, dtype=float)
        displacement_shared[:,:,:] = self.post_nodefield.dict['DISPL']

        # # Debug
        # ab_delineation = np.zeros(self.node_fields.dict['ab'].shape[0])
        # ab_delineation[basal_mesh_nodes] = 1
        # ab_delineation[apical_mesh_nodes] = 2
        # self.post_nodefield.add_field(data=ab_delineation, data_name='ABDLN', field_type='postnodefield')
        # self.post_nodefield.save_to_ensight(output_dir=self.alya_output_dir, casename=self.name, geometry=self.geometry, fieldname='ABDLN', fieldtype='postnodefield')
        #
        threadsNum = int(multiprocessing.cpu_count() * 0.7)
        with pymp.Parallel(min(threadsNum, time_idx.shape[0])) as p1:
            for time_i in p1.range(time_idx.shape[0]):
                avpd[time_i] = np.mean(np.dot(displacement_shared[
                                              basal_mesh_nodes, :, time_i],
                                              mean_ab_vector))
                apical_displacement[time_i] = np.mean(np.dot(displacement_shared[
                                                  apical_mesh_nodes, :, time_i],
                                                  mean_ab_vector))
        # Wall thickness
        print('Evaluating wall thickness')
        lv_nodes = np.nonzero((self.node_fields.dict['tv'] == self.geometry.lv)
                              & (self.node_fields.dict['rvlv'] > 0.2))[0]  # Exclude also the septum from this.
        rv_nodes = np.nonzero(self.node_fields.dict['tv'] == self.geometry.rv)[0]
        lv_endo_nodes = lv_nodes[np.nonzero(self.node_fields.dict['tm'][lv_nodes] == self.geometry.tm_endo)[0]]
        lv_epi_nodes = lv_nodes[np.nonzero(self.node_fields.dict['tm'][lv_nodes] == self.geometry.tm_epi)[0]]
        rv_endo_nodes = rv_nodes[np.nonzero(self.node_fields.dict['tm'][rv_nodes] == self.geometry.tm_endo)[0]]
        rv_epi_nodes = rv_nodes[np.nonzero(self.node_fields.dict['tm'][rv_nodes] == self.geometry.tm_epi)[0]]
        lv_mapped_epi_nodes = lv_epi_nodes[mapIndices(points_to_map_xyz=self.geometry.nodes_xyz[lv_endo_nodes, :],
                                                      reference_points_xyz=self.geometry.nodes_xyz[lv_epi_nodes, :])]
        rv_mapped_epi_nodes = rv_epi_nodes[mapIndices(points_to_map_xyz=self.geometry.nodes_xyz[rv_endo_nodes, :],
                                                      reference_points_xyz=self.geometry.nodes_xyz[rv_epi_nodes, :])]
        lv_wall_thickness = pymp.shared.array(time_idx.shape[0])
        rv_wall_thickness = pymp.shared.array(time_idx.shape[0])
        nodes_xyz = pymp.shared.array(self.geometry.nodes_xyz.shape, dtype=float)
        nodes_xyz[:,:] = self.geometry.nodes_xyz
        with pymp.Parallel(min(threadsNum, time_idx.shape[0])) as p1:
            for time_i in p1.range(time_idx.shape[0]):
                lv_update_epi_coords = nodes_xyz[lv_mapped_epi_nodes, :] + displacement_shared[lv_mapped_epi_nodes, :, time_i]
                lv_update_endo_coords = nodes_xyz[lv_endo_nodes, :] + displacement_shared[lv_endo_nodes, :, time_i]
                lv_transmural_vector = lv_update_epi_coords - lv_update_endo_coords
                lv_wall_thickness[time_i] = np.mean(np.linalg.norm(lv_transmural_vector, axis=1))

                rv_update_epi_coords = nodes_xyz[rv_mapped_epi_nodes, :] + displacement_shared[rv_mapped_epi_nodes, :, time_i]
                rv_update_endo_coords = nodes_xyz[rv_endo_nodes, :] + displacement_shared[rv_endo_nodes, :, time_i]
                rv_transmural_vector = rv_update_epi_coords - rv_update_endo_coords
                rv_wall_thickness[time_i] = np.mean(np.linalg.norm(rv_transmural_vector, axis=1))

        # Save transients
        self.deformation_transients['deformation_t'] = deformation_t
        self.deformation_transients['avpd'] = avpd
        self.deformation_transients['apical_displacement'] = apical_displacement
        self.deformation_transients['lv_wall_thickness'] = lv_wall_thickness
        self.deformation_transients['rv_wall_thickness'] = rv_wall_thickness
        # Save QoIs
        qoi = {}
        qoi['peak_avpd'] = np.amax(avpd)
        qoi['min_avpd'] = np.amin(avpd)
        qoi['es_ed_avpd'] = np.amax(avpd) - np.amin(avpd)
        qoi['peak_apical_displacement'] = np.amax(apical_displacement)
        qoi['min_apical_displacement'] = np.amin(apical_displacement)
        qoi['es_ed_apical_displacement'] = np.amax(apical_displacement) - np.amin(apical_displacement)
        qoi['peak_lv_wall_thickness'] = np.amax(lv_wall_thickness)
        qoi['min_lv_wall_thickness'] = np.amin(lv_wall_thickness)
        qoi['peak_rv_wall_thickness'] = np.amax(rv_wall_thickness)
        qoi['min_rv_wall_thickness'] = np.amin(rv_wall_thickness)
        qoi['diff_lv_wall_thickness'] = np.amax(lv_wall_thickness) - np.amin(lv_wall_thickness)
        qoi['diff_rv_wall_thickness'] = np.amax(rv_wall_thickness) - np.amin(rv_wall_thickness)
        self.qoi.update(qoi)

    def evaluate_cube_deformation_ta_biomarkers(self):
        self.deformation_transients = {}
        nodes_xyz = np.loadtxt(self.alyacsv_dir + '3D_point_coordinates.csv', delimiter=',')
        number_of_nodes = nodes_xyz.shape[0]
        # self.read_csv_fields(read_field_name='DISPL', read_field_type='vector', nodes_xyz=nodes_xyz)
        # self.read_csv_fields(read_field_name='ACTST', read_field_type='vector', nodes_xyz=nodes_xyz)
        self.read_binary_outputs(read_field_name='DISPL', read_field_type='vector', nodes_xyz=nodes_xyz)
        self.read_binary_outputs(read_field_name='ACTST', read_field_type='vector', nodes_xyz=nodes_xyz)

        # Select nodes on x=1 face
        x1nodes = []
        tol = 1e-5
        for node_i in range(number_of_nodes):
            if abs(nodes_xyz[node_i,0] - 2.0) < tol:
                x1nodes.append(node_i)

        # Get displacement in the x axis as transient over time.
        t = self.post_nodefield.dict['time']
        displ_x = pymp.shared.array((len(x1nodes), len(t)))
        ta = pymp.shared.array((len(x1nodes), len(t)))
        displacement_shared = pymp.shared.array(self.post_nodefield.dict['DISPL'].shape, dtype=float)
        displacement_shared[:, :, :] = self.post_nodefield.dict['DISPL']
        ta_shared = pymp.shared.array(self.post_nodefield.dict['ACTST'].shape, dtype=float)
        ta_shared[:, :, :] = self.post_nodefield.dict['ACTST']
        x1nodes_shared = pymp.shared.array(len(x1nodes), dtype='int')
        x1nodes_shared[:] = x1nodes
        threadsNum = int(multiprocessing.cpu_count()*0.7)
        with pymp.Parallel(min(threadsNum, t.shape[0])) as p1:
            for time_i in p1.range(t.shape[0]):
                for node_i in range(x1nodes_shared.shape[0]):
                    displ_x[node_i, time_i] = displacement_shared[x1nodes_shared[node_i],0, time_i]
                    ta[node_i, time_i] = ta_shared[x1nodes_shared[node_i], 0, time_i]/10000 # convert from Barye to kPa
        mean_displ_x = np.mean(displ_x, axis=0)
        mean_ta = np.mean(ta, axis=0)
        self.deformation_transients['deformation_t'] = t
        self.deformation_transients['mean_displ_x'] = mean_displ_x
        self.deformation_transients['mean_ta'] = mean_ta

        # Get gradients of x displacement
        peak_displ_x = np.amax(abs(mean_displ_x))
        dxdt = self.get_first_derivative(t, mean_displ_x)
        rise_dxdt = np.amax(dxdt)
        decay_dxdt = np.amin(dxdt)

        # Get gradients of Ta transient
        peak_ta = np.amax(abs(mean_ta))
        dtadt = self.get_first_derivative(t, mean_ta)
        rise_dtadt = np.amax(dtadt)
        decay_dtadt = np.amin(dtadt)

        qoi = {}
        qoi['rise_dxdt'] = rise_dxdt
        qoi['decay_dxdt'] = decay_dxdt
        qoi['peak_displ_x'] = peak_displ_x
        qoi['peak_ta'] = peak_ta
        qoi['rise_dtadt'] = rise_dtadt
        qoi['decay_dtadt'] = decay_dtadt
        self.qoi.update(qoi)

    def evaluate_fibre_work_biomarkers(self, beat):
        self.fibre_work = {}
        # self.read_csv_fields(read_field_name='LAMBD', read_field_type='vector')
        # self.read_csv_fields(read_field_name='ACTST', read_field_type='vector')
        self.read_binary_outputs(read_field_name='LAMBD', read_field_type='vector')
        self.read_binary_outputs(read_field_name='ACTST', read_field_type='vector')
        # Select time segment according to specified beat # TODO assuming CL is always 1.0 s
        CL = 1.0
        time_idx = []
        t = []
        for time_i in range(self.post_nodefield.dict['time'].shape[0]):
            if (self.post_nodefield.dict['time'][time_i] > (beat - 1) * CL) & (
                    self.post_nodefield.dict['time'][time_i] < beat * CL):
                t.append(self.post_nodefield.dict['time'][time_i] - (beat - 1) * CL)
                time_idx.append(time_i)
        fibrework_t = np.array(t)
        time_idx = np.array(time_idx, dtype=int)

        # Separate mesh into three longitudinal chunks
        base_cutoff = 0.66
        basal_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] >= base_cutoff)[0]
        apical_cutoff = 0.33
        apical_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] <= apical_cutoff)[0]
        midv_mesh_nodes = np.nonzero((self.node_fields.dict['ab'] > apical_cutoff) & (self.node_fields.dict['ab'] < base_cutoff))[0]

        endo_cutoff = 0.7
        epi_cutoff = 0.3
        endo_mesh_nodes = np.nonzero(self.node_fields.dict['tm'] >= endo_cutoff)[0]
        epi_mesh_nodes = np.nonzero(self.node_fields.dict['tm'] <= epi_cutoff)[0]
        midw_mesh_nodes = np.nonzero((self.node_fields.dict['tm'] < endo_cutoff) & (self.node_fields.dict['tm'] > epi_cutoff))[0]

        mid_short_axis_nodes = np.nonzero(self.node_fields.dict['short-axis-slices'] == 2)[0]
        endo_midshort_nodes = mid_short_axis_nodes[np.nonzero(self.node_fields.dict['tm'][mid_short_axis_nodes]>=endo_cutoff)[0]]
        epi_midshort_nodes = mid_short_axis_nodes[np.nonzero(self.node_fields.dict['tm'][mid_short_axis_nodes]<= epi_cutoff)[0]]
        mid_midshort_nodes = mid_short_axis_nodes[np.nonzero((self.node_fields.dict['tm'][mid_short_axis_nodes] < endo_cutoff)
                                                             & (self.node_fields.dict['tm'][mid_short_axis_nodes] > epi_cutoff))[0]]


        # Apical lambda and Ta
        apical_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)
        midv_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)
        basal_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)
        endo_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)
        midw_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)
        epi_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)

        apical_Ta = pymp.shared.array(time_idx.shape[0], dtype=float)
        midv_Ta = pymp.shared.array(time_idx.shape[0], dtype=float)
        basal_Ta = pymp.shared.array(time_idx.shape[0], dtype=float)
        endo_Ta = pymp.shared.array(time_idx.shape[0], dtype=float)
        midw_Ta = pymp.shared.array(time_idx.shape[0], dtype=float)
        epi_Ta = pymp.shared.array(time_idx.shape[0], dtype=float)

        endo_cross_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)
        midw_cross_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)
        epi_cross_lambda = pymp.shared.array(time_idx.shape[0], dtype=float)

        endo_midshort_lambda = pymp.shared.array((endo_midshort_nodes.shape[0], time_idx.shape[0]))
        epi_midshort_lambda = pymp.shared.array((epi_midshort_nodes.shape[0], time_idx.shape[0]))
        mid_midshort_lambda = pymp.shared.array((mid_midshort_nodes.shape[0], time_idx.shape[0]))

        ta_shared = pymp.shared.array(self.post_nodefield.dict['ACTST'].shape, dtype=float)
        lambda_shared = pymp.shared.array(self.post_nodefield.dict['LAMBD'].shape, dtype=float)
        ta_shared[:,:] = self.post_nodefield.dict['ACTST']
        lambda_shared[:,:] = self.post_nodefield.dict['LAMBD']
        threadsNum = int(multiprocessing.cpu_count()*0.7)
        print('Regionally dividing lambda and Ta')
        with pymp.Parallel(min(threadsNum, time_idx.shape[0])) as p1:
            for time_i in p1.range(time_idx.shape[0]):
                # apical_lambda[time_i] = np.mean(lambda_shared[apical_mesh_nodes, 0, time_i])
                # midv_lambda[time_i] = np.mean(lambda_shared[midv_mesh_nodes, 0, time_i])
                # basal_lambda[time_i] = np.mean(lambda_shared[basal_mesh_nodes, 0, time_i])
                # endo_lambda[time_i] = np.mean(lambda_shared[endo_mesh_nodes, 0, time_i])
                # midw_lambda[time_i] = np.mean(lambda_shared[midw_mesh_nodes, 0, time_i])
                # epi_lambda[time_i] = np.mean(lambda_shared[epi_mesh_nodes, 0, time_i])
                # endo_cross_lambda[time_i] = np.mean(lambda_shared[endo_mesh_nodes, 1, time_i])
                # midw_cross_lambda[time_i] = np.mean(lambda_shared[midw_mesh_nodes, 1, time_i])
                # epi_cross_lambda[time_i] = np.mean(lambda_shared[epi_mesh_nodes, 1, time_i])

                apical_Ta[time_i] = np.mean(ta_shared[apical_mesh_nodes, 0, time_i])
                midv_Ta[time_i] = np.mean(ta_shared[midv_mesh_nodes, 0, time_i])
                basal_Ta[time_i] = np.mean(ta_shared[basal_mesh_nodes, 0, time_i])
                endo_Ta[time_i] = np.mean(ta_shared[endo_mesh_nodes, 0, time_i])
                midw_Ta[time_i] = np.mean(ta_shared[midw_mesh_nodes, 0, time_i])
                epi_Ta[time_i] = np.mean(ta_shared[epi_mesh_nodes, 0, time_i])

                endo_midshort_lambda[:, time_i] = lambda_shared[endo_midshort_nodes, 0, time_i]
                epi_midshort_lambda[:, time_i] = lambda_shared[epi_midshort_nodes, 0, time_i]
                mid_midshort_lambda[:, time_i] = lambda_shared[mid_midshort_nodes, 0, time_i]
        systole = np.min(endo_midshort_lambda, axis=1)
        tol = 0.01 * abs(np.min(systole))
        endo_midshort_lambda_median = endo_midshort_lambda[np.where(abs(systole - np.percentile(systole, 50)) < tol)[0][0], :]
        endo_midshort_lambda_uq = endo_midshort_lambda[np.where(abs(systole - np.percentile(systole, 75)) < tol)[0][0], :]
        endo_midshort_lambda_lq = endo_midshort_lambda[np.where(abs(systole - np.percentile(systole, 25)) < tol)[0][0], :]

        systole = np.min(epi_midshort_lambda, axis=1)
        tol = 0.01 * abs(np.min(systole))
        epi_midshort_lambda_median = epi_midshort_lambda[
                                      np.where(abs(systole - np.percentile(systole, 50)) < tol)[0][0], :]
        epi_midshort_lambda_uq = epi_midshort_lambda[np.where(abs(systole - np.percentile(systole, 75)) < tol)[0][0],
                                  :]
        epi_midshort_lambda_lq = epi_midshort_lambda[np.where(abs(systole - np.percentile(systole, 25)) < tol)[0][0],
                                  :]

        systole = np.min(mid_midshort_lambda, axis=1)
        tol = 0.01 * abs(np.min(systole))
        mid_midshort_lambda_median = mid_midshort_lambda[
                                      np.where(abs(systole - np.percentile(systole, 50)) < tol)[0][0], :]
        mid_midshort_lambda_uq = mid_midshort_lambda[np.where(abs(systole - np.percentile(systole, 75)) < tol)[0][0],
                                  :]
        mid_midshort_lambda_lq = mid_midshort_lambda[np.where(abs(systole - np.percentile(systole, 25)) < tol)[0][0],
                                  :]

        self.fibre_work_transients = {}
        self.fibre_work_transients['fibrework_t'] = fibrework_t
        # self.fibre_work['apical_lambda'] = apical_lambda
        # self.fibre_work['midv_lambda'] = midv_lambda
        # self.fibre_work['basal_lambda'] = basal_lambda
        # self.fibre_work['endo_lambda'] = endo_lambda
        # self.fibre_work['midw_lambda'] = midw_lambda
        # self.fibre_work['epi_lambda'] = epi_lambda
        # self.fibre_work['endo_cross_lambda'] = endo_cross_lambda
        # self.fibre_work['midw_cross_lambda'] = midw_cross_lambda
        # self.fibre_work['epi_cross_lambda'] = epi_cross_lambda

        self.fibre_work_transients['endo_midshort_lambda']= endo_midshort_lambda
        self.fibre_work_transients['endo_midshort_lambda_median'] = endo_midshort_lambda_median
        self.fibre_work_transients['endo_midshort_lambda_uq'] = endo_midshort_lambda_uq
        self.fibre_work_transients['endo_midshort_lambda_lq'] = endo_midshort_lambda_lq

        self.fibre_work_transients['epi_midshort_lambda'] = epi_midshort_lambda
        self.fibre_work_transients['epi_midshort_lambda_median'] = epi_midshort_lambda_median
        self.fibre_work_transients['epi_midshort_lambda_uq'] = epi_midshort_lambda_uq
        self.fibre_work_transients['epi_midshort_lambda_lq'] = epi_midshort_lambda_lq

        self.fibre_work_transients['mid_midshort_lambda'] = mid_midshort_lambda
        self.fibre_work_transients['mid_midshort_lambda_median'] = mid_midshort_lambda_median
        self.fibre_work_transients['mid_midshort_lambda_uq'] = mid_midshort_lambda_uq
        self.fibre_work_transients['mid_midshort_lambda_lq'] = mid_midshort_lambda_lq

        self.fibre_work_transients['mean_lambda'] = np.mean(lambda_shared[:,0,:], axis=0)[time_idx]
        self.fibre_work_transients['mean_Ta'] = np.mean(ta_shared[:,0,:], axis=0)[time_idx]

        self.fibre_work_transients['apical_Ta'] = apical_Ta
        self.fibre_work_transients['midv_Ta'] = midv_Ta
        self.fibre_work_transients['basal_Ta'] = basal_Ta
        self.fibre_work_transients['endo_Ta'] = endo_Ta
        self.fibre_work_transients['midw_Ta'] = midw_Ta
        self.fibre_work_transients['epi_Ta'] = epi_Ta
        # Save QoIs
        qoi = {}
        # qoi['peak_Ta'] = np.amax(np.mean(ta_shared[:,0, :], axis=0))
        # qoi['min_Ta'] = np.amin(np.mean(ta_shared[:,0, :], axis=0))
        qoi['peak_lambda'] = np.amax(np.mean(lambda_shared[:,0, :], axis=0))
        qoi['min_lambda'] = np.amin(np.mean(lambda_shared[:, 0, :], axis=0))
        qoi['peak_ta'] = np.amax(np.mean(ta_shared[:,0,:], axis=0))
        mean_ta = np.mean(ta_shared[:,0,:], axis=0)
        qoi['diastolic_ta'] = np.amin(mean_ta[np.nonzero(mean_ta)])
        self.qoi.update(qoi)

    def calculate_ecg_biomarkers_HolmesSmith(self, T, V, width=3, show=False):
        # Resample V and T to 1000 Hz
        V = resample(V, 1000)
        T = resample(T, 1000)

        def get_window(signal, i, width):
            window = signal[(i - width):(i + width + 1)]
            return np.mean(window)

        dV = np.gradient(V)
        ddV = np.gradient(dV)

        T_ex = np.concatenate([[-3, -2, -1, 0, 1], T])
        V_ex = np.concatenate([[0, 0, 0], V])
        dV_ex = np.gradient(V_ex)
        ddV_ex = np.gradient(dV_ex)
        dV_windowed = np.zeros(300)
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
        QRS_window2 = np.zeros(501 - 30)
        for i in range(30, 501):
            QRS_window2[i - 30] = get_window(abs(ddV_ex), i, 2)
        QRS_end_idx_lst = 30 + np.where(QRS_window2 < QRS_end_tol_ddV)[0] - 1

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
        peak_idx = np.where(abs(segment) == t_magnitude)[-1][0]
        t_peak_idx = QRS_end_idx + 100 + peak_idx
        t_sign = np.sign(V[t_peak_idx])
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
            t_wave_start_idx = np.nan
            t_wave_start_time = np.nan

        t_wave_duration = t_wave_end_time - t_wave_start_time
        QT_duration = t_wave_end_time - QRS_start_time
        t_peak_end = t_wave_end_time - t_wave_peak_time
        t_start_peak = t_wave_peak_time - t_wave_start_time
        QTpeak_duration = t_wave_peak_time - QRS_start_time

        segment = V[t_wave_end_idx - 10:t_wave_end_idx]
        if max(abs(segment)) > t_magnitude * 0.1:
            if max(segment) > V(t_wave_end_idx):
                inverse = True
            else:
                inverse = False
        else:
            inverse = False
        landmarks = np.array(
            [[T[QRS_start_idx], V[QRS_start_idx]], [T[QRS_end_idx], V[QRS_end_idx]], [T[t_peak_idx], V[t_peak_idx]],
             [T[t_wave_end_idx], V[t_wave_end_idx]]])

        if show:

            plt.plot(T, V , landmarks[:, 0], landmarks[:, 1], '*')
            plt.title(str(int(QRS_duration * 1000)) + ' ' + str(int(QT_duration * 1000)) + '\n' + str(
                int(t_wave_duration * 1000)))
            plt.xlabel('Time (s)')
            plt.ylabel('Normalised ECG')
            plt.show()

        return QRS_duration, QT_duration, QTpeak_duration, t_wave_duration, t_peak_end, t_start_peak, t_magnitude_true,  landmarks

    def calculate_ecg_biomarkers(self, time, V, LAT=None, qrs_end_t=None):
        """This is not to be used with clinical data because it assumes that the signal returns to baseline after the end
            of the T wave."""
        # TODO: Assumes ECGs are at 1000 Hz.
        dV = abs(np.gradient(V))
        ddV = abs(np.gradient(dV))
        dV[0:2] = 0.0  # remove gradient artefacts
        ddV[0:2] = 0.0
        # Find Q start
        dVTOL_end_of_Twave = 0.0002  # mV/ms
        # dVTOL_start_of_QRS = 0.0002
        ddVTOL_start_of_Twave = 0.0002  # mV/ms^2

        # for i in range(len(V)):
        #     if (dV[i] > dVTOL_start_of_QRS) & (i > 10):
        #         break
        # q_start_idx = i


        # Set QRS end
        if LAT:
            q_start_idx = np.nanargmin(LAT)
            qrs_end_idx = np.nanargmax(LAT).astype(int)
        elif qrs_end_t:
            q_start_idx = 0
            qrs_end_idx = np.argmin(abs(time-0.1)) # TODO Use actual LAT to calculate this!
        qrs_dur = time[qrs_end_idx] - time[q_start_idx]  # TODO code how to find end of QRS

        # Find T peak and amplitude
        segment = V[qrs_end_idx:]
        t_amplitude = abs(segment).max()
        t_peak_idx = np.where(abs(segment) == t_amplitude)[0][0] + qrs_end_idx
        t_sign = np.sign(segment[t_peak_idx - qrs_end_idx])
        t_peak = t_sign * t_amplitude
        t_min = np.amin(segment)
        t_max = abs(np.amax(segment))
        # t_polarity = t_max/t_min * 1/(max(abs(t_max),abs(t_min))) # Value close to 1 is positive monophasic, close to 0 is negative monophasic, around 0.5 is biphasic.
        t_polarity = (t_max + t_min) / (max(abs(t_max), abs(t_min)))
        # Find T-wave end
        i = len(V) - 1
        for i in range(len(V) - 1, t_peak_idx, -1):
            if (dV[i] > dVTOL_end_of_Twave):
                break
        t_end_idx = i

        # Find T start
        # segment = ddV[qrs_end_idx:t_peak_idx]
        # min_dd_idx= np.where(segment == segment.min())[0][0] + qrs_end_idx
        # for i in range(min_dd_idx, t_peak_idx):
        #     if (abs(ddV[i] > ddVTOL_start_of_Twave)):
        #         break
        # t_start_idx = i

        # t_dur = t_end_idx - t_start_idx
        qt_dur = time[t_end_idx] - time[q_start_idx]
        t_pe = time[t_end_idx] - time[t_peak_idx]
        qtpeak_dur = time[t_peak_idx] - time[q_start_idx]
        # t_op = t_start_idx - t_peak_idx
        landmarks = np.array([[q_start_idx, V[q_start_idx]], [qrs_end_idx, V[qrs_end_idx]], [t_peak_idx, V[t_peak_idx]],
                              [t_end_idx, V[t_end_idx]]])
        return qrs_dur, qt_dur, t_pe, t_peak, qtpeak_dur, t_polarity, landmarks

    def read_binary_outputs(self, read_field_name, read_field_type, nodes_xyz=None):
        inputfolder = self.alya_output_dir
        print('Reading ' + read_field_name + ' alyabins from : ' + self.alya_output_dir)
        project_name = 'heart'
        alya_id_type = identify_alya_id_type(inputfolder, project_name)
        file_suffix = '.post.alyabin'
        # read the partitioning info
        partition_filename = os.path.join(inputfolder, project_name + '.post.alyapar')
        with open(partition_filename) as f:
            partitions = np.fromstring(f.read(), dtype=alya_id_type, sep=' ')

        # describes the mesh partitions
        # partition_id,  NumberOfElementsInPartition,  NumberOfPointsInPartition, NumberOfBoundariesInPartition
        partitions = np.reshape(partitions[1:], (partitions[0], 4))
        number_of_blocks = partitions.shape[0]

        LNINV = read_alyabin_array(os.path.join(inputfolder, project_name + '-LNINV' + file_suffix),
                                   number_of_blocks=number_of_blocks, alya_id_type=alya_id_type)
        inverse_pt_correspondence = (LNINV['tuples'] - 1).ravel()  # convert ids to python
        # verify the point correspondence
        max_id = inverse_pt_correspondence.max()
        pt_ids = np.zeros(max_id + 1)
        pt_ids[inverse_pt_correspondence] = 1
        assert not (LNINV['tuples'] < 0).any(), "Negative elements in LNINV, wrong mesh"
        assert (pt_ids > 0).all(), "Some points in the mesh do not have a correspondence in the parittions"
        pt_ids = None  # free memory

        with open(os.path.join(inputfolder, project_name + '.post.alyafil'), 'r') as f:
            field_filelist = f.read().splitlines()

        # remove spaces and empty lines
        field_filelist = [x.strip() for x in field_filelist if x.strip() != '']
        new_field_filelist = []
        fields = []
        iteration_numbers = []
        for filename in field_filelist:
            s1 = filename[len(project_name):].split('-')
            fields = fields + [s1[1]]
            iteration_numbers = iteration_numbers + [int(s1[2].split('.')[0])]  # this will be long in python 3
            new_field_filelist = new_field_filelist + [filename]

        variable_info = pd.DataFrame(
            {'field': fields, 'iteration': iteration_numbers, 'filename': new_field_filelist})
        variable_info['time_int'] = 0
        variable_info['time_real'] = 0
        variable_info['variabletype'] = ''
        variable_info['association'] = ''

        alya_id_type = identify_alya_id_type(inputfolder, project_name)
        if nodes_xyz:
            num_nodes = nodes_xyz.shape[0]
        else:
            num_nodes = self.geometry.number_of_nodes
        if read_field_name:
            if read_field_type == 'scalar':
                # temp = np.zeros((num_nodes, 1, variable_info[variable_info.field == read_field_name].shape[0]))
                # time = np.zeros(variable_info[variable_info.field == read_field_name].shape[0])
                temp = pymp.shared.array((num_nodes, variable_info[variable_info.field == read_field_name].shape[0]), dtype=float)
                time = pymp.shared.array(variable_info[variable_info.field == read_field_name].shape[0], dtype=float)
                result = read_alyabin_array(
                    filename=self.alya_output_dir + variable_info[variable_info.field == read_field_name].iloc[0][
                        'filename'],
                    number_of_blocks=number_of_blocks, alya_id_type=alya_id_type)
                threadsNum = int(multiprocessing.cpu_count()*0.7)
                with pymp.Parallel(min(threadsNum, variable_info[variable_info.field == read_field_name].shape[0])) as p1:
                    for i in p1.range(variable_info[variable_info.field == read_field_name].shape[0]):
                # for i in range(variable_info[variable_info.field == read_field_name].shape[0]):
                        result = read_alyabin_array(
                            filename=self.alya_output_dir + variable_info[variable_info.field == read_field_name].iloc[i][
                                'filename'],
                            number_of_blocks=number_of_blocks, alya_id_type=alya_id_type)
                        temp[inverse_pt_correspondence, i] = result['tuples'].ravel()
                        time[i] = result['header']['Time']

            elif read_field_type == 'vector':
                # temp = np.zeros((num_nodes, 3, variable_info[variable_info.field == read_field_name].shape[0]))
                # time = np.zeros(variable_info[variable_info.field == read_field_name].shape[0])
                # for i in range(variable_info[variable_info.field == read_field_name].shape[0]):
                temp = pymp.shared.array(
                    (num_nodes, 3, variable_info[variable_info.field == read_field_name].shape[0]), dtype=float)
                time = pymp.shared.array(variable_info[variable_info.field == read_field_name].shape[0],
                                         dtype=float)
                threadsNum = int(multiprocessing.cpu_count()*0.7)
                with pymp.Parallel(
                        min(threadsNum, variable_info[variable_info.field == read_field_name].shape[0])) as p1:
                    for i in p1.range(variable_info[variable_info.field == read_field_name].shape[0]):
                        result = read_alyabin_array(
                            filename=self.alya_output_dir + variable_info[variable_info.field == read_field_name].iloc[i][
                                'filename'],
                            number_of_blocks=number_of_blocks, alya_id_type=alya_id_type)
                        ncomponents = result['tuples'].shape[1]
                        temp[inverse_pt_correspondence, 0:ncomponents, i] = result['tuples']
                        time[i] = result['header']['Time']
            else:
                print('postprocessing.py: Read field type ' + read_field_type + ' not recognised. ')
                quit()
        else:
            print('Read_CSV_FIELD: Must specify field and type')
        self.post_nodefield.add_field(data=temp, data_name=read_field_name, field_type='postnodefield')
        self.post_nodefield.add_field(data=time, data_name='time', field_type='postnodefield')
        temp = None
        time = None

    def read_csv_fields(self, read_field_name, read_field_type, nodes_xyz=None):
        name = self.simulation_dict['name']
        self.post_nodefield.add_field(data=pd.read_csv(self.alyacsv_dir + 'timeset_1.csv', delimiter=',',header=None).values.transpose()[0].astype(float),
                                      data_name='time', field_type='postnodefield')
        time = self.post_nodefield.dict['time']
        time_index = pd.read_csv(self.alyacsv_dir + 'timeindices_1.csv', delimiter=',',header=None).values.transpose()[0].astype(int)
        assert time.shape == time_index.shape
        time_index_shared = pymp.shared.array(time_index.shape, dtype=int)
        time_index_shared[:] = time_index[:]
        filenames = np.array([f for f in os.listdir(self.alyacsv_dir) if
                              os.path.isfile(os.path.join(self.alyacsv_dir, f)) and ('scalar' in f or 'vector' in f)])
        field_names_types = []
        for filename in filenames:
            field_names_types.append([filename.split('.')[2].split('-')[0], filename.split('.')[1]])
        field_names_types = np.unique(field_names_types, axis=0)
        if nodes_xyz is None:
            number_of_nodes = self.geometry.number_of_nodes
        else:
            number_of_nodes = nodes_xyz.shape[0]
        if read_field_name:
            temp = np.zeros((number_of_nodes, time.shape[0]))
            field_name = read_field_name
            field_type = read_field_type
            print('Reading csv for field: '  + 'from : '+self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-******.csv')
            if field_type == 'scalar':
                temp = pymp.shared.array((number_of_nodes, time.shape[0]), dtype=float)
                threadsNum = int(multiprocessing.cpu_count()*0.7)
                with pymp.Parallel(min(threadsNum, time_index_shared.shape[0])) as p1:
                    for time_i in p1.range(time_index_shared.shape[0]):
                # temp = np.zeros((number_of_nodes, time.shape[0]))
                # print(temp.shape)
                # if True:
                #     for time_i in range(time_index.shape[0]):
                        index = '{:06d}'.format(time_index_shared[time_i])
                        filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
                        temp[:, time_i] = pd.read_csv(filename, delimiter=',', header=None).values.transpose()[0] #.astype(float)
                self.post_nodefield.dict[field_name] = temp
            elif field_type == 'vector':
                temp = pymp.shared.array((number_of_nodes, 3, time.shape[0]), dtype=float)
                threadsNum = int(multiprocessing.cpu_count()*0.7)
                with pymp.Parallel(min(threadsNum, time_index_shared.shape[0])) as p1:
                    for time_i in p1.range(time_index_shared.shape[0]):
                        index = '{:06d}'.format(time_index_shared[time_i])
                        filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
                        temp[:, :, time_i] = pd.read_csv(filename, delimiter=',', header=None).values #.astype(float)
                self.post_nodefield.dict[field_name] = temp
            self.post_nodefield.add_field(data=np.array(temp), data_name=field_name, field_type='postnodefield')

        else:
            print('Read_CSV_FIELD: Must specify field and type')
        # else:
        #     for field_i in range(field_names_types.shape[0]):
        #         temp = np.zeros((number_of_nodes, time.shape[0]))
        #         field_name = field_names_types[field_i, 0]
        #         field_type = field_names_types[field_i, 1]
        #         print('Reading csv for field: '+field_name + ', of type: '+field_type)
        #         if field_type == 'scalar':
        #             temp = pymp.shared.array((self.geometry.number_of_nodes, time.shape[0]), dtype=float)
        #             threadsNum = multiprocessing.cpu_count()
        #             with pymp.Parallel(min(threadsNum, time_index.shape[0])) as p1:
        #                 for time_i in p1.range(time_index.shape[0]):
        #                     index = '{:06d}'.format(time_index[time_i])
        #                     filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
        #                     temp[:, time_i] = np.loadtxt(filename, delimiter=',').astype(float)
        #             self.post_nodefield.dict[field_name] = temp
        #         elif field_type == 'vector':
        #             temp = pymp.shared.array((self.geometry.number_of_nodes, 3, time.shape[0]), dtype=float)
        #             threadsNum = multiprocessing.cpu_count()
        #             with pymp.Parallel(min(threadsNum, time_index.shape[0])) as p1:
        #                 for time_i in p1.range(time_index.shape[0]):
        #                     index = '{:06d}'.format(time_index[time_i])
        #                     filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
        #                     temp[:, :, time_i] = np.loadtxt(filename, delimiter=',').astype(float)
        #             self.post_nodefield.dict[field_name] = temp
        #         self.post_nodefield.add_field(data=np.array(temp), data_name=field_name, field_type='postnodefield')

    def evaluate_ep_maps(self):
        print('Evaluating LAT and RT')
        lat = evaluate_lat(time=self.post_nodefield.dict['time'], vm=self.post_nodefield.dict['INTRA'], percentage=0.7,
                           time_window=[float(self.simulation_dict['exmedi_delay_time']),
                                        float(self.simulation_dict['cycle_length'])])
        rt = evaluate_rt(time=self.post_nodefield.dict['time'], vm=self.post_nodefield.dict['INTRA'], percentage=0.9,
                         time_window=[float(self.simulation_dict['exmedi_delay_time']),
                                      float(self.simulation_dict['cycle_length'])])
        rt[self.node_fields.dict['tv'] == -10] = np.nan
        lat[self.node_fields.dict['tv'] == -10] = np.nan
        print('Max LAT: ', str(np.nanmax(lat)))
        print('Max RT: ', str(np.nanmax(rt)))
        self.post_nodefield.add_field(data=lat, data_name='lat', field_type='postnodefield')
        self.post_nodefield.add_field(data=rt, data_name='rt', field_type='postnodefield')

    def evaluate_ventricular_cvs(self):
        print('Evaluating mean transmural conduction velocity')
        lv_nodes = np.nonzero((self.node_fields.dict['tv'] == self.geometry.lv)
                              & (self.node_fields.dict['rvlv'] > 0.2))[0] # Exclude also the septum from this.
        rv_nodes = np.nonzero(self.node_fields.dict['tv'] == self.geometry.rv)[0]
        lv_endo_nodes = lv_nodes[np.nonzero(self.node_fields.dict['tm'][lv_nodes] == self.geometry.tm_endo)[0]]
        lv_epi_nodes = lv_nodes[np.nonzero(self.node_fields.dict['tm'][lv_nodes] == self.geometry.tm_epi)[0]]
        rv_endo_nodes = rv_nodes[np.nonzero(self.node_fields.dict['tm'][rv_nodes] == self.geometry.tm_endo)[0]]
        rv_epi_nodes = rv_nodes[np.nonzero(self.node_fields.dict['tm'][rv_nodes] == self.geometry.tm_epi)[0]]
        lv_mapped_epi_nodes = lv_epi_nodes[mapIndices(points_to_map_xyz=self.geometry.nodes_xyz[lv_endo_nodes, :],
                                         reference_points_xyz=self.geometry.nodes_xyz[lv_epi_nodes, :])]
        rv_mapped_epi_nodes = rv_epi_nodes[mapIndices(points_to_map_xyz=self.geometry.nodes_xyz[rv_endo_nodes, :],
                                         reference_points_xyz=self.geometry.nodes_xyz[rv_epi_nodes, :])]
        lv_transmural_vector = self.geometry.nodes_xyz[lv_mapped_epi_nodes, :] - self.geometry.nodes_xyz[lv_endo_nodes, :]
        rv_transmural_vector = self.geometry.nodes_xyz[rv_mapped_epi_nodes, :] - self.geometry.nodes_xyz[rv_endo_nodes, :]

        # transmural_csv_mapping = np.zeros((self.geometry.number_of_nodes, 3))
        # transmural_csv_mapping[lv_endo_nodes, :] = lv_transmural_vector
        # transmural_csv_mapping[rv_endo_nodes, :] = rv_transmural_vector
        # self.post_nodefield.add_field(data=transmural_csv_mapping, data_name='transmural_cvs_mapping',
        #                               field_type='postnodefield')
        dlat_lv = self.post_nodefield.dict['lat'][lv_mapped_epi_nodes] - self.post_nodefield.dict['lat'][lv_endo_nodes]
        dlat_rv = self.post_nodefield.dict['lat'][rv_mapped_epi_nodes] - self.post_nodefield.dict['lat'][rv_endo_nodes]
        lv_transmural_cvs = np.linalg.norm(lv_transmural_vector, axis=1) / dlat_lv # [cm/s]
        rv_transmural_cvs = np.linalg.norm(rv_transmural_vector, axis=1) / dlat_rv
        cvs = np.zeros(self.geometry.number_of_nodes)
        cvs[lv_endo_nodes] = lv_transmural_cvs
        cvs[rv_endo_nodes] = rv_transmural_cvs
        wall_thicnkess = np.zeros(self.geometry.number_of_nodes)
        wall_thicnkess[lv_endo_nodes] = np.linalg.norm(lv_transmural_vector, axis=1)
        wall_thicnkess[rv_endo_nodes] = np.linalg.norm(rv_transmural_vector, axis=1)
        self.post_nodefield.add_field(data=cvs, data_name='transmural-cv', field_type='postnodefield')
        self.post_nodefield.add_field(data=wall_thicnkess, data_name='wall-thickness', field_type='postnodefield')
        qoi = {}
        qoi['mean_lv_transmural_cv'] = np.ma.masked_invalid(lv_transmural_cvs).mean()
        qoi['mean_rv_transmural_cv'] = np.ma.masked_invalid(rv_transmural_cvs).mean()
        self.qoi.update(qoi)



    def evaluate_basal_displacement(self):
        print('Basal displacement in the longitudinal direction')
        base_cutoff = 0.7
        truncated_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] < base_cutoff)[0]
        basal_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] >= base_cutoff) [0]
        mean_ab_vector = np.mean(self.node_fields.dict['longitudinal-vector'][truncated_mesh_nodes,:], axis=0)
        apical_cutoff = 0.2
        apical_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] <= apical_cutoff)[0]
        reference_apex_to_base_length = np.amax(np.dot(self.geometry.nodes_xyz[basal_mesh_nodes, :], mean_ab_vector)) - \
                                        np.amin(np.dot(self.geometry.nodes_xyz[apical_mesh_nodes, :], mean_ab_vector))

        mean_basal_ab_displacement_transient = np.zeros(self.post_nodefield.dict['time'].shape[0])
        mean_apical_ab_displacement_transient = np.zeros(self.post_nodefield.dict['time'].shape[0])
        mean_apicobasal_sum_displacement_transient = np.zeros(self.post_nodefield.dict['time'].shape[0])
        for time_i in range(self.post_nodefield.dict['time'].shape[0]):
            mean_basal_ab_displacement_transient[time_i] = np.mean(np.dot(self.post_nodefield.dict['DISPL'][
                                                                          basal_mesh_nodes, :, time_i],
                                                                          mean_ab_vector))
            mean_apical_ab_displacement_transient[time_i] = np.mean(np.dot(self.post_nodefield.dict['DISPL'][
                                                                           apical_mesh_nodes, :, time_i],
                                                                           mean_ab_vector))
            mean_apicobasal_sum_displacement_transient[time_i] = mean_basal_ab_displacement_transient[time_i] + mean_apical_ab_displacement_transient[time_i]
        qoi = {}
        qoi['max_basal_ab_displacement'] = np.amax(mean_basal_ab_displacement_transient)
        qoi['min_basal_ab_displacement'] = np.amin(mean_basal_ab_displacement_transient)
        qoi['max_apical_ab_displacement'] = np.amax(mean_apical_ab_displacement_transient)
        qoi['min_apical_ab_displacement'] = np.amin(mean_apical_ab_displacement_transient)
        qoi['max_longitudinal_strain'] = np.amax(mean_apicobasal_sum_displacement_transient)/reference_apex_to_base_length
        qoi['mean_longitudinal_strain'] = np.mean(mean_apicobasal_sum_displacement_transient)/reference_apex_to_base_length
        qoi['min_longitudinal_strain'] = np.amin(mean_apicobasal_sum_displacement_transient)/reference_apex_to_base_length
        self.qoi.update(qoi)
        np.savetxt(self.results_dir + 'mean_basal_ab_displacement_transient.txt', mean_basal_ab_displacement_transient, delimiter=',')
        np.savetxt(self.results_dir + 'mean_apical_ab_displacement_transient.txt', mean_apical_ab_displacement_transient, delimiter=',')
        np.savetxt(self.results_dir + 'mean_apicobasal_sum_displacement_transient.txt',
                   mean_apicobasal_sum_displacement_transient, delimiter=',')
        plt.figure()
        plt.plot(self.post_nodefield.dict['time'], mean_basal_ab_displacement_transient,
                 self.post_nodefield.dict['time'], mean_apical_ab_displacement_transient)
        plt.xlabel('Time (s)')
        plt.ylabel('Longitudinal displacement (cm)')
        plt.legend(['Base', 'Apex'])
        plt.savefig(self.results_dir + 'apicobasal_displacement_transient.png')
        plt.show()

    def evaluate_torsion(self):
        print('Evaluate LV torsion using two short axis slices ')
        # TODO Need to extract radial vectors and calculate rotation angles for basala nd apical slices separately
        # https://link.springer.com/article/10.1186/1532-429X-14-49 using the Russel et al formula
        # thetaCL = (phi_apex * r_apex - phi_base * r_base) / D
        lv_nodes = np.nonzero((self.node_fields.dict['tv'] == self.geometry.lv))[0]
        rv_nodes = np.nonzero(self.node_fields.dict['tv'] == self.geometry.rv)[0]
        apical_cutoff = 0.3
        apical_slice_nodes = lv_nodes[np.nonzero((abs(self.node_fields.dict['ab'][lv_nodes]-apical_cutoff)) < 0.01)[0]]
        basal_slice_nodes = lv_nodes[np.nonzero((abs(self.node_fields.dict['ab'][lv_nodes]-base_cutoff)) < 0.01)[0]]
        basal_mapped_nodes = basal_slice_nodes[mapIndices(points_to_map_xyz=self.geometry.nodes_xyz[apical_slice_nodes,:],
                                         reference_points_xyz=self.geometry.nodes_xyz[basal_slice_nodes,:])]
        apical_centre = np.mean(self.geometry.nodes_xyz[apical_slice_nodes, :], axis=0)
        basal_centre = np.mean(self.geometry.nodes_xyz[basal_mapped_nodes, :], axis=0)
        longitudinal_slice_vectors = self.geometry.nodes_xyz[basal_mapped_nodes,:] - self.geometry.nodes_xyz[apical_slice_nodes,:]
        print(apical_centre)
        print(basal_centre)
        global_slices_nodes = np.zeros(self.geometry.number_of_nodes)
        global_longitudinal_slice_vectors = np.zeros((self.geometry.number_of_nodes, 3))
        global_longitudinal_slice_vectors[apical_slice_nodes, :] = longitudinal_slice_vectors
        global_slices_nodes[apical_slice_nodes] = 1
        global_slices_nodes[basal_slice_nodes] = 2
        self.post_nodefield.add_field(data=global_slices_nodes, data_name='torsion-shortaxis-slices',
                                      field_type='postnodefield')
        self.post_nodefield.add_field(data=global_longitudinal_slice_vectors, data_name='torsion-longitudinal-vectors',
                                      field_type='postnodefield')
        # resting_longitudinal_vectors = self.geometry.nodes_xyz[apical_mapped_nodes, :] - \
        #                        self.geometry.nodes_xyz[basal_slice_nodes, :]
        # mean_torsion_angle_transient = np.zeros(self.post_nodefield.dict['time'].shape[0])
        # for time_i in range(self.post_nodefield.dict['time'].shape[0]):
        #     updated_apical_slice_coord = self.geometry.nodes_xyz[apical_mapped_nodes, :] + \
        #                                  self.post_nodefield.dict['DISPL'][apical_mapped_nodes, :, time_i]
        #     updated_basal_slice_coord = self.geometry.nodes_xyz[basal_slice_nodes, :] + \
        #                                 self.post_nodefield.dict['DISPL'][basal_slice_nodes, :, time_i]
        #     updated_longitudinal_vectors = updated_apical_slice_coord - updated_basal_slice_coord
        #     angle = np.zeros(updated_longitudinal_vectors.shape[0])
        #     for node_i in range(updated_longitudinal_vectors.shape[0]):
        #         dot_product = np.dot(updated_longitudinal_vectors[node_i,:], resting_longitudinal_vectors[node_i,:])
        #         norm = np.linalg.norm(resting_longitudinal_vectors[node_i, :], axis=0) * \
        #                np.linalg.norm(updated_longitudinal_vectors[node_i, :], axis=0)
        #         angle[node_i] = np.degrees(np.arccos(dot_product/norm))
        #     mean_torsion_angle_transient[time_i] = np.mean(angle)
        # np.savetxt(self.results_dir +'torsion_angle_transient.txt', mean_torsion_angle_transient, delimiter=',')
        # plt.figure()
        # plt.plot(self.post_nodefield.dict['time'], mean_torsion_angle_transient)
        # plt.xlabel('Time (s)')
        # plt.ylabel('Torsion angle in degrees')
        # plt.savefig(self.results_dir + 'torsion_angle_transient.png')
        #
        # print('Wall thickness at diastasis, end diastole, and end systole')
        # lv_nodes = np.nonzero((self.node_fields.dict['tv'] == self.geometry.lv))[0]  # Exclude also the septum from this.
        # rv_nodes = np.nonzero(self.node_fields.dict['tv'] == self.geometry.rv)[0]
        # lv_endo_nodes = lv_nodes[np.nonzero(self.node_fields.dict['tm'][lv_nodes] == self.geometry.tm_endo)[0]]
        # lv_epi_nodes = lv_nodes[np.nonzero(self.node_fields.dict['tm'][lv_nodes] == self.geometry.tm_epi)[0]]
        # rv_endo_nodes = rv_nodes[np.nonzero(self.node_fields.dict['tm'][rv_nodes] == self.geometry.tm_endo)[0]]
        # rv_epi_nodes = rv_nodes[np.nonzero(self.node_fields.dict['tm'][rv_nodes] == self.geometry.tm_epi)[0]]
        # lv_mapped_epi_nodes = lv_epi_nodes[mapIndices(points_to_map_xyz=self.geometry.nodes_xyz[lv_endo_nodes, :],
        #                                               reference_points_xyz=self.geometry.nodes_xyz[lv_epi_nodes, :])]
        # rv_mapped_epi_nodes = rv_epi_nodes[mapIndices(points_to_map_xyz=self.geometry.nodes_xyz[rv_endo_nodes, :],
        #                                               reference_points_xyz=self.geometry.nodes_xyz[rv_epi_nodes, :])]
        # lv_wall_thickness_transient = np.zeros(self.post_nodefield.dict['time'].shape[0])
        # rv_wall_thickness_transient = np.zeros(self.post_nodefield.dict['time'].shape[0])
        # for time_i in range(self.post_nodefield.dict['time'].shape[0]):
        #     updated_node_coords = self.geometry.nodes_xyz[:, :] + self.post_nodefield.dict['DISPL'][:, :, time_i]
        #     lv_wall_thickness_transient[time_i] = np.mean(np.linalg.norm(updated_node_coords[lv_mapped_epi_nodes,:] -
        #                                                                  updated_node_coords[lv_endo_nodes, :], axis=1))
        #     rv_wall_thickness_transient[time_i] = np.mean(np.linalg.norm(updated_node_coords[rv_mapped_epi_nodes, :] -
        #                                                                  updated_node_coords[rv_endo_nodes, :], axis=1))
        # np.savetxt(self.results_dir + 'lv_wall_thickness_transient.txt', lv_wall_thickness_transient, delimiter=',')
        # np.savetxt(self.results_dir + 'rv_wall_thickness_transient.txt', rv_wall_thickness_transient, delimiter=',')
        # plt.figure()
        # plt.plot(self.post_nodefield.dict['time'], lv_wall_thickness_transient, self.post_nodefield.dict['time'], rv_wall_thickness_transient)
        # plt.xlabel('Time (s)')
        # plt.ylabel('Averaged wall thickness (cm)')
        # plt.legend(['LV', 'RV'])
        # plt.savefig(self.results_dir + 'wall_thicknesss_transient.png')
        # diastasis_time_idx = np.nonzero(self.post_nodefield.dict['time'] == self.simulation_dict['diastasis_t'])
        # self.qoi['mean_lv_thickness_diastasis'] = lv_wall_thickness_transient[diastasis_time_idx]
        # self.qoi['mean_rv_thickness_diastasis'] = rv_wall_thickness_transient[diastasis_time_idx]
        # end_diastole_time_idx = np.nonzero(self.post_nodefield.dict['time'] == self.simulation_dict['end_diastole_t'])
        # self.qoi['mean_lv_thickness_end_diastole'] = lv_wall_thickness_transient[end_diastole_time_idx]
        # self.qoi['mean_rv_thickness_end_diastole'] = rv_wall_thickness_transient[end_diastole_time_idx]
        # end_systole_time_idx = np.nonzero(self.post_nodefield.dict['lv_phase'] == 3)[0][0] # Index at which phase first changes to 3 IVR.
        # self.qoi['mean_lv_thickness_end_systole'] = lv_wall_thickness_transient[end_systole_time_idx]
        # self.qoi['mean_rv_thickness_end_systole'] = rv_wall_thickness_transient[end_systole_time_idx]
        #
        #
        # print('Fibre stretch ratio transient')
        # lv_lambda_transients = np.mean(self.post_nodefield.dict['LAMBD'][lv_nodes, 0, :], axis=0)
        # rv_lambda_transients = np.mean(self.post_nodefield.dict['LAMBD'][rv_nodes, 0, :], axis=0)
        # np.savetxt(self.results_dir+'lv_fibre_stretch_ratio_transient.txt', lv_lambda_transients, delimiter=',')
        # np.savetxt(self.results_dir + 'rv_fibre_stretch_ratio_transient.txt', rv_lambda_transients, delimiter=',')
        # plt.figure()
        # plt.plot(self.post_nodefield.dict['time'], lv_lambda_transients, self.post_nodefield.dict['time'], rv_lambda_transients)
        # plt.xlabel('Time (s)')
        # plt.ylabel('Fibre stretch ratio')
        # plt.legend(['LV', 'RV'])
        # plt.savefig(self.results_dir + 'fibre_stretch_ratio_transient.png')


    def evaluate_strain_biomarkers(self, beat):
        print('Strain evaluations in ventricular coordinates')
        self.read_binary_outputs(read_field_name='EPSXX', read_field_type='scalar') # to get time information
        # self.read_csv_fields(read_field_name='EPSXX', read_field_type='scalar') # to get time information
        # Select time segment according to specified beat # TODO assuming CL is always 1.0 s
        CL = 1.0
        time_idx = []
        strain_t = []
        for time_i in range(self.post_nodefield.dict['time'].shape[0]):
            if (self.post_nodefield.dict['time'][time_i] > (beat - 1) * CL) & (
                    self.post_nodefield.dict['time'][time_i] < beat * CL):
                strain_t.append(self.post_nodefield.dict['time'][time_i] - (beat - 1) * CL)
                time_idx.append(time_i)
        strain_t = np.array(strain_t)
        time_idx = np.array(time_idx, dtype=int)

        # Evaluate longitudinal, circumferential, and radial strains
        if 'E_cc' not in self.post_nodefield.dict.keys():
            # self.read_csv_fields(read_field_name='EPSYY', read_field_type='scalar')
            # self.read_csv_fields(read_field_name='EPSZZ', read_field_type='scalar')
            # self.read_csv_fields(read_field_name='EPSXY', read_field_type='scalar')
            # self.read_csv_fields(read_field_name='EPSXZ', read_field_type='scalar')
            # self.read_csv_fields(read_field_name='EPSYZ', read_field_type='scalar')
            self.read_binary_outputs(read_field_name='EPSYY', read_field_type='scalar')
            self.read_binary_outputs(read_field_name='EPSZZ', read_field_type='scalar')
            self.read_binary_outputs(read_field_name='EPSXY', read_field_type='scalar')
            self.read_binary_outputs(read_field_name='EPSXZ', read_field_type='scalar')
            self.read_binary_outputs(read_field_name='EPSYZ', read_field_type='scalar')
            print('Evaluating E_cc, E_ll, E_rr...')
            E = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0], 3, 3)) # Strain tensors at each node over time
            E_cc = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            E_ll = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            E_rr = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            exx = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            eyy = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            ezz = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            exy = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            exz = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            eyz = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
            local_l = pymp.shared.array((self.geometry.number_of_nodes, 3))
            local_c = pymp.shared.array((self.geometry.number_of_nodes, 3))
            local_r = pymp.shared.array((self.geometry.number_of_nodes, 3))
            local_l[:, :] = self.node_fields.dict['local_l']
            local_c[:, :] = self.node_fields.dict['local_c']
            local_r[:, :] = self.node_fields.dict['local_r']
            exx[:, :] = self.post_nodefield.dict['EPSXX'][:, time_idx]
            eyy[:, :] = self.post_nodefield.dict['EPSYY'][:, time_idx]
            ezz[:, :] = self.post_nodefield.dict['EPSZZ'][:, time_idx]
            exy[:, :] = self.post_nodefield.dict['EPSXY'][:, time_idx]
            exz[:, :] = self.post_nodefield.dict['EPSXZ'][:, time_idx]
            eyz[:, :] = self.post_nodefield.dict['EPSYZ'][:, time_idx]
            threadsNum = int(multiprocessing.cpu_count()*0.7)
            with pymp.Parallel(min(threadsNum, time_idx.shape[0])) as p1:
                for time_i in p1.range(time_idx.shape[0]):
                    E[:, time_i, 0, 0] = exx[:, time_i]
                    E[:, time_i, 0, 1] = exy[:, time_i]
                    E[:, time_i, 0, 2] = exz[:, time_i]
                    E[:, time_i, 1, 0] = exy[:, time_i]
                    E[:, time_i, 1, 1] = eyy[:, time_i]
                    E[:, time_i, 1, 2] = eyz[:, time_i]
                    E[:, time_i, 2, 0] = exz[:, time_i]
                    E[:, time_i, 2, 1] = eyz[:, time_i]
                    E[:, time_i, 2, 2] = ezz[:, time_i]
                    for node_i in range(self.geometry.number_of_nodes):
                        E_cc[node_i, time_i] = np.dot(np.dot(E[node_i, time_i, :, :], local_c[node_i, :]), local_c[node_i, :])
                        E_ll[node_i, time_i] = np.dot(np.dot(E[node_i, time_i, :, :], local_l[node_i, :]), local_l[node_i, :])
                        E_rr[node_i, time_i] = np.dot(np.dot(E[node_i, time_i, :, :], local_r[node_i, :]), local_r[node_i, :])
            # Delete raw inputs
            temp = self.post_nodefield.dict.pop('EPSXX')
            temp = self.post_nodefield.dict.pop('EPSYY')
            temp = self.post_nodefield.dict.pop('EPSZZ')
            temp = self.post_nodefield.dict.pop('EPSXY')
            temp = self.post_nodefield.dict.pop('EPSXZ')
            temp = self.post_nodefield.dict.pop('EPSYZ')
            # Save local strain values
            self.post_nodefield.add_field(data=E_cc, data_name='E_cc', field_type='postnodefield')
            self.post_nodefield.add_field(data=E_rr, data_name='E_rr', field_type='postnodefield')
            self.post_nodefield.add_field(data=E_ll, data_name='E_ll', field_type='postnodefield')
            # self.save_postprocessing_fields()
        else:
            print('E_ll, E_cc, E_rr have already been calculated')

        # Save strains in short and long axis slices for comparison with https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
        print('Extract strain transients in short and long axis slices')
        base_short_axis_nodes = np.nonzero(self.node_fields.dict['short-axis-slices'] == 3)[0]
        mid_short_axis_nodes = np.nonzero(self.node_fields.dict['short-axis-slices'] == 2)[0]
        apex_short_axis_nodes = np.nonzero(self.node_fields.dict['short-axis-slices'] == 1)[0]
        two_chamber_long_axis_nodes = np.nonzero(self.node_fields.dict['long-axis-slices'] == 1)[0]
        four_chamber_long_axis_nodes = np.nonzero(self.node_fields.dict['long-axis-slices'] == 2)[0]
        three_chamber_long_axis_nodes = np.nonzero(self.node_fields.dict['long-axis-slices'] == 3)[0]
        base_E_cc = pymp.shared.array((base_short_axis_nodes.shape[0], time_idx.shape[0], ))
        mid_E_cc = pymp.shared.array((mid_short_axis_nodes.shape[0], time_idx.shape[0]))
        apex_E_cc = pymp.shared.array((apex_short_axis_nodes.shape[0], time_idx.shape[0]))
        base_E_rr = pymp.shared.array((base_short_axis_nodes.shape[0], time_idx.shape[0]))
        mid_E_rr = pymp.shared.array((mid_short_axis_nodes.shape[0], time_idx.shape[0]))
        apex_E_rr = pymp.shared.array((apex_short_axis_nodes.shape[0], time_idx.shape[0]))
        two_chamber_E_ll = pymp.shared.array((two_chamber_long_axis_nodes.shape[0], time_idx.shape[0]))
        four_chamber_E_ll = pymp.shared.array((four_chamber_long_axis_nodes.shape[0], time_idx.shape[0], ))
        three_chamber_E_ll = pymp.shared.array((three_chamber_long_axis_nodes.shape[0], time_idx.shape[0]))
        E_cc = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
        E_ll = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
        E_rr = pymp.shared.array((self.geometry.number_of_nodes, time_idx.shape[0]))
        E_cc[:,:] = self.post_nodefield.dict['E_cc']
        E_ll[:, :] = self.post_nodefield.dict['E_ll']
        E_rr[:, :] = self.post_nodefield.dict['E_rr']
        threadsNum = int(multiprocessing.cpu_count()*0.7)
        with pymp.Parallel(min(threadsNum, time_idx.shape[0])) as p1:
            for time_i in p1.range(time_idx.shape[0]):
                # Circumferential strains
                base_E_cc[:, time_i] = E_cc[base_short_axis_nodes, time_i]
                mid_E_cc[:, time_i] = E_cc[mid_short_axis_nodes, time_i]
                apex_E_cc[:, time_i] = E_cc[apex_short_axis_nodes, time_i]
                # Radial strains
                base_E_rr[:, time_i] = E_rr[base_short_axis_nodes, time_i]
                mid_E_rr[:, time_i] = E_rr[mid_short_axis_nodes, time_i]
                apex_E_rr[:, time_i] = E_rr[apex_short_axis_nodes, time_i]
                # Longitudinal strains
                two_chamber_E_ll[:, time_i] = E_ll[two_chamber_long_axis_nodes, time_i]
                four_chamber_E_ll[:, time_i] = E_ll[four_chamber_long_axis_nodes, time_i]
                three_chamber_E_ll[:, time_i] = E_ll[three_chamber_long_axis_nodes, time_i]

        # Evaluate median and upper and lower quartiles of traces
        print('Evaluate median and interquartile ranges of strain traces')
        # Circumferential
        systole = np.min(base_E_cc, axis=1)
        tol = 0.01 * abs(np.min(systole))
        # base_E_cc_median = base_E_cc[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # base_E_cc_uq = base_E_cc[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # base_E_cc_lq = base_E_cc[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # systole = np.min(mid_E_cc, axis=1)
        # mid_E_cc_median = mid_E_cc[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # mid_E_cc_uq = mid_E_cc[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # mid_E_cc_lq = mid_E_cc[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # systole = np.min(apex_E_cc, axis=1)
        # apex_E_cc_median = apex_E_cc[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # apex_E_cc_uq = apex_E_cc[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # apex_E_cc_lq = apex_E_cc[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # # Radial
        # systole = np.min(base_E_rr, axis=1)
        # base_E_rr_median = base_E_rr[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # base_E_rr_uq = base_E_rr[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # base_E_rr_lq = base_E_rr[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # systole = np.min(mid_E_rr, axis=1)
        # mid_E_rr_median = mid_E_rr[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # mid_E_rr_uq = mid_E_rr[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # mid_E_rr_lq = mid_E_rr[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # systole = np.min(apex_E_rr, axis=1)
        # apex_E_rr_median = apex_E_rr[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # apex_E_rr_uq = apex_E_rr[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # apex_E_rr_lq = apex_E_rr[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # # Longitudinal
        # systole = np.min(two_chamber_E_ll, axis=1)
        # two_chamber_E_ll_median = two_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # two_chamber_E_ll_uq = two_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # two_chamber_E_ll_lq = two_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # systole = np.min(four_chamber_E_ll, axis=1)
        # four_chamber_E_ll_median = four_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # four_chamber_E_ll_uq = four_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # four_chamber_E_ll_lq = four_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]
        # systole = np.min(three_chamber_E_ll, axis=1)
        # three_chamber_E_ll_median = three_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 50)) < tol)[0][0], :]
        # three_chamber_E_ll_uq = three_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 75)) < tol)[0][0], :]
        # three_chamber_E_ll_lq = three_chamber_E_ll[np.where(abs(systole-np.percentile(systole, 25)) < tol)[0][0], :]

        self.strain_transients = {}
        mean_mid_E_cc = np.mean(mid_E_cc, axis=0)
        mean_mid_E_rr = np.mean(mid_E_rr, axis=0)
        mean_four_chamber_E_ll = np.mean(four_chamber_E_ll, axis=0)
        self.strain_transients['strain_t'] = strain_t
        self.strain_transients['mean_mid_E_cc'] = mean_mid_E_cc
        self.strain_transients['mean_mid_E_rr'] = mean_mid_E_rr
        self.strain_transients['mean_four_chamber_E_ll'] = mean_four_chamber_E_ll

        qoi = {}
        qoi['max_mid_Ecc'] = np.max(mean_mid_E_cc)
        qoi['min_mid_Ecc'] = np.min(mean_mid_E_cc)
        qoi['max_mid_Err'] = np.max(mean_mid_E_rr)
        qoi['min_mid_Err'] = np.min(mean_mid_E_rr)
        qoi['max_four_chamber_Ell'] = np.max(mean_four_chamber_E_ll)
        qoi['min_four_chamber_Ell'] = np.min(mean_four_chamber_E_ll)
        self.qoi.update(qoi)

        # self.strain_transients['base_E_cc'] = base_E_cc
        # self.strain_transients['base_E_cc_median'] = base_E_cc_median
        # self.strain_transients['base_E_cc_uq'] = base_E_cc_uq
        # self.strain_transients['base_E_cc_lq'] = base_E_cc_lq
        #
        # self.strain_transients['mid_E_cc'] = mid_E_cc
        # self.strain_transients['mean_mid_Ecc'] = np.mean(mid_E_cc, axis=0)
        # self.strain_transients['mid_E_cc_median'] = mid_E_cc_median
        # self.strain_transients['mid_E_cc_uq'] = mid_E_cc_uq
        # self.strain_transients['mid_E_cc_lq'] = mid_E_cc_lq
        #
        # self.strain_transients['apex_E_cc'] = apex_E_cc
        # self.strain_transients['apex_E_cc_median'] = apex_E_cc_median
        # self.strain_transients['apex_E_cc_uq'] = apex_E_cc_uq
        # self.strain_transients['apex_E_cc_lq'] = apex_E_cc_lq
        #
        # self.strain_transients['base_E_rr'] = base_E_rr
        # self.strain_transients['base_E_rr_median'] = base_E_rr_median
        # self.strain_transients['base_E_rr_uq'] = base_E_rr_uq
        # self.strain_transients['base_E_rr_lq'] = base_E_rr_lq
        #
        # self.strain_transients['mid_E_rr'] = mid_E_rr
        # self.strain_transients['mean_mid_Err'] = np.mean(mid_E_rr, axis=0)
        # self.strain_transients['mid_E_rr_median'] = mid_E_rr_median
        # self.strain_transients['mid_E_rr_uq'] = mid_E_rr_uq
        # self.strain_transients['mid_E_rr_lq'] = mid_E_rr_lq
        #
        # self.strain_transients['apex_E_rr'] = apex_E_rr
        # self.strain_transients['apex_E_rr_median'] = apex_E_rr_median
        # self.strain_transients['apex_E_rr_uq'] = apex_E_rr_uq
        # self.strain_transients['apex_E_rr_lq'] = apex_E_rr_lq
        #
        # self.strain_transients['two_chamber_E_ll'] = two_chamber_E_ll
        # self.strain_transients['two_chamber_E_ll_median'] = two_chamber_E_ll_median
        # self.strain_transients['two_chamber_E_ll_uq'] = two_chamber_E_ll_uq
        # self.strain_transients['two_chamber_E_ll_lq'] = two_chamber_E_ll_lq
        #
        # self.strain_transients['four_chamber_E_ll'] = four_chamber_E_ll
        # self.strain_transients['mean_four_chamber_E_ll'] = np.mean(four_chamber_E_ll, axis=0)
        # self.strain_transients['four_chamber_E_ll_median'] = four_chamber_E_ll_median
        # self.strain_transients['four_chamber_E_ll_uq'] = four_chamber_E_ll_uq
        # self.strain_transients['four_chamber_E_ll_lq'] = four_chamber_E_ll_lq
        #
        # self.strain_transients['three_chamber_E_ll'] = three_chamber_E_ll
        # self.strain_transients['three_chamber_E_ll_median'] = three_chamber_E_ll_median
        # self.strain_transients['three_chamber_E_ll_uq'] = three_chamber_E_ll_uq
        # self.strain_transients['three_chamber_E_ll_lq'] = three_chamber_E_ll_lq


        # del base_E_cc
        # del mid_E_cc
        # del apex_E_cc
        # del base_E_rr
        # del mid_E_rr
        # del apex_E_rr
        # del two_chamber_E_ll
        # del four_chamber_E_ll
        # del three_chamber_E_ll


    def evaluate_baseline_qoi_against_healthy_ranges(self, qoi_names):
        simulated_qois = self.qoi
        self.baseline_qoi_differences = {}
        for qoi_i in range(len(qoi_names)):
            qoi_name = qoi_names[qoi_i]
            if (simulated_qois[qoi_name] <= self.healthy_ranges[qoi_name][1]) & (
                    simulated_qois[qoi_name] >= self.healthy_ranges[qoi_name][0]):
                self.baseline_qoi_differences[qoi_name] = 0.0 # Within range
            elif (simulated_qois[qoi_name] > self.healthy_ranges[qoi_name][1]):
                self.baseline_qoi_differences[qoi_name] = self.healthy_ranges[qoi_name][1] - simulated_qois[qoi_name]
            elif (simulated_qois[qoi_name] < self.healthy_ranges[qoi_name][0]):
                self.baseline_qoi_differences[qoi_name] = self.healthy_ranges[qoi_name][0] - simulated_qois[qoi_name]
        print(self.baseline_qoi_differences)

    def shift_to_start_at_ED(self, t, trace):
        t_tol = 1e-3
        ed_idx = np.argmin(abs(t-self.simulation_dict['end_diastole_t'][0]))
        return np.concatenate((trace[ed_idx:], trace[:ed_idx]))


    def get_first_derivative(self, t, y):
        dydt = [0]
        for i in range(len(y)-1):
            i += 1
            dydt.append( (y[i] - y[i - 1]) / (t[i] - t[i - 1]))
        assert len(dydt) == len(y)
        return np.array(dydt)

    def visualise_strain_volume_loops(self, beat):
        fig = plt.figure(tight_layout=True, figsize=(18, 10))
        gs = GridSpec(1, 1)
        ax = fig.add_subplot(gs[0, 0])
        ax.set_title('Longitudinal strain vs LV volume')

        # Get Ell
        Ell = self.strains['four_chamber_E_ll_median']
        Ell_t = self.strains['strain_t']

        # Get LVV
        LVV = self.pvs['vls'][beat - 1]
        LVV_t = self.pvs['ts'][beat - 1]

        Ell_resampled = np.zeros(LVV.shape[0])
        for i in range(0, len(LVV_t)):
            map_idx = np.argmin(abs(Ell_t - LVV_t[i]))
            print(map_idx)
            Ell_resampled[i] = Ell[map_idx]

        ax.plot(LVV, Ell_resampled)
        plt.show()

    def visualise_qoi_comparisons(self, qoi_names, simulated_qois=None):

        qoi_units = ['ms', 'ms', 'mV', 'mL', 'mL', 'kPa', '%', 'mL', 'mL/s', 'mL/s', 'kPa/s', 'mL', 'mL',
                     'kPa', 'mL',  'cm', 'cm', 'cm', '', '', '', '']
        simulated_qois = self.qoi
        fig = plt.figure(tight_layout=True, figsize=(5, 10))
        gs = GridSpec(len(qoi_names), 1)
        for qoi_i in range(len(qoi_names)):
            qoi_name = qoi_names[qoi_i]
            qoi_unit = qoi_units[qoi_i]
            ax = fig.add_subplot(gs[qoi_i, 0])
            if (simulated_qois[qoi_name] <= self.healthy_ranges[qoi_name][1]) & (
                    simulated_qois[qoi_name] >= self.healthy_ranges[qoi_name][0]):
                ax.barh(0, simulated_qois[qoi_name], align='center', color='g')
            elif (simulated_qois[qoi_name] > self.healthy_ranges[qoi_name][1]):
                ax.barh(0, simulated_qois[qoi_name], align='center', color='r')
            elif (simulated_qois[qoi_name] < self.healthy_ranges[qoi_name][0]):
                ax.barh(0, simulated_qois[qoi_name], align='center', color='b')
            ax.set_yticks([0], [qoi_name])
            ax.set_xticks([])
            plt.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on', labelbottom='off')
            plt.text(self.healthy_ranges[qoi_name][1], 0, qoi_unit, ha='right', va='center')
            ax.axvspan(self.healthy_ranges[qoi_name][0], self.healthy_ranges[qoi_name][1], alpha=0.3, color='green')
            for spine in plt.gca().spines.values():
                spine.set_visible(False)
        plt.show()

    def visualise_calibration_comparisons_global(self, beat):
        fig = plt.figure(tight_layout=True, figsize=(8, 12))
        gs = GridSpec(5, 6)
        ax_lv_pv = fig.add_subplot(gs[1:3, 0:3])
        ax_lv_pv.set_title('PV loops')
        # ax_rv_pv = fig.add_subplot(gs[0:2, 2:4])
        # ax_rv_pv.set_title('Right PV')
        # ax_rv_pv.set_ylim([0, 5])
        ax_vt = fig.add_subplot(gs[1, 3:6])
        ax_vt.set_title('Volume transients')
        ax_pt = fig.add_subplot(gs[2, 3:6])
        ax_pt.set_title('Pressure transients')
        ax_thickness = fig.add_subplot(gs[3, 0:2])
        ax_thickness.set_title('Midventricular Wall thickness')
        ax_torsion = fig.add_subplot(gs[3, 2:4])
        ax_torsion.set_title('Torsion')
        ax_avpd = fig.add_subplot(gs[3, 4:6])
        ax_avpd.set_title('AVPD')
        ax_avpd.set_xlim(0, 1)
        ax_lambda = fig.add_subplot(gs[4, 0:2])
        ax_lambda.set_title('Fibre stretch ratio')
        ax_ell = fig.add_subplot(gs[4, 2:4])
        ax_ell.set_title('Longitudinal strain')
        ax_err = fig.add_subplot(gs[4, 4:6])
        ax_err.set_title('Radial strain')

        ax_V1 = fig.add_subplot(gs[0, 0])
        ax_V1.set_title('V1')
        ax_V2 = fig.add_subplot(gs[0, 1])
        ax_V2.set_title('V2')
        ax_V3 = fig.add_subplot(gs[0, 2])
        ax_V3.set_title('V3')
        ax_V4 = fig.add_subplot(gs[0, 3])
        ax_V4.set_title('V4')
        ax_V5 = fig.add_subplot(gs[0, 4])
        ax_V5.set_title('V5')
        ax_V6 = fig.add_subplot(gs[0, 5])
        ax_V6.set_title('V6')

        # ECGs
        ax_V1.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V1s'][beat - 1] / self.ecgs['max_all_leads'])
        ax_V1.axvspan(self.healthy_ranges['QT'][0], self.healthy_ranges['QT'][1], alpha=0.3, color='green')
        ax_V1.axvspan(self.healthy_ranges['QRS_duration'][0], self.healthy_ranges['QRS_duration'][1], alpha=0.3,
                      color='green')
        ax_V2.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V2s'][beat - 1] / self.ecgs['max_all_leads'])
        ax_V3.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V3s'][beat - 1] / self.ecgs['max_all_leads'])
        ax_V4.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V4s'][beat - 1] / self.ecgs['max_all_leads'])
        ax_V5.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V5s'][beat - 1] / self.ecgs['max_all_leads'])
        ax_V6.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V6s'][beat - 1] / self.ecgs['max_all_leads'])
        ax_V2.axvspan(self.healthy_ranges['QT'][0], self.healthy_ranges['QT'][1], alpha=0.3, color='green')
        ax_V2.axvspan(self.healthy_ranges['QRS_duration'][0], self.healthy_ranges['QRS_duration'][1], alpha=0.3,
                      color='green')
        ax_V3.axvspan(self.healthy_ranges['QT'][0], self.healthy_ranges['QT'][1], alpha=0.3, color='green')
        ax_V3.axvspan(self.healthy_ranges['QRS_duration'][0], self.healthy_ranges['QRS_duration'][1], alpha=0.3,
                      color='green')
        ax_V4.axvspan(self.healthy_ranges['QT'][0], self.healthy_ranges['QT'][1], alpha=0.3, color='green')
        ax_V4.axvspan(self.healthy_ranges['QRS_duration'][0], self.healthy_ranges['QRS_duration'][1], alpha=0.3,
                      color='green')
        ax_V5.axvspan(self.healthy_ranges['QT'][0], self.healthy_ranges['QT'][1], alpha=0.3, color='green')
        ax_V5.axvspan(self.healthy_ranges['QRS_duration'][0], self.healthy_ranges['QRS_duration'][1], alpha=0.3,
                      color='green')
        ax_V6.axvspan(self.healthy_ranges['QT'][0], self.healthy_ranges['QT'][1], alpha=0.3, color='green')
        ax_V6.axvspan(self.healthy_ranges['QRS_duration'][0], self.healthy_ranges['QRS_duration'][1], alpha=0.3,
                      color='green')
        ax_V1.set_xlim([0, 500])
        ax_V2.set_xlim([0, 500])
        ax_V3.set_xlim([0, 500])
        ax_V4.set_xlim([0, 500])
        ax_V5.set_xlim([0, 500])
        ax_V6.set_xlim([0, 500])
        ax_V1.set_ylim([-1, 1])
        ax_V2.set_ylim([-1, 1])
        ax_V3.set_ylim([-1, 1])
        ax_V4.set_ylim([-1, 1])
        ax_V5.set_ylim([-1, 1])
        ax_V6.set_ylim([-1, 1])

        # Pressure volume loops
        ax_lv_pv.plot(self.pvs['vls'][beat - 1], self.pvs['pls'][beat - 1] / 10000)
        ax_lv_pv.axvspan(self.healthy_ranges['LVEDV'][0], self.healthy_ranges['LVEDV'][1], alpha=0.3, color='C0')
        ax_lv_pv.axvspan(self.healthy_ranges['LVESV'][0], self.healthy_ranges['LVESV'][1], alpha=0.3, color='C0')
        ax_lv_pv.axhspan(self.healthy_ranges['LVESP'][0], self.healthy_ranges['LVESP'][1], alpha=0.3, color='C0')
        ax_lv_pv.axhspan(self.healthy_ranges['LVEDP'][0], self.healthy_ranges['LVEDP'][1], alpha=0.3, color='C0')


        ax_lv_pv.plot(self.pvs['vrs'][beat - 1], self.pvs['prs'][beat - 1] / 10000)
        ax_lv_pv.axvspan(self.healthy_ranges['RVEDV'][0], self.healthy_ranges['RVEDV'][1], alpha=0.3, color='C1')
        ax_lv_pv.axvspan(self.healthy_ranges['RVESV'][0], self.healthy_ranges['RVESV'][1], alpha=0.3, color='C1')
        ax_lv_pv.axhspan(self.healthy_ranges['RVESP'][0], self.healthy_ranges['RVESP'][1], alpha=0.3, color='C1')
        # ax_lv_pv.set_title('RVEF: ' + str((np.amax(self.pvs['vrs'][beat - 1]) - np.amin(self.pvs['vrs'][beat - 1])) /
        #                                   np.amax(self.pvs['vrs'][beat - 1]) * 100) + '%')
        ax_lv_pv.set_title('LVEF: {lvef:d} %, RVEF: {rvef:d}'.format(
            lvef=int((np.amax(self.pvs['vls'][beat - 1]) - np.amin(self.pvs['vls'][beat - 1])) /
                 np.amax(self.pvs['vls'][beat - 1]) * 100),
                           rvef=int((np.amax(self.pvs['vrs'][beat - 1]) - np.amin(self.pvs['vrs'][beat - 1])) /
                                np.amax(self.pvs['vrs'][beat - 1]) * 100)))
        ax_lv_pv.set_xlabel('Volume (mL)')
        ax_lv_pv.set_ylabel('Pressure (kPa)')

        # Pressure transients
        ax_pt.plot(self.pvs['ts'][beat-1], self.pvs['pls'][beat-1]/10000, color='C0', label='LV')
        ax_pt.plot(self.pvs['ts'][beat-1], self.pvs['prs'][beat-1]/10000, color='C1', label='RV')
        ax_pt.legend()
        t = self.pvs['ts'][beat - 1]
        p = self.pvs['pls'][beat - 1] / 10000
        dpdt = self.get_first_derivative(t=t, y=p)
        peak = np.argmax(dpdt)
        peak_dpdt = np.amax(dpdt)
        ax_pt.axline((t[peak], p[peak]), slope=peak_dpdt, color='grey', linestyle='--')
        intercept_0 = p[peak] - self.healthy_ranges['dpdt_max'][0] * t[peak]
        intercept_1 = p[peak] - self.healthy_ranges['dpdt_max'][1] * t[peak]
        ax_pt.fill_between(t, self.healthy_ranges['dpdt_max'][0] * t + intercept_0,
                           self.healthy_ranges['dpdt_max'][1] * t + intercept_1, alpha=0.3, facecolor='green')
        ax_pt.set_ylim(0, 17)

        # Volume transients
        # Diastolic fill
        end_systole_idx = np.argmin(self.pvs['vls'][beat-1])
        ldvdt = self.get_first_derivative(t=self.pvs['ts'][beat-1], y=self.pvs['vls'][beat-1])
        dvdt_ejection = np.amin(ldvdt[10:end_systole_idx])
        idx_ejection = np.argmin(ldvdt[10:end_systole_idx]) + 10
        dvdt_filling = np.amax(ldvdt[(end_systole_idx+400):])
        idx_filling = np.argmax(ldvdt[(end_systole_idx+400):]) + end_systole_idx + 400
        landmarks = np.zeros(len(self.pvs['ts'][beat-1]))
        landmarks[idx_ejection] = 1
        landmarks[idx_filling] = 2
        landmarks_shifted = self.shift_to_start_at_ED(self.pvs['ts'][beat-1], landmarks)
        t = self.pvs['ts'][beat - 1]
        vl = self.pvs['vls'][beat-1]
        vr = self.pvs['vrs'][beat-1]
        # vl = self.shift_to_start_at_ED(self.pvs['ts'][beat-1], self.pvs['vls'][beat-1])
        # vr = self.shift_to_start_at_ED(self.pvs['ts'][beat-1], self.pvs['vrs'][beat-1])
        ax_vt.plot(self.pvs['ts'][beat-1], vl, label='LV', color='C0')
        ax_vt.plot(self.pvs['ts'][beat-1], vr, label='RV', color='C1')
        ax_vt.legend()
        # ax_vt.axhspan(self.healthy_ranges['LVEDV'][0], self.healthy_ranges['LVEDV'][1], alpha=0.3, color='C0')
        # ax_vt.axhspan(self.healthy_ranges['LVESV'][0], self.healthy_ranges['LVESV'][1], alpha=0.3, color='C0')
        # ax_vt.axhspan(self.healthy_ranges['RVEDV'][0], self.healthy_ranges['RVEDV'][1], alpha=0.3, color='C1')
        # ax_vt.axhspan(self.healthy_ranges['RVESV'][0], self.healthy_ranges['RVESV'][1], alpha=0.3, color='C1')
        ax_vt.axline((t[np.where(landmarks==1)[0][0]],
                      vl[np.where(landmarks==1)[0][0]]),
                     slope=dvdt_ejection, color='grey', linestyle='--')
        intercept_0 = vl[np.where(landmarks==1)[0][0]] +  self.healthy_ranges['dvdt_ejection'][0] * t[np.where(landmarks==1)[0][0]]
        intercept_1 = vl[np.where(landmarks==1)[0][0]] +  self.healthy_ranges['dvdt_ejection'][1] * t[np.where(landmarks==1)[0][0]]
        ax_vt.fill_between(t, -self.healthy_ranges['dvdt_ejection'][0] * t + intercept_0, -self.healthy_ranges['dvdt_ejection'][1] * t + intercept_1, alpha=0.3, facecolor='green')
        # Diastolic filling rate
        ax_vt.axline((t[np.where(landmarks==2)[0][0]],
                      vl[np.where(landmarks==2)[0][0]]),
                     slope=dvdt_filling, color='grey', linestyle='--')
        intercept_0 = vl[np.where(landmarks==2)[0][0]] - self.healthy_ranges['dvdt_filling'][0] * t[np.where(landmarks==2)[0][0]]
        intercept_1 = vl[np.where(landmarks==2)[0][0]] - self.healthy_ranges['dvdt_filling'][1] * t[np.where(landmarks==2)[0][0]]
        ax_vt.fill_between(t, self.healthy_ranges['dvdt_filling'][0] * t + intercept_0,
                           self.healthy_ranges['dvdt_filling'][1] * t + intercept_1, alpha=0.3, facecolor='green')
        ax_vt.set_ylim(self.healthy_ranges['LVESV'][0], self.healthy_ranges['LVEDV'][1])
        ax_vt.set_xlim(0, 1)

        # Displacements
        avpd = self.shift_to_start_at_ED(self.deformation_transients['deformation_t'], self.deformation_transients['avpd'])
        apical_displacement = self.shift_to_start_at_ED(self.deformation_transients['deformation_t'],
                                         self.deformation_transients['apical_displacement'])
        ax_avpd.plot(self.deformation_transients['deformation_t'], avpd, label='AVPD',
                     color='C0')
        # ax_avpd.plot(self.deformation_transients['deformation_t'], apical_displacement, label='Apex displacement', color='C0')
        ax_avpd.axhspan(np.amax(avpd)-self.healthy_ranges['AVPD'][0],np.amax(avpd)-self.healthy_ranges['AVPD'][1], alpha=0.3, color='C0')
        # ax_avpd.axhspan(self.healthy_ranges['apical_displacement'][0],
        #                 self.healthy_ranges['apical_displacement'][1], alpha=0.3,
        #                 color='C0')

        thickness = self.shift_to_start_at_ED(self.deformation_transients['deformation_t'], self.deformation_transients['lv_wall_thickness'])
        ax_thickness.plot(self.deformation_transients['deformation_t'], thickness, label='AVPD',
                     color='C0')
        ax_thickness.axhspan(self.healthy_ranges['ED_wall_thickness'][0], self.healthy_ranges['ED_wall_thickness'][1],
                             alpha=0.3, color='green')
        ax_thickness.axhspan(self.healthy_ranges['ES_wall_thickness'][0], self.healthy_ranges['ES_wall_thickness'][1],
                             alpha=0.3, color='green')


        # strain_transients
        ell_median = self.shift_to_start_at_ED(self.strain_transients['strain_t'], self.strain_transients['four_chamber_E_ll_median'])
        ell_uq = self.shift_to_start_at_ED(self.strain_transients['strain_t'], self.strain_transients['four_chamber_E_ll_uq'])
        ell_lq = self.shift_to_start_at_ED(self.strain_transients['strain_t'], self.strain_transients['four_chamber_E_ll_lq'])
        ax_ell.plot(self.strain_transients['strain_t'], ell_median, label='Median', color='C0', linestyle='-')
        ax_ell.plot(self.strain_transients['strain_t'], ell_uq, label='IQR', color='C0', linestyle='--')
        ax_ell.plot(self.strain_transients['strain_t'], ell_lq, color='C0', linestyle='--')
        ax_ell.axhspan(self.healthy_ranges['peak_E_ll'][0], self.healthy_ranges['peak_E_ll'][1], alpha=0.3, color='green')
        ax_ell.legend()

        err_median = self.shift_to_start_at_ED(self.strain_transients['strain_t'], self.strain_transients['mid_E_rr_median'])
        err_uq = self.shift_to_start_at_ED(self.strain_transients['strain_t'], self.strain_transients['mid_E_rr_uq'])
        err_lq = self.shift_to_start_at_ED(self.strain_transients['strain_t'], self.strain_transients['mid_E_rr_lq'])
        ax_err.plot(self.strain_transients['strain_t'], err_median, color='C0', linestyle='-')
        ax_err.plot(self.strain_transients['strain_t'], err_uq, color='C0', linestyle='--')
        ax_err.plot(self.strain_transients['strain_t'], err_lq, color='C0', linestyle='--')
        ax_err.axhspan(self.healthy_ranges['peak_E_rr'][0], self.healthy_ranges['peak_E_rr'][1], alpha=0.3, color='green')

        lambda_endo = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'], self.fibre_work['endo_midshort_lambda_median'])
        lambda_endo_uq = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'], self.fibre_work['endo_midshort_lambda_uq'])
        lambda_endo_lq = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'],
                                                   self.fibre_work['endo_midshort_lambda_lq'])
        lambda_epi = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'], self.fibre_work['epi_midshort_lambda_median'])
        lambda_epi_uq = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'],
                                                   self.fibre_work['epi_midshort_lambda_uq'])
        lambda_epi_lq = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'],
                                                   self.fibre_work['epi_midshort_lambda_lq'])

        lambda_mid = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'], self.fibre_work['mid_midshort_lambda_median'])
        lambda_mid_uq = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'],
                                                  self.fibre_work['mid_midshort_lambda_uq'])
        lambda_mid_lq = self.shift_to_start_at_ED(self.fibre_work['fibrework_t'],
                                                  self.fibre_work['mid_midshort_lambda_lq'])

        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_endo, label='Endo', color='C2')
        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_endo_uq, color='C2', linestyle='--')
        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_endo_lq, color='C2', linestyle='--')

        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_epi, label='Epi', color='C3')
        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_epi_uq, color='C2', linestyle='--')
        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_epi_lq, color='C2', linestyle='--')

        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_mid, label='Mid', color='C4')
        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_mid_uq, color='C2', linestyle='--')
        ax_lambda.plot(self.fibre_work['fibrework_t'], lambda_mid_lq, color='C2', linestyle='--')

        ax_lambda.axhspan(self.healthy_ranges['lambda_endo'][0], self.healthy_ranges['lambda_endo'][1], alpha=0.3, color='C2')
        ax_lambda.axhspan(self.healthy_ranges['lambda_epi'][0], self.healthy_ranges['lambda_endo'][1], alpha=0.3,
                          color='C3')
        ax_lambda.axhspan(self.healthy_ranges['lambda_mid'][0], self.healthy_ranges['lambda_endo'][1], alpha=0.3,
                          color='C4')
        ax_lambda.legend()
        plt.show()


    def visualise_calibration_comparisons_strain(self):
        fig = plt.figure(tight_layout=True, figsize=(10, 12))
        gs = GridSpec(4, 3)

        ax_E_ff_endo = fig.add_subplot(gs[0, 0])
        ax_E_ff_endo.set_title('E_ff endo')
        ax_E_ff_endo.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_ff_mid = fig.add_subplot(gs[0, 1])
        ax_E_ff_mid.set_title('E_ff mid')
        ax_E_ff_mid.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_ff_epi = fig.add_subplot(gs[0, 2])
        ax_E_ff_epi.set_title('E_ff epi')
        ax_E_ff_epi.axhline(y=-0.15, color='red', linestyle='--')

        ax_E_cc_base = fig.add_subplot(gs[1, 0])
        ax_E_cc_base.set_title('E_cc Base')
        ax_E_cc_base.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_cc_mid = fig.add_subplot(gs[2, 0])
        ax_E_cc_mid.set_title('E_cc Mid')
        ax_E_cc_mid.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_cc_apex = fig.add_subplot(gs[3, 0])
        ax_E_cc_apex.set_title('E_cc Apex')
        ax_E_cc_apex.axhline(y=-0.15, color='red', linestyle='--')

        ax_E_rr_base = fig.add_subplot(gs[1, 1])
        ax_E_rr_base.set_title('E_rr Base')
        ax_E_rr_base.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_rr_mid = fig.add_subplot(gs[2, 1])
        ax_E_rr_mid.set_title('E_rr Mid')
        ax_E_rr_mid.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_rr_apex = fig.add_subplot(gs[3, 1])
        ax_E_rr_apex.set_title('E_rr Apex')
        ax_E_rr_apex.axhline(y=-0.15, color='red', linestyle='--')

        ax_E_ll_two = fig.add_subplot(gs[1, 2])
        ax_E_ll_two.set_title('E_ll 2-chamber')
        ax_E_ll_two.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_ll_four = fig.add_subplot(gs[2, 2])
        ax_E_ll_four.set_title('E_ll 4-chamber')
        ax_E_ll_four.axhline(y=-0.15, color='red', linestyle='--')
        ax_E_ll_three = fig.add_subplot(gs[3, 2])
        ax_E_ll_three.set_title('E_ll 3-chamber')
        ax_E_ll_three.axhline(y=-0.15, color='red', linestyle='--')

        # Fibre strain
        endo_E_ff = 0.5*(self.fibre_work['endo_midshort_lambda']**2-1)
        endo_E_ff_median = 0.5*(self.fibre_work['endo_midshort_lambda_median']**2 - 1)
        endo_E_ff_lq = 0.5 * (self.fibre_work['endo_midshort_lambda_lq'] ** 2 - 1)
        endo_E_ff_uq = 0.5 * (self.fibre_work['endo_midshort_lambda_uq'] ** 2 - 1)
        for node_i in range(self.fibre_work['endo_midshort_lambda'].shape[0]):
            ax_E_ff_endo.plot(self.fibre_work['fibrework_t'], endo_E_ff[node_i, :], color='m', alpha=0.3, linewidth=0.1)
        ax_E_ff_endo.plot(self.fibre_work['fibrework_t'], endo_E_ff_median, color='k', linestyle='--')
        ax_E_ff_endo.plot(self.fibre_work['fibrework_t'], endo_E_ff_lq, color='k', linestyle='-')
        ax_E_ff_endo.plot(self.fibre_work['fibrework_t'], endo_E_ff_uq, color='k', linestyle='-')

        mid_E_ff = 0.5 * (self.fibre_work['mid_midshort_lambda'] ** 2 - 1)
        mid_E_ff_median = 0.5 * (self.fibre_work['mid_midshort_lambda_median'] ** 2 - 1)
        mid_E_ff_lq = 0.5 * (self.fibre_work['mid_midshort_lambda_lq'] ** 2 - 1)
        mid_E_ff_uq = 0.5 * (self.fibre_work['mid_midshort_lambda_uq'] ** 2 - 1)
        for node_i in range(self.fibre_work['mid_midshort_lambda'].shape[0]):
            ax_E_ff_mid.plot(self.fibre_work['fibrework_t'], mid_E_ff[node_i, :], color='g', alpha=0.3, linewidth=0.1)
        ax_E_ff_mid.plot(self.fibre_work['fibrework_t'], mid_E_ff_median, color='k', linestyle='--')
        ax_E_ff_mid.plot(self.fibre_work['fibrework_t'], mid_E_ff_lq, color='k', linestyle='-')
        ax_E_ff_mid.plot(self.fibre_work['fibrework_t'], mid_E_ff_uq, color='k', linestyle='-')

        epi_E_ff = 0.5 * (self.fibre_work['epi_midshort_lambda'] ** 2 - 1)
        epi_E_ff_median = 0.5 * (self.fibre_work['epi_midshort_lambda_median'] ** 2 - 1)
        epi_E_ff_lq = 0.5 * (self.fibre_work['epi_midshort_lambda_lq'] ** 2 - 1)
        epi_E_ff_uq = 0.5 * (self.fibre_work['epi_midshort_lambda_uq'] ** 2 - 1)
        for node_i in range(self.fibre_work['epi_midshort_lambda'].shape[0]):
            ax_E_ff_epi.plot(self.fibre_work['fibrework_t'], epi_E_ff[node_i, :], color='b', alpha=0.3, linewidth=0.1)
        ax_E_ff_epi.plot(self.fibre_work['fibrework_t'], epi_E_ff_median, color='k', linestyle='--')
        ax_E_ff_epi.plot(self.fibre_work['fibrework_t'], epi_E_ff_lq, color='k', linestyle='-')
        ax_E_ff_epi.plot(self.fibre_work['fibrework_t'], epi_E_ff_uq, color='k', linestyle='-')

        # Strains
        for node_i in range(self.strain_transients['base_E_cc'].shape[0]):
            ax_E_cc_base.plot(self.strain_transients['strain_t'], self.strain_transients['base_E_cc'][node_i, :], color='m', alpha=0.3, linewidth=0.1)
        ax_E_cc_base.plot(self.strain_transients['strain_t'],self.strain_transients['base_E_cc_median'], color='k', linestyle='--')
        ax_E_cc_base.plot(self.strain_transients['strain_t'], self.strain_transients['base_E_cc_uq'], color='k', linestyle='-')
        ax_E_cc_base.plot(self.strain_transients['strain_t'], self.strain_transients['base_E_cc_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['mid_E_cc'].shape[0]):
            ax_E_cc_mid.plot(self.strain_transients['strain_t'], self.strain_transients['mid_E_cc'][node_i, :], color='m', alpha=0.3, linewidth=0.1)
        ax_E_cc_mid.plot(self.strain_transients['strain_t'],self.strain_transients['mid_E_cc_median'], color='k', linestyle='--')
        ax_E_cc_mid.plot(self.strain_transients['strain_t'], self.strain_transients['mid_E_cc_uq'], color='k', linestyle='-')
        ax_E_cc_mid.plot(self.strain_transients['strain_t'], self.strain_transients['mid_E_cc_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['apex_E_cc'].shape[0]):
            ax_E_cc_apex.plot(self.strain_transients['strain_t'], self.strain_transients['apex_E_cc'][node_i, :], color='m', alpha=0.3, linewidth=0.1)
        ax_E_cc_apex.plot(self.strain_transients['strain_t'],self.strain_transients['apex_E_cc_median'], color='k', linestyle='--')
        ax_E_cc_apex.plot(self.strain_transients['strain_t'], self.strain_transients['apex_E_cc_uq'], color='k', linestyle='-')
        ax_E_cc_apex.plot(self.strain_transients['strain_t'], self.strain_transients['apex_E_cc_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['base_E_rr'].shape[0]):
            ax_E_rr_base.plot(self.strain_transients['strain_t'], self.strain_transients['base_E_rr'][node_i, :], color='g', alpha=0.3, linewidth=0.1)
        ax_E_rr_base.plot(self.strain_transients['strain_t'],self.strain_transients['base_E_rr_median'], color='k', linestyle='--')
        ax_E_rr_base.plot(self.strain_transients['strain_t'], self.strain_transients['base_E_rr_uq'], color='k', linestyle='-')
        ax_E_rr_base.plot(self.strain_transients['strain_t'], self.strain_transients['base_E_rr_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['mid_E_rr'].shape[0]):
            ax_E_rr_mid.plot(self.strain_transients['strain_t'], self.strain_transients['mid_E_rr'][node_i, :], color='g', alpha=0.3, linewidth=0.1)
        ax_E_rr_mid.plot(self.strain_transients['strain_t'],self.strain_transients['mid_E_rr_median'], color='k', linestyle='--')
        ax_E_rr_mid.plot(self.strain_transients['strain_t'], self.strain_transients['mid_E_rr_uq'], color='k', linestyle='-')
        ax_E_rr_mid.plot(self.strain_transients['strain_t'], self.strain_transients['mid_E_rr_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['apex_E_rr'].shape[0]):
            ax_E_rr_apex.plot(self.strain_transients['strain_t'], self.strain_transients['apex_E_rr'][node_i, :], color='g', alpha=0.3, linewidth=0.1)
        ax_E_rr_apex.plot(self.strain_transients['strain_t'],self.strain_transients['apex_E_rr_median'], color='k', linestyle='--')
        ax_E_rr_apex.plot(self.strain_transients['strain_t'], self.strain_transients['apex_E_rr_uq'], color='k', linestyle='-')
        ax_E_rr_apex.plot(self.strain_transients['strain_t'], self.strain_transients['apex_E_rr_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['two_chamber_E_ll'].shape[0]):
            ax_E_ll_two.plot(self.strain_transients['strain_t'], self.strain_transients['two_chamber_E_ll'][node_i, :], color='b', alpha=0.3, linewidth=0.1)
        ax_E_ll_two.plot(self.strain_transients['strain_t'],self.strain_transients['two_chamber_E_ll_median'], color='k', linestyle='--')
        ax_E_ll_two.plot(self.strain_transients['strain_t'], self.strain_transients['two_chamber_E_ll_uq'], color='k', linestyle='-')
        ax_E_ll_two.plot(self.strain_transients['strain_t'], self.strain_transients['two_chamber_E_ll_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['four_chamber_E_ll'].shape[0]):
            ax_E_ll_four.plot(self.strain_transients['strain_t'], self.strain_transients['four_chamber_E_ll'][node_i, :], color='b', alpha=0.3, linewidth=0.1)
        ax_E_ll_four.plot(self.strain_transients['strain_t'],self.strain_transients['four_chamber_E_ll_median'], color='k', linestyle='--')
        ax_E_ll_four.plot(self.strain_transients['strain_t'], self.strain_transients['four_chamber_E_ll_uq'], color='k', linestyle='-')
        ax_E_ll_four.plot(self.strain_transients['strain_t'], self.strain_transients['four_chamber_E_ll_lq'], color='k', linestyle='-')

        for node_i in range(self.strain_transients['three_chamber_E_ll'].shape[0]):
            ax_E_ll_three.plot(self.strain_transients['strain_t'], self.strain_transients['three_chamber_E_ll'][node_i, :], color='b', alpha=0.3, linewidth=0.1)
        ax_E_ll_three.plot(self.strain_transients['strain_t'],self.strain_transients['three_chamber_E_ll_median'], color='k', linestyle='--')
        ax_E_ll_three.plot(self.strain_transients['strain_t'], self.strain_transients['three_chamber_E_ll_uq'], color='k', linestyle='-')
        ax_E_ll_three.plot(self.strain_transients['strain_t'], self.strain_transients['three_chamber_E_ll_lq'], color='k', linestyle='-')

        plt.show()




    def visualise_compare_with_clinical_ranges(self, beat):
        healthy_ranges = {'LVEDV': [88, 161], # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3 95% confidence interval
                          'LVESV': [31, 68], # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                          'LVSV': [49, 100], # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                          'LVEF': [51, 70], # % female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                          'LVESP': [14.3, 15.3], # [kPa] +- standard error of the mean https://link.springer.com/article/10.1007/s12265-018-9816-y/tables/1
                          'RVEDV': [85, 166], # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                          'RVESV': [27, 77], # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                          'RVSV': [48, 99], # [mL] female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                          'RVEF': [47,68], # % female UKB https://jcmr-online.biomedcentral.com/articles/10.1186/s12968-017-0327-9/tables/3
                          'RVESP': [2.8, 3.6], # [kPa] +- standard deviation https://www.sciencedirect.com/science/article/pii/016752739190221A
                          'AVPD': [1.4, 1.9], # [cm] ES - ED for controls https://journals.physiology.org/doi/full/10.1152/ajpheart.01148.2006 range
                          'apical_displacement': [-0.001, 0.51], # [cm] ES - ED for controls ttps://journals.physiology.org/doi/full/10.1152/ajpheart.01148.2006 range
                          'ED_wall_thickness': [0.571, 1.381], # [cm] min and max mean value from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                          'ES_wall_thickness': [1.07, 1.749], # [cm] min and max mean value from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                          'wall_thickening': [0.084, 0.429], # [cm] min and max mean values from 17 AHA regions https://www.sciencedirect.com/science/article/pii/S1361841519300489
                          'lambda_endo': [0.85, 0.88], # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                          'lambda_mid' : [0.85, 0.88], # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                          'lambda_epi': [0.85, 0.86], # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                          'cross_lambda_endo': [0.79, 0.81], # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                          'cross_lambda_mid' : [0.83, 0.86], #range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                          'cross_lambda_epi': [0.88, 0.9], # range https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
                          'QTc': [369.2, 399.2], # [ms] UKB normal LV mass interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                          'Tpe': [50, 72], # [ms] LTVA interquartile range UKB https://www.ahajournals.org/action/downloadSupplement?doi=10.1161%2FJAHA.121.025897&file=jah37786-sup-0001-data-tabs-figs.pdf
                          'T_duration': [100, 114], # [ms] UKB normal LV mass interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                          'QRS_duration': [81, 97], # [ms] UKB normal LV mass interquartile range https://academic.oup.com/ehjdh/article/4/4/316/7188159?login=true
                          'diastoic_duration': [0.45]}

        fig = plt.figure(tight_layout=True, figsize=(18, 10))
        gs = GridSpec(3, 6)
        ax_lv_pv = fig.add_subplot(gs[0:2, 0:2])
        ax_lv_pv.set_title('Left PV')
        ax_rv_pv = fig.add_subplot(gs[0:2, 2:4])
        ax_rv_pv.set_title('Right PV')
        ax_rv_pv.set_ylim([0, 5])
        ax_lambda = fig.add_subplot(gs[0, 4])
        ax_lambda.set_title('lambda')
        ax_E_ll = fig.add_subplot(gs[0, 5])
        ax_E_ll.set_title('E_ll')
        ax_E_rr = fig.add_subplot(gs[1, 4])
        ax_E_rr.set_title('E_rr')
        ax_E_cc = fig.add_subplot(gs[1, 5])
        ax_E_cc.set_title('E_cc')
        ax_V1 = fig.add_subplot(gs[2, 0])
        ax_V1.set_title('V1')
        ax_V2 = fig.add_subplot(gs[2, 1])
        ax_V1.set_title('V2')
        ax_V3 = fig.add_subplot(gs[2, 2])
        ax_V1.set_title('V3')
        ax_V4 = fig.add_subplot(gs[2, 3])
        ax_V1.set_title('V4')
        ax_V5 = fig.add_subplot(gs[2, 4])
        ax_V1.set_title('V5')
        ax_V6 = fig.add_subplot(gs[2, 5])
        ax_V1.set_title('V6')

        # Pressure volume loops
        ax_lv_pv.plot(self.pvs['vls'][beat-1], self.pvs['pls'][beat-1]/10000)
        ax_lv_pv.axvspan(healthy_ranges['LVEDV'][0], healthy_ranges['LVEDV'][1], alpha=0.3, color='green')
        ax_lv_pv.axvspan(healthy_ranges['LVESV'][0], healthy_ranges['LVESV'][1], alpha=0.3, color='green')
        ax_lv_pv.axhspan(healthy_ranges['LVESP'][0], healthy_ranges['LVESP'][1], alpha=0.3, color='green')

        ax_rv_pv.plot(self.pvs['vrs'][beat-1], self.pvs['prs'][beat-1]/10000)
        ax_rv_pv.axvspan(healthy_ranges['RVEDV'][0], healthy_ranges['RVEDV'][1], alpha=0.3, color='green')
        ax_rv_pv.axvspan(healthy_ranges['RVESV'][0], healthy_ranges['RVESV'][1], alpha=0.3, color='green')
        ax_rv_pv.axhspan(healthy_ranges['RVESP'][0], healthy_ranges['RVESP'][1], alpha=0.3, color='green')

        # ECGs
        ax_V1.plot((self.ecgs['ts'][beat-1]-self.simulation_dict['end_diastole_t'][0])*1000,
                   self.ecgs['V1s'][beat-1]/self.ecgs['max_all_leads'])
        ax_V1.axvspan(healthy_ranges['QTc'][0], healthy_ranges['QTc'][1], alpha=0.3, color='green')
        ax_V1.axvspan(healthy_ranges['QRS_duration'][0], healthy_ranges['QRS_duration'][1], alpha=0.3, color='green')
        ax_V2.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V2s'][beat - 1]/self.ecgs['max_all_leads'])
        ax_V3.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V3s'][beat - 1]/self.ecgs['max_all_leads'])
        ax_V4.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V4s'][beat - 1]/self.ecgs['max_all_leads'])
        ax_V5.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V5s'][beat - 1]/self.ecgs['max_all_leads'])
        ax_V6.plot((self.ecgs['ts'][beat - 1] - self.simulation_dict['end_diastole_t'][0]) * 1000,
                   self.ecgs['V6s'][beat - 1]/self.ecgs['max_all_leads'])
        ax_V1.set_xlim([0, 500])
        ax_V2.set_xlim([0, 500])
        ax_V3.set_xlim([0, 500])
        ax_V4.set_xlim([0, 500])
        ax_V5.set_xlim([0, 500])
        ax_V6.set_xlim([0, 500])
        ax_V1.set_ylim([-1, 1])
        ax_V2.set_ylim([-1, 1])
        ax_V3.set_ylim([-1, 1])
        ax_V4.set_ylim([-1, 1])
        ax_V5.set_ylim([-1, 1])
        ax_V6.set_ylim([-1, 1])


        # Fibre work
        ax_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['endo_lambda'], label='endo', color='C0')
        ax_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['midw_lambda'], label='mid', color='C1')
        ax_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['epi_lambda'], label='epi', color='C2')
        ax_lambda.legend()
        ax_lambda.axhspan(healthy_ranges['lambda_endo'][0], healthy_ranges['lambda_endo'][1], alpha=0.3, color='C0')
        ax_lambda.axhspan(healthy_ranges['lambda_mid'][0], healthy_ranges['lambda_mid'][1], alpha=0.3, color='C1')
        ax_lambda.axhspan(healthy_ranges['lambda_epi'][0], healthy_ranges['lambda_epi'][1], alpha=0.3, color='C2')

        # Displacements deformations
        ax_E_ll.plot(self.strain_transients['strain_t'], self.strain_transients['two_chamber_E_ll_median'], label='Two chamber', color='C0')
        ax_E_ll.plot(self.strain_transients['strain_t'], self.strain_transients['four_chamber_E_ll_median'], label='Four chamber', color='C1')
        ax_E_ll.plot(self.strain_transients['strain_t'], self.strain_transients['three_chamber_E_ll_median'], label='Three chamber',
                     color='C1')

        ax_E_ll.legend()
        ax_E_ll.axhspan(healthy_ranges['AVPD'][0], healthy_ranges['AVPD'][1], alpha=0.3, color='C0')
        ax_E_ll.axhspan(healthy_ranges['apical_displacement'][0], healthy_ranges['apical_displacement'][1], alpha=0.3, color='C1')

        ax_wall.plot(self.deformation_transients['deformation_t'], self.deformation_transients['lv_wall_thickness'])
        ax_wall.axhspan(healthy_ranges['ED_wall_thickness'][0], healthy_ranges['ED_wall_thickness'][1], alpha=0.3, color='green')
        ax_wall.axhspan(healthy_ranges['ES_wall_thickness'][0], healthy_ranges['ES_wall_thickness'][1], alpha=0.3, color='green')


        ax_cross_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['endo_cross_lambda'], label='endo', color='C0')
        ax_cross_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['midw_cross_lambda'], label='mid', color='C1')
        ax_cross_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['epi_cross_lambda'], label='epi', color='C2')
        ax_cross_lambda.legend()
        ax_lambda.axhspan(healthy_ranges['cross_lambda_endo'][0], healthy_ranges['cross_lambda_endo'][1], alpha=0.3, color='C0')
        ax_lambda.axhspan(healthy_ranges['cross_lambda_mid'][0], healthy_ranges['cross_lambda_mid'][1], alpha=0.3, color='C1')
        ax_lambda.axhspan(healthy_ranges['cross_lambda_epi'][0], healthy_ranges['cross_lambda_epi'][1], alpha=0.3, color='C2')

        plt.show()



def evaluate_lat(time, vm, percentage, time_window):
    print('evaluate Lat: vm shape: ', vm.shape)
    window_idx = np.nonzero((time > time_window[0]) & (time < time_window[1]))[0]
    vm = vm[:, window_idx]
    time = time[window_idx] - time_window[0]
    vm_range = np.amax(vm, axis=1) - np.amin(vm, axis=1)
    vm_threshold = vm_range * (1.0 - percentage / 100.0) + np.amin(vm, axis=1)
    activation_map = pymp.shared.array(vm.shape[0])
    time_shared = pymp.shared.array(time.shape)
    time_shared[:] = time
    threadsNum = int(multiprocessing.cpu_count()*0.7)
    with pymp.Parallel(min(threadsNum, vm.shape[0])) as p1:
        for node_i in p1.range(vm.shape[0]):
            local_vm = vm[node_i, :]
            if any(np.nonzero(local_vm > vm_threshold[node_i])[0]):
                index = np.nonzero(local_vm > vm_threshold[node_i])[0][0]  # np.searchsorted(local_vm, vm_threshold[node_i])
                activation_map[node_i] = time_shared[index]
            else:
                activation_map[node_i] = np.nan # Has not found activation within the time window.
    return activation_map


def evaluate_rt(time, vm, percentage, time_window):
    window_idx = np.nonzero((time > time_window[0]) & (time < time_window[1]))[0]
    vm = vm[:, window_idx]
    time = time[window_idx] - time_window[0]  # Offset by beginning of time window.
    repolarisation_map = np.ones(vm.shape[0])
    vm_range = np.amax(vm, axis=1) - np.amin(vm, axis=1)
    vm_threshold = vm_range * (1.0 - percentage / 100.0) + np.amin(vm, axis=1)
    vm_max_idx = np.argmax(vm, axis=1)
    for node_i in range(vm.shape[0]):  # Loop through every node in mesh
        local_vm = vm[node_i, vm_max_idx[node_i]:]
        local_vm_fliped = np.flip(local_vm)
        fliped_index = np.searchsorted(local_vm_fliped, vm_threshold[node_i])
        if fliped_index == 0:
            repolarisation_map[node_i] = np.nan # Has not found repolarisation within the time window.
        else:
            index = local_vm.shape[0] - fliped_index - 1
            repolarisation_map[node_i] = time[index + vm_max_idx[node_i]]
    return repolarisation_map


def mapIndices(points_to_map_xyz, reference_points_xyz,
               return_unique_only=False):  # TODO the unique should be done after this function
    mapped_indexes = pymp.shared.array((points_to_map_xyz.shape[0]), dtype=int)
    threadsNum = int(multiprocessing.cpu_count()*0.7)
    with pymp.Parallel(min(threadsNum, points_to_map_xyz.shape[0])) as p1:
        for conf_i in p1.range(points_to_map_xyz.shape[0]):
            mapped_indexes[conf_i] = np.argmin(
                np.linalg.norm(reference_points_xyz - points_to_map_xyz[conf_i, :], ord=2, axis=1)).astype(int)
    if return_unique_only:  # use the unique function without sorting the contents of the array (meta_indexes)
        unique_meta_indexes = np.unique(mapped_indexes, axis=0, return_index=True)[
            1]  # indexes to the indexes (meta_indexes) that are unique
        mapped_indexes = mapped_indexes[sorted(unique_meta_indexes)]  # TODO this could just be one line of code
    return mapped_indexes


def identify_alya_id_type(inputfolder, project_name):
    # read the header where element ids are stored and see if it's int8 or int4
    file_suffix = '.post.alyabin'
    filename = os.path.join(inputfolder, project_name + '-LNODS' + str(file_suffix))

    with open(filename, 'rb') as f:
        header = read_header_alyabin(f)
        if '32' in header['DataType']:
            alya_id_type = np.int32
        elif '64' in header['DataType']:
            alya_id_type = np.int64
        else:
            assert False, 'Alya id type ' + header[6] + ' is not supported'
    return alya_id_type


def read_alyabin_array(filename, number_of_blocks, alya_id_type):
    with open(filename, 'rb') as f:
        header = read_header_alyabin(f)
        number_of_dimensions = header['Columns']
        number_of_tuples_total = header['Lines']
        time_instant_int = header['TimeStepNo']
        time_instant_real = header['Time']
        # print(f'Reading array: {number_of_dimensions} dim, {number_of_tuples_total} tuples\n')

        datatype = np.dtype(header['DataType'])

        tuples = np.zeros((number_of_tuples_total, number_of_dimensions), dtype=datatype)

        c = 0

        tuples_per_block = np.zeros(number_of_blocks, dtype=np.int32)
        for i in range(number_of_blocks):
            number_of_tuples_in_block = read_one_fp90_record(f, 1, alya_id_type)[0]  # stored by alya
            tuples_per_block[i] = number_of_tuples_in_block

            # print(f'Block {i}/{number_of_blocks}: {(number_of_tuples_in_block)} tuples\n')
            tuples_temp = read_one_fp90_record(f, number_of_dimensions * number_of_tuples_in_block, datatype)

            tuples[c:c + number_of_tuples_in_block, :] = np.reshape(tuples_temp, (
                number_of_tuples_in_block, number_of_dimensions))
            c = c + number_of_tuples_in_block
    return {'tuples': tuples, 'header': header, 'tuples_per_block': tuples_per_block}


def read_one_fp90_record(file_object, number_of_elements, datatype):
    # fortran stores length with every block, in the beginning and the end
    count_read = 0
    record = []
    while count_read < number_of_elements:
        # in case the record is stored as several blocks join them
        block_len = np.fromfile(file_object, dtype=np.int32, count=1)
        # block_len is in bytes
        block = np.fromfile(file_object, dtype=datatype, count=block_len[0] // np.dtype(datatype).itemsize)
        block_len = np.fromfile(file_object, dtype=np.int32, count=1)
        count_read = count_read + block_len
        record = record + [block]

    return np.concatenate(record)


def read_header_alyabin(file_object):
    # a sanity check
    assert hasattr(file_object, 'read'), "read_header: argument is not a file object"

    ihead = read_one_fp90_record(file_object, 1, np.int32)  # ! Header: 1234, int32
    assert ihead[0] == 1234, "Header is not 1234"
    strings = []
    integers = []
    for i in range(0, 9):
        strings = strings + [read_one_fp90_record(file_object, 8, np.uint8).tobytes().decode().strip()]

    for i in range(0, 5):
        integers = integers + [read_one_fp90_record(file_object, 1, np.int32)]

    if (strings[1][0:5] != 'V0001'):
        integers = integers + [read_one_fp90_record(file_object, 1, np.int32)]
        # read(ii) integers(4) ! Time step, int32

    reals = read_one_fp90_record(file_object, 1, np.float64)
    # read(ii) reals(1)    ! Time, float64

    if (strings[1][0:5] == 'V0001'):
        integers[3] = int(reals)  # ! floor()?

    number_of_dimensions = integers[0][0]
    number_of_tuples_total = integers[1][0]
    time_instant_int = integers[3][0]
    time_instant_real = reals[0]

    if (strings[5] == 'REAL'):
        field_dtype = 'float'
    if (strings[5] == 'INTEG'):
        field_dtype = 'int'

    if (strings[4] == 'NPOIN'):
        association = 'node'
    else:
        association = 'element'

    if strings[6] == '4BYTE':
        field_dtype = field_dtype + '32'
    elif strings[6] == '8BYTE':
        field_dtype = field_dtype + '64'
    else:
        assert False, 'Alya id type ' + str(field_dtype) + ' is not supported'

    if (strings[3] == 'SCALA'):
        variabletype = 'scalar'
    elif (strings[3] == 'VECTO'):
        variabletype = 'vector'
    else:
        assert False, "unsupported type of variable"

    assert (strings[8] == 'NOFIL'), "Filtered types not supported"

    header = {'DataType': field_dtype, 'Lines': number_of_tuples_total, 'Columns': number_of_dimensions,
              'TimeStepNo': time_instant_int, \
              'Time': time_instant_real, 'Association': association, 'VariableType': variabletype}

    return header


