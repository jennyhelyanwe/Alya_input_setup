import matplotlib
matplotlib.use('TKagg')
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

class PostProcessing(MeshStructure):
    def __init__(self, alya, simulation_json_file, alya_output_dir, verbose):
        # super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.geometry = alya.geometry
        self.name = alya.name
        self.node_fields = alya.node_fields
        self.element_fields = alya.element_fields
        self.alya_output_dir = alya_output_dir
        self.alyacsv_dir = self.alya_output_dir+'/results_csv/'
        self.results_dir = self.alya_output_dir + '/results_analysis/'
        if not os.path.exists(self.results_dir):
            os.mkdir(self.results_dir)
        self.post_nodefield = Fields(self.name, field_type='postnodefield', verbose=verbose)
        print ('Reading in post_nodefield CSV from: ', self.results_dir)
        self.post_nodefield.read_csv_to_attributes(input_dir=self.results_dir, field_type='postnodefield')
        self.post_elementfield = Fields(self.name, field_type='postelementfield', verbose=verbose)
        print('Reading in post_elementfield CSV from: ', self.results_dir)
        self.post_elementfield.read_csv_to_attributes(input_dir=self.results_dir, field_type='postelementfield')
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.verbose = verbose
        self.qoi = {}
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

    def evaluate_ecg_biomarkers(self, beat):
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
        # lead_names = ['Is', 'IIs', 'IIIs', 'aVLs', 'aVRs', 'aVFs', 'V1s', 'V2s', 'V3s', 'V4s', 'V5s', 'V6s']
        # qrs_dur = np.zeros((len(lead_names)))
        # qt_dur = np.zeros((len(lead_names)))
        # t_pe = np.zeros((len(lead_names)))
        # t_peak = np.zeros((len(lead_names)))
        # qtpeak_dur = np.zeros((len(lead_names)))
        # t_polarity = np.zeros((len(lead_names)))
        # for i, lead_i in enumerate(lead_names):
        #     qrs_dur[i], qt_dur[i], t_pe[i], t_peak[i], qtpeak_dur[i], t_polarity[i], landmarks_temp = \
        #         self.calculate_ecg_biomarkers(time=self.ecgs['ts'][beat], V=self.ecgs[lead_i][beat],
        #                                       LAT=self.post_nodefield.dict['lat'])
        # qoi = {}
        # qoi['qrs_dur_mean'] = np.mean(qrs_dur)
        # qoi['qt_dur_mean'] = np.mean(qt_dur)
        # qoi['t_pe_mean'] = np.mean(t_pe)
        # qoi['t_peak_mean'] = np.mean(t_peak)
        # qoi['qtpeak_dur_mean'] = np.mean(qtpeak_dur)
        # qoi['t_polarity_mean'] = np.mean(t_polarity)
        # qoi['qrs_dur_std'] = np.std(qrs_dur)
        # qoi['qt_dur_std'] = np.std(qt_dur)
        # qoi['t_pe_std'] = np.std(t_pe)
        # qoi['t_peak_std'] = np.std(t_peak)
        # qoi['qtpeak_dur_std']  = np.std(qtpeak_dur)
        # qoi['t_polarity_std'] = np.std(t_polarity)
        # self.qoi.update(qoi)

    def evaluate_pv_biomarkers(self, beat):
        pv_filename = self.alya_output_dir + self.simulation_dict['name']
        ecgpv = ECGPV_visualisation(CL=1.0)
        print('Reading PV from: ', pv_filename)
        self.pvs = ecgpv._read_PV(pv_filename)
        pv_analysis = ecgpv.analysis_PV(pvs=self.pvs, beat=beat)
        qoi = {}
        qoi['EDVL'] = pv_analysis['EDVL']
        qoi['EDVR'] = pv_analysis['EDVR']
        qoi['ESVL'] = pv_analysis['ESVL']
        qoi['ESVR'] = pv_analysis['ESVR']
        qoi['LVEF'] = pv_analysis['LVEF']
        qoi['RVEF'] = pv_analysis['RVEF']
        qoi['PmaxL'] = pv_analysis['PmaxL']
        qoi['PmaxR'] = pv_analysis['PmaxR']
        qoi['SVL'] = pv_analysis['SVL']
        qoi['SVR'] = pv_analysis['SVR']
        self.qoi.update(qoi)

    def evaluate_deformation_biomarkers(self, beat):
        self.deformation_transients = {}
        self.read_csv_fields(read_field_name='DISPL', read_field_type='vector')
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
        print('Evaluating AVPD')
        base_cutoff = 0.7
        truncated_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] < base_cutoff)[0]
        basal_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] >= base_cutoff)[0]
        mean_ab_vector = np.mean(self.node_fields.dict['longitudinal-vector'][truncated_mesh_nodes, :], axis=0)
        apical_cutoff = 0.2
        apical_mesh_nodes = np.nonzero(self.node_fields.dict['ab'] <= apical_cutoff)[0]
        avpd = pymp.shared.array(time_idx.shape[0], dtype=float)
        apical_displacement = pymp.shared.array(time_idx.shape[0], dtype=float)
        displacement_shared = pymp.shared.array(self.post_nodefield.dict['DISPL'].shape, dtype=float)
        displacement_shared[:,:,:] = self.post_nodefield.dict['DISPL']
        threadsNum = multiprocessing.cpu_count()
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
        qoi['peak_apical_displacement'] = np.amax(apical_displacement)
        qoi['min_apical_displacement'] = np.amin(apical_displacement)
        qoi['peak_lv_wall_thickness'] = np.amax(lv_wall_thickness)
        qoi['min_lv_wall_thickness'] = np.amin(lv_wall_thickness)
        qoi['peak_rv_wall_thickness'] = np.amax(rv_wall_thickness)
        qoi['min_rv_wall_thickness'] = np.amin(rv_wall_thickness)
        self.qoi.update(qoi)

    def evaluate_fibre_work_biomarkers(self, beat):
        self.fibre_work = {}
        self.read_csv_fields(read_field_name='LAMBD', read_field_type='vector')
        self.read_csv_fields(read_field_name='ACTST', read_field_type='vector')
        # Select time segment according to specified beat # TODO assuming CL is always 1.0 s
        CL = 1.0
        time_idx = []
        deformation_t = []
        for time_i in range(self.post_nodefield.dict['time'].shape[0]):
            if (self.post_nodefield.dict['time'][time_i] > (beat - 1) * CL) & (
                    self.post_nodefield.dict['time'][time_i] < beat * CL):
                deformation_t.append(self.post_nodefield.dict['time'][time_i] - (beat - 1) * CL)
                time_idx.append(time_i)
        fibrework_t = np.array(deformation_t)
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

        ta_shared = pymp.shared.array(self.post_nodefield.dict['ACTST'].shape, dtype=float)
        lambda_shared = pymp.shared.array(self.post_nodefield.dict['LAMBD'].shape, dtype=float)
        ta_shared[:,:] = self.post_nodefield.dict['ACTST']
        lambda_shared[:,:] = self.post_nodefield.dict['LAMBD']
        threadsNum = multiprocessing.cpu_count()
        with pymp.Parallel(min(threadsNum, time_idx.shape[0])) as p1:
            for time_i in p1.range(time_idx.shape[0]):
                apical_lambda[time_i] = np.mean(lambda_shared[apical_mesh_nodes, 0, time_i])
                midv_lambda[time_i] = np.mean(lambda_shared[midv_mesh_nodes, 0, time_i])
                basal_lambda[time_i] = np.mean(lambda_shared[basal_mesh_nodes, 0, time_i])
                endo_lambda[time_i] = np.mean(lambda_shared[endo_mesh_nodes, 0, time_i])
                midw_lambda[time_i] = np.mean(lambda_shared[midw_mesh_nodes, 0, time_i])
                epi_lambda[time_i] = np.mean(lambda_shared[epi_mesh_nodes, 0, time_i])
                endo_cross_lambda[time_i] = np.mean(lambda_shared[endo_mesh_nodes, 1, time_i])
                midw_cross_lambda[time_i] = np.mean(lambda_shared[midw_mesh_nodes, 1, time_i])
                epi_cross_lambda[time_i] = np.mean(lambda_shared[epi_mesh_nodes, 1, time_i])

                apical_Ta[time_i] = np.mean(ta_shared[apical_mesh_nodes, 0, time_i])
                midv_Ta[time_i] = np.mean(ta_shared[midv_mesh_nodes, 0, time_i])
                basal_Ta[time_i] = np.mean(ta_shared[basal_mesh_nodes, 0, time_i])
                endo_Ta[time_i] = np.mean(ta_shared[endo_mesh_nodes, 0, time_i])
                midw_Ta[time_i] = np.mean(ta_shared[midw_mesh_nodes, 0, time_i])
                epi_Ta[time_i] = np.mean(ta_shared[epi_mesh_nodes, 0, time_i])
        self.fibre_work = {}
        self.fibre_work['fibrework_t'] = fibrework_t
        self.fibre_work['apical_lambda'] = apical_lambda
        self.fibre_work['midv_lambda'] = midv_lambda
        self.fibre_work['basal_lambda'] = basal_lambda
        self.fibre_work['endo_lambda'] = endo_lambda
        self.fibre_work['midw_lambda'] = midw_lambda
        self.fibre_work['epi_lambda'] = epi_lambda
        self.fibre_work['endo_cross_lambda'] = endo_cross_lambda
        self.fibre_work['midw_cross_lambda'] = midw_cross_lambda
        self.fibre_work['epi_cross_lambda'] = epi_cross_lambda

        self.fibre_work['apical_Ta'] = apical_Ta
        self.fibre_work['midv_Ta'] = midv_Ta
        self.fibre_work['basal_Ta'] = basal_Ta
        self.fibre_work['endo_Ta'] = endo_Ta
        self.fibre_work['midw_Ta'] = midw_Ta
        self.fibre_work['epi_Ta'] = epi_Ta
        # Save QoIs
        qoi = {}
        qoi['peak_Ta'] = np.amax(np.mean(ta_shared[:,0, :], axis=0))
        qoi['min_Ta'] = np.amin(np.mean(ta_shared[:,0, :], axis=0))
        qoi['peak_lambda'] = np.amax(np.mean(lambda_shared[:,0, :], axis=0))
        qoi['min_lambda'] = np.amin(np.mean(lambda_shared[:, 0, :], axis=0))
        self.qoi.update(qoi)


    def calculate_ecg_biomarkers(self, time, V, LAT):
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
        q_start_idx = np.nanargmin(LAT)

        # Set QRS end
        qrs_end_idx = np.nanargmax(LAT).astype(int)
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

    def read_csv_fields(self, read_field_name, read_field_type):
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
        if read_field_name:
            temp = np.zeros((self.geometry.number_of_nodes, time.shape[0]))
            field_name = read_field_name
            field_type = read_field_type
            print('Reading csv for field: '  + 'from : '+self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-******.csv')
            if field_type == 'scalar':
                temp = pymp.shared.array((self.geometry.number_of_nodes, time.shape[0]), dtype=float)
                threadsNum = 50 # multiprocessing.cpu_count()
                with pymp.Parallel(min(threadsNum, time_index_shared.shape[0])) as p1:
                    for time_i in p1.range(time_index_shared.shape[0]):
                # temp = np.zeros((self.geometry.number_of_nodes, time.shape[0]))
                # print(temp.shape)
                # if True:
                #     for time_i in range(time_index.shape[0]):
                        index = '{:06d}'.format(time_index_shared[time_i])
                        filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
                        temp[:, time_i] = pd.read_csv(filename, delimiter=',', header=None).values.transpose()[0] #.astype(float)
                self.post_nodefield.dict[field_name] = temp
            elif field_type == 'vector':
                temp = pymp.shared.array((self.geometry.number_of_nodes, 3, time.shape[0]), dtype=float)
                threadsNum = multiprocessing.cpu_count()
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
        #         temp = np.zeros((self.geometry.number_of_nodes, time.shape[0]))
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
        self.read_csv_fields(read_field_name='EPSXX', read_field_type='scalar')
        self.read_csv_fields(read_field_name='EPSYY', read_field_type='scalar')
        self.read_csv_fields(read_field_name='EPSZZ', read_field_type='scalar')
        self.read_csv_fields(read_field_name='EPSXY', read_field_type='scalar')
        self.read_csv_fields(read_field_name='EPSXZ', read_field_type='scalar')
        self.read_csv_fields(read_field_name='EPSYZ', read_field_type='scalar')

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
        threadsNum = multiprocessing.cpu_count()
        for time_i in range(time_idx.shape[0]):
            E[:, time_i, 0, 0] = exx[:, time_i]
            E[:, time_i, 0, 1] = exy[:, time_i]
            E[:, time_i, 0, 2] = exz[:, time_i]
            E[:, time_i, 1, 0] = exy[:, time_i]
            E[:, time_i, 1, 1] = eyy[:, time_i]
            E[:, time_i, 1, 2] = eyz[:, time_i]
            E[:, time_i, 2, 0] = exz[:, time_i]
            E[:, time_i, 2, 1] = eyz[:, time_i]
            E[:, time_i, 2, 2] = ezz[:, time_i]
            with pymp.Parallel(min(threadsNum, self.geometry.number_of_nodes)) as p1:
                for node_i in p1.range(self.geometry.number_of_nodes):
                    E_cc[node_i, time_i] = np.dot(np.dot(E[node_i, time_i, :, :], local_c[node_i, :]), local_c[node_i, :])
                    E_ll[node_i, time_i] = np.dot(np.dot(E[node_i, time_i, :, :], local_l[node_i, :]), local_l[node_i, :])
                    E_rr[node_i, time_i] = np.dot(np.dot(E[node_i, time_i, :, :], local_r[node_i, :]), local_r[node_i, :])

        # Save strains in short and long axis slices for comparison with https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28724
        base_short_axis_nodes = np.nonzero(self.node_fields.dict['short-axis-slices'] == 3)[0]
        mid_short_axis_nodes = np.nonzero(self.node_fields.dict['short-axis-slices'] == 2)[0]
        apex_short_axis_nodes = np.nonzero(self.node_fields.dict['short-axis-slices'] == 1)[0]
        two_chamber_long_axis_nodes = np.nonzero(self.node_fields.dict['long-axis-slices'] == 1)[0]
        four_chamber_long_axis_nodes = np.nonzero(self.node_fields.dict['long-axis-slices'] == 2)[0]
        three_chamber_long_axis_nodes = np.nonzero(self.node_fields.dict['long-axis-slices'] == 3)[0]
        base_E_cc = pymp.shared.array((time_idx.shape[0], base_short_axis_nodes.shape[0]))
        mid_E_cc = pymp.shared.array((time_idx.shape[0], mid_short_axis_nodes.shape[0]))
        apex_E_cc = pymp.shared.array((time_idx.shape[0], apex_short_axis_nodes.shape[0]))
        base_E_rr = pymp.shared.array((time_idx.shape[0], base_short_axis_nodes.shape[0]))
        mid_E_rr = pymp.shared.array((time_idx.shape[0], mid_short_axis_nodes.shape[0]))
        apex_E_rr = pymp.shared.array((time_idx.shape[0], apex_short_axis_nodes.shape[0]))
        two_chamber_E_ll = pymp.shared.array((time_idx.shape[0], two_chamber_long_axis_nodes.shape[0]))
        four_chamber_E_ll = pymp.shared.array((time_idx.shape[0], four_chamber_long_axis_nodes.shape[0]))
        three_chamber_E_ll = pymp.shared.array((time_idx.shape[0], three_chamber_long_axis_nodes.shape[0]))
        with pymp.Parallel(min(threadsNum, time_idx.shape[0])) as p1:
            for time_i in p1.range(time_idx.shape[0]):
                # Circumferential strains
                base_E_cc[time_i, :] = E_cc[time_i, base_short_axis_nodes]
                mid_E_cc[time_i, :] = E_cc[time_i, mid_short_axis_nodes]
                apex_E_cc[time_i, :] = E_cc[time_i, apex_short_axis_nodes]
                # Radial strains
                base_E_rr[time_i, :] = E_rr[time_i, base_short_axis_nodes]
                mid_E_rr[time_i, :] = E_rr[time_i, mid_short_axis_nodes]
                apex_E_rr[time_i, :] = E_rr[time_i, apex_short_axis_nodes]
                # Longitudinal strains
                two_chamber_E_ll[time_i, :] = E_ll[time_i, two_chamber_long_axis_nodes]
                four_chamber_E_ll[time_i, :] = E_ll[time_i, four_chamber_long_axis_nodes]
                three_chamber_E_ll[time_i, :] = E_ll[time_i, three_chamber_long_axis_nodes]

        # Evaluate median and upper and lower quartiles of traces
        # Circumferential
        systole = np.min(base_E_cc, axis=0)
        base_E_cc_median = base_E_cc[:, np.where(systole == np.percentile(systole, 50))[0]]
        base_E_cc_uq = base_E_cc[:, np.where(systole == np.percentile(systole, 75))[0]]
        base_E_cc_lq = base_E_cc[:, np.where(systole == np.percentile(systole, 25))[0]]
        systole = np.min(mid_E_cc, axis=0)
        mid_E_cc_median = mid_E_cc[:, np.where(systole == np.percentile(systole, 50))[0]]
        mid_E_cc_uq = mid_E_cc[:, np.where(systole == np.percentile(systole, 75))[0]]
        mid_E_cc_lq = mid_E_cc[:, np.where(systole == np.percentile(systole, 25))[0]]
        systole = np.min(apex_E_cc, axis=0)
        apex_E_cc_median = apex_E_cc[:, np.where(systole == np.percentile(systole, 50))[0]]
        apex_E_cc_uq = apex_E_cc[:, np.where(systole == np.percentile(systole, 75))[0]]
        apex_E_cc_lq = apex_E_cc[:, np.where(systole == np.percentile(systole, 25))[0]]
        # Radial
        systole = np.min(base_E_rr, axis=0)
        base_E_rr_median = base_E_rr[:, np.where(systole == np.percentile(systole, 50))[0]]
        base_E_rr_uq = base_E_rr[:, np.where(systole == np.percentile(systole, 75))[0]]
        base_E_rr_lq = base_E_rr[:, np.where(systole == np.percentile(systole, 25))[0]]
        systole = np.min(mid_E_rr, axis=0)
        mid_E_rr_median = mid_E_rr[:, np.where(systole == np.percentile(systole, 50))[0]]
        mid_E_rr_uq = mid_E_rr[:, np.where(systole == np.percentile(systole, 75))[0]]
        mid_E_rr_lq = mid_E_rr[:, np.where(systole == np.percentile(systole, 25))[0]]
        systole = np.min(apex_E_rr, axis=0)
        apex_E_rr_median = apex_E_rr[:, np.where(systole == np.percentile(systole, 50))[0]]
        apex_E_rr_uq = apex_E_rr[:, np.where(systole == np.percentile(systole, 75))[0]]
        apex_E_rr_lq = apex_E_rr[:, np.where(systole == np.percentile(systole, 25))[0]]
        # Longitudinal
        systole = np.min(two_chamber_E_ll, axis=0)
        two_chamber_E_ll_median = two_chamber_E_ll[:, np.where(systole == np.percentile(systole, 50))[0]]
        two_chamber_E_ll_uq = two_chamber_E_ll[:, np.where(systole == np.percentile(systole, 75))[0]]
        two_chamber_E_ll_lq = two_chamber_E_ll[:, np.where(systole == np.percentile(systole, 25))[0]]
        systole = np.min(four_chamber_E_ll, axis=0)
        four_chamber_E_ll_median = four_chamber_E_ll[:, np.where(systole == np.percentile(systole, 50))[0]]
        four_chamber_E_ll_uq = four_chamber_E_ll[:, np.where(systole == np.percentile(systole, 75))[0]]
        four_chamber_E_ll_lq = four_chamber_E_ll[:, np.where(systole == np.percentile(systole, 25))[0]]
        systole = np.min(three_chamber_E_ll, axis=0)
        three_chamber_E_ll_median = three_chamber_E_ll[:, np.where(systole == np.percentile(systole, 50))[0]]
        three_chamber_E_ll_uq = three_chamber_E_ll[:, np.where(systole == np.percentile(systole, 75))[0]]
        three_chamber_E_ll_lq = three_chamber_E_ll[:, np.where(systole == np.percentile(systole, 25))[0]]

        self.strains = {}
        self.strains['base_E_cc'] = base_E_cc
        self.strains['base_E_cc_median'] = base_E_cc_median
        self.strains['base_E_cc_uq'] = base_E_cc_uq
        self.strains['base_E_cc_lq'] = base_E_cc_lq

        self.strains['mid_E_cc'] = mid_E_cc
        self.strains['mid_E_cc_median'] = mid_E_cc_median
        self.strains['mid_E_cc_uq'] = mid_E_cc_uq
        self.strains['mid_E_cc_lq'] = mid_E_cc_lq

        self.strains['apex_E_cc'] = apex_E_cc
        self.strains['apex_E_cc_median'] = apex_E_cc_median
        self.strains['apex_E_cc_uq'] = apex_E_cc_uq
        self.strains['apex_E_cc_lq'] = apex_E_cc_lq

        self.strains['base_E_rr'] = base_E_rr
        self.strains['base_E_rr_median'] = base_E_rr_median
        self.strains['base_E_rr_uq'] = base_E_rr_uq
        self.strains['base_E_rr_lq'] = base_E_rr_lq

        self.strains['mid_E_rr'] = mid_E_rr
        self.strains['mid_E_rr_median'] = mid_E_rr_median
        self.strains['mid_E_rr_uq'] = mid_E_rr_uq
        self.strains['mid_E_rr_lq'] = mid_E_rr_lq

        self.strains['apex_E_rr'] = apex_E_rr
        self.strains['apex_E_rr_median'] = apex_E_rr_median
        self.strains['apex_E_rr_uq'] = apex_E_rr_uq
        self.strains['apex_E_rr_lq'] = apex_E_rr_lq

        self.strains['two_chamber_E_ll'] = two_chamber_E_ll
        self.strains['two_chamber_E_ll_median'] = two_chamber_E_ll_median
        self.strains['two_chamber_E_ll_uq'] = two_chamber_E_ll_uq
        self.strains['two_chamber_E_ll_lq'] = two_chamber_E_ll_lq

        self.strains['four_chamber_E_ll'] = four_chamber_E_ll
        self.strains['four_chamber_E_ll_median'] = four_chamber_E_ll_median
        self.strains['four_chamber_E_ll_uq'] = four_chamber_E_ll_uq
        self.strains['four_chamber_E_ll_lq'] = four_chamber_E_ll_lq

        self.strains['three_chamber_E_ll'] = three_chamber_E_ll
        self.strains['three_chamber_E_ll_median'] = three_chamber_E_ll_median
        self.strains['three_chamber_E_ll_uq'] = three_chamber_E_ll_uq
        self.strains['three_chamber_E_ll_lq'] = three_chamber_E_ll_lq

        qoi = {}
        qoi['max_median_mid_Ecc'] = np.max(mid_E_cc_median)
        qoi['min_median_mid_Ecc'] = np.min(mid_E_cc_median)
        qoi['max_median_mid_Err'] = np.max(mid_E_rr_median)
        qoi['min_median_mid_Err'] = np.min(mid_E_rr_median)
        qoi['max_median_four_chamber_Ell'] = np.max(four_chamber_E_ll_median)
        qoi['min_median_four_chamber_Ell'] = np.min(four_chamber_E_ll_median)
        self.qoi.update(qoi)



    def visualise_compare_pv_with_clinical_ranges(self, beat):
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
        ax_avpd = fig.add_subplot(gs[0, 4])
        ax_avpd.set_title('AVPD')
        ax_lambda = fig.add_subplot(gs[0, 5])
        ax_lambda.set_title('Lambda')
        ax_wall = fig.add_subplot(gs[1, 4])
        ax_wall.set_title('Wall thickness')
        ax_cross_lambda = fig.add_subplot(gs[1, 5])
        ax_cross_lambda.set_title('Cross lambda')
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

        # Displacements deformations
        ax_avpd.plot(self.deformation_transients['deformation_t'], self.deformation_transients['avpd'], label='AVPD', color='C0')
        ax_avpd.plot(self.deformation_transients['deformation_t'], self.deformation_transients['apical_displacement'], label='Apex displacemet', color='C1')
        ax_avpd.legend()
        ax_avpd.axhspan(healthy_ranges['AVPD'][0], healthy_ranges['AVPD'][1], alpha=0.3, color='C0')
        ax_avpd.axhspan(healthy_ranges['apical_displacement'][0], healthy_ranges['apical_displacement'][1], alpha=0.3, color='C1')

        ax_wall.plot(self.deformation_transients['deformation_t'], self.deformation_transients['lv_wall_thickness'])
        ax_wall.axhspan(healthy_ranges['ED_wall_thickness'][0], healthy_ranges['ED_wall_thickness'][1], alpha=0.3, color='green')
        ax_wall.axhspan(healthy_ranges['ES_wall_thickness'][0], healthy_ranges['ES_wall_thickness'][1], alpha=0.3, color='green')

        # Fibre work
        ax_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['endo_lambda'], label='endo',color='C0')
        ax_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['midw_lambda'], label='mid',color='C1')
        ax_lambda.plot(self.fibre_work['fibrework_t'], self.fibre_work['epi_lambda'], label='epi',color='C2')
        ax_lambda.legend()
        ax_lambda.axhspan(healthy_ranges['lambda_endo'][0], healthy_ranges['lambda_endo'][1], alpha=0.3, color='C0')
        ax_lambda.axhspan(healthy_ranges['lambda_mid'][0], healthy_ranges['lambda_mid'][1], alpha=0.3, color='C1')
        ax_lambda.axhspan(healthy_ranges['lambda_epi'][0], healthy_ranges['lambda_epi'][1], alpha=0.3, color='C2')

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
    threadsNum = multiprocessing.cpu_count()
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
    threadsNum = multiprocessing.cpu_count()
    with pymp.Parallel(min(threadsNum, points_to_map_xyz.shape[0])) as p1:
        for conf_i in p1.range(points_to_map_xyz.shape[0]):
            mapped_indexes[conf_i] = np.argmin(
                np.linalg.norm(reference_points_xyz - points_to_map_xyz[conf_i, :], ord=2, axis=1)).astype(int)
    if return_unique_only:  # use the unique function without sorting the contents of the array (meta_indexes)
        unique_meta_indexes = np.unique(mapped_indexes, axis=0, return_index=True)[
            1]  # indexes to the indexes (meta_indexes) that are unique
        mapped_indexes = mapped_indexes[sorted(unique_meta_indexes)]  # TODO this could just be one line of code
    return mapped_indexes



