import numpy as np
import os
import json
# from mpi4py import MPI #needs install
# from multiprocessing import Queue
import pymp, multiprocessing
from matplotlib import pyplot as plt

from meshstructure import MeshStructure
from myformat import Fields


class PostProcessing(MeshStructure):
    def __init__(self, name, geometric_data_dir, simulation_json_file, alyacsv_dir, alya_output_dir, results_dir, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.alyacsv_dir = alyacsv_dir
        self.alya_output_dir = alya_output_dir
        self.results_dir = results_dir
        self.simulation_biomarkers = {}
        self.simulation_data = {}
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.post_nodefield = Fields(name, field_type='postnodefield', verbose=verbose)
        self.post_elementfield = Fields(name, field_type='postelementfield', verbose=verbose)
        print('Read Alya log outputs')
        self.read_log_files()
        print('Reading postprocessing csv fields')
        self.post_nodefield.read_csv_to_attributes(input_dir=results_dir, field_type='postnodefield')
        # self.post_elementfield.read_csv_to_attributes(input_dir=results_dir, field_type='postelementfield')
        # Evaluate various biomarkers
        self.read_csv_fields()
        # self.evaluate_ep_maps()
        # self.evaluate_ventricular_cvs()
        if 'SOLIDZ' in self.simulation_dict['physics']:
            self.evaluate_basal_displacement()

        self.save_postprocessing()
        print('Writing out simulation biomarkers to '+self.results_dir+'simulation_biomarkers.json')
        print(self.simulation_biomarkers)
        json.dump(self.simulation_biomarkers, open(self.results_dir+'simulation_biomarkers.json', 'w'))

    def save_postprocessing(self):
        # Save to CSV
        print('Saving geometry and fields to CSV at '+self.results_dir)
        self.geometry.save_to_csv(self.results_dir)
        self.post_nodefield.save_to_csv(self.results_dir)
        self.post_elementfield.save_to_csv(self.results_dir)
        # Save to ensight
        print('Saving geometry and fields to Ensight at '+self.results_dir+'ensight/')
        self.geometry.save_to_ensight(self.results_dir + 'ensight/')
        self.post_nodefield.save_to_ensight(output_dir=self.results_dir + 'ensight/',
                                            casename=self.name + '_postnodefield',
                                            geometry=self.geometry)
        self.post_elementfield.save_to_ensight(output_dir=self.results_dir + 'ensight/',
                                               casename=self.name + '_postelementfield',
                                               geometry=self.geometry)


    def read_log_files(self):
        rawdata = np.loadtxt(self.alya_output_dir+self.simulation_dict['name']+'.exm.vin', skiprows=8)
        self.ecg_time = rawdata[:, 1]
        LA = rawdata[:, 3]
        RA = rawdata[:, 4]
        LL = rawdata[:, 5]
        RL = rawdata[:, 6]
        V1 = rawdata[:, 7]
        V2 = rawdata[:, 8]
        V3 = rawdata[:, 9]
        V4 = rawdata[:, 10]
        V5 = rawdata[:, 11]
        V6 = rawdata[:, 12]
        VW = 1/3 * (RA + LA + LL)
        V1 = V1 - VW
        V2 = V2 - VW
        V3 = V3 - VW
        V4 = V4 - VW
        V5 = V5 - VW
        V6 = V6 - VW
        I = LA - RA
        II = LL - RA
        III = LL - LA
        aVL = LA - (RA + LL)/2
        aVF = LL - (LA + RA)/2
        aVR = RA - (LA + LL)/2
        max_precordial_leads = np.amax(abs(np.concatenate((V1, V2, V3, V4, V5, V6))))
        max_limb_leads = np.amax(abs(np.concatenate((I, II, III, aVR, aVL, aVF))))
        self.I = I/max_limb_leads
        self.II = II/max_limb_leads
        self.III = III/max_limb_leads
        self.aVL = aVL/max_limb_leads
        self.aVF = aVF/max_limb_leads
        self.aVR = aVR/max_limb_leads
        self.V1 = V1/max_precordial_leads
        self.V2 = V2/max_precordial_leads
        self.V3 = V3/max_precordial_leads
        self.V4 = V4/max_precordial_leads
        self.V5 = V5/max_precordial_leads
        self.V6 = V6/max_precordial_leads

        if 'SOLIDZ' in self.simulation_dict['physics']:
            rawdata = np.loadtxt(self.alya_output_dir + self.simulation_dict['name'] + '-cardiac-cycle.sld.res', skiprows=18)
            self.pv_time = rawdata[:, 1]
            self.lv_cycle = rawdata[:, 3]
            self.lv_phase = rawdata[:, 4]
            self.lv_volume = rawdata[:, 5]
            self.lv_pressure = rawdata[:, 6]
            self.lv_arterial_pressure = rawdata[:, 7]
            self.rv_cycle= rawdata[:, 9]
            self.rv_phase = rawdata[:, 10]
            self.rv_volume = rawdata[:, 11]
            self.rv_pressure = rawdata[:, 12]
            self.rv_arterial_pressure = rawdata[:, 13]

    def read_csv_fields(self):
        name = self.simulation_dict['name']
        self.simulation_data['time'] =np.loadtxt(self.alyacsv_dir + 'timeset_1.csv', delimiter=',').astype(float)
        time = self.simulation_data['time']
        time_index = np.loadtxt(self.alyacsv_dir + 'timeindices_1.csv', delimiter=',').astype(int)
        assert time.shape == time_index.shape
        filenames = np.array([f for f in os.listdir(self.alyacsv_dir) if
                              os.path.isfile(os.path.join(self.alyacsv_dir, f)) and ('scalar' in f or 'vector' in f)])
        field_names_types = []
        for filename in filenames:
            field_names_types.append([filename.split('.')[2].split('-')[0], filename.split('.')[1]])
        field_names_types = np.unique(field_names_types, axis=0)
        for field_i in range(field_names_types.shape[0]):
            field_name = field_names_types[field_i, 0]
            field_type = field_names_types[field_i, 1]
            print('Reading csv for field: '+field_name + ', of type: '+field_type)
            if field_type == 'scalar':
                temp = pymp.shared.array((self.geometry.number_of_nodes, time.shape[0]), dtype=float)
                threadsNum = multiprocessing.cpu_count()
                with pymp.Parallel(min(threadsNum, time_index.shape[0])) as p1:
                    for time_i in p1.range(time_index.shape[0]):
                        index = '{:06d}'.format(time_index[time_i])
                        filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
                        temp[:, time_i] = np.loadtxt(filename, delimiter=',').astype(float)
                self.simulation_data[field_name] = temp
            elif field_type == 'vector':
                temp = pymp.shared.array((self.geometry.number_of_nodes, 3, time.shape[0]), dtype=float)
                threadsNum = multiprocessing.cpu_count()
                with pymp.Parallel(min(threadsNum, time_index.shape[0])) as p1:
                    for time_i in p1.range(time_index.shape[0]):
                        index = '{:06d}'.format(time_index[time_i])
                        filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
                        temp[:, :, time_i] = np.loadtxt(filename, delimiter=',').astype(float)
                self.simulation_data[field_name] = temp
            # self.post_nodefield.add_field(data=np.array(temp), data_name=field_name, field_type='postnodefield')
    def evaluate_ep_maps(self):
        print('Evaluating LAT')
        lat = evaluate_lat(time=self.simulation_data['time'], vm=self.simulation_data['INTRA'], percentage=20,
                           time_window=[float(self.simulation_dict['exmedi_delay_time']),
                                        float(self.simulation_dict['cycle_length'])])
        print('Evaluating RT')
        rt = evaluate_rt(time=self.simulation_data['time'], vm=self.simulation_data['INTRA'], percentage=90,
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
        self.simulation_biomarkers['mean_lv_transmural_cv'] = np.ma.masked_invalid(lv_transmural_cvs).mean()
        self.simulation_biomarkers['mean_rv_transmural_cv'] = np.ma.masked_invalid(rv_transmural_cvs).mean()

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

        mean_basal_ab_displacement_transient = np.zeros(self.simulation_data['time'].shape[0])
        mean_apical_ab_displacement_transient = np.zeros(self.simulation_data['time'].shape[0])
        mean_apicobasal_sum_displacement_transient = np.zeros(self.simulation_data['time'].shape[0])
        for time_i in range(self.simulation_data['time'].shape[0]):
            mean_basal_ab_displacement_transient[time_i] = np.mean(np.dot(self.simulation_data['DISPL'][
                                                                          basal_mesh_nodes, :, time_i],
                                                                          mean_ab_vector))
            mean_apical_ab_displacement_transient[time_i] = np.mean(np.dot(self.simulation_data['DISPL'][
                                                                           apical_mesh_nodes, :, time_i],
                                                                           mean_ab_vector))
            mean_apicobasal_sum_displacement_transient[time_i] = mean_basal_ab_displacement_transient[time_i] + mean_apical_ab_displacement_transient[time_i]
        self.simulation_biomarkers['max_basal_ab_displacement'] = np.amax(mean_basal_ab_displacement_transient)
        self.simulation_biomarkers['min_basal_ab_displacement'] = np.amin(mean_basal_ab_displacement_transient)
        self.simulation_biomarkers['max_apical_ab_displacement'] = np.amax(mean_apical_ab_displacement_transient)
        self.simulation_biomarkers['min_apical_ab_displacement'] = np.amin(mean_apical_ab_displacement_transient)
        self.simulation_biomarkers['max_longitudinal_strain'] = np.amax(mean_apicobasal_sum_displacement_transient)/reference_apex_to_base_length
        self.simulation_biomarkers['mean_longitudinal_strain'] = np.mean(mean_apicobasal_sum_displacement_transient)/reference_apex_to_base_length
        self.simulation_biomarkers['min_longitudinal_strain'] = np.amin(mean_apicobasal_sum_displacement_transient)/reference_apex_to_base_length
        np.savetxt(self.results_dir + 'mean_basal_ab_displacement_transient.txt', mean_basal_ab_displacement_transient, delimiter=',')
        np.savetxt(self.results_dir + 'mean_apical_ab_displacement_transient.txt', mean_apical_ab_displacement_transient, delimiter=',')
        np.savetxt(self.results_dir + 'mean_apicobasal_sum_displacement_transient.txt',
                   mean_apicobasal_sum_displacement_transient, delimiter=',')
        plt.figure()
        plt.plot(self.simulation_data['time'], mean_basal_ab_displacement_transient,
                 self.simulation_data['time'], mean_apical_ab_displacement_transient,
                 self.simulation_data['time'], mean_apicobasal_sum_displacement_transient)
        plt.xlabel('Time (s)')
        plt.ylabel('Longitudinal displacement (cm)')
        plt.legend(['Base', 'Apex'])
        plt.savefig(self.results_dir + 'apicobasal_displacement_transient.png')

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
        # mean_torsion_angle_transient = np.zeros(self.simulation_data['time'].shape[0])
        # for time_i in range(self.simulation_data['time'].shape[0]):
        #     updated_apical_slice_coord = self.geometry.nodes_xyz[apical_mapped_nodes, :] + \
        #                                  self.simulation_data['DISPL'][apical_mapped_nodes, :, time_i]
        #     updated_basal_slice_coord = self.geometry.nodes_xyz[basal_slice_nodes, :] + \
        #                                 self.simulation_data['DISPL'][basal_slice_nodes, :, time_i]
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
        # plt.plot(self.simulation_data['time'], mean_torsion_angle_transient)
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
        # lv_wall_thickness_transient = np.zeros(self.simulation_data['time'].shape[0])
        # rv_wall_thickness_transient = np.zeros(self.simulation_data['time'].shape[0])
        # for time_i in range(self.simulation_data['time'].shape[0]):
        #     updated_node_coords = self.geometry.nodes_xyz[:, :] + self.simulation_data['DISPL'][:, :, time_i]
        #     lv_wall_thickness_transient[time_i] = np.mean(np.linalg.norm(updated_node_coords[lv_mapped_epi_nodes,:] -
        #                                                                  updated_node_coords[lv_endo_nodes, :], axis=1))
        #     rv_wall_thickness_transient[time_i] = np.mean(np.linalg.norm(updated_node_coords[rv_mapped_epi_nodes, :] -
        #                                                                  updated_node_coords[rv_endo_nodes, :], axis=1))
        # np.savetxt(self.results_dir + 'lv_wall_thickness_transient.txt', lv_wall_thickness_transient, delimiter=',')
        # np.savetxt(self.results_dir + 'rv_wall_thickness_transient.txt', rv_wall_thickness_transient, delimiter=',')
        # plt.figure()
        # plt.plot(self.simulation_data['time'], lv_wall_thickness_transient, self.simulation_data['time'], rv_wall_thickness_transient)
        # plt.xlabel('Time (s)')
        # plt.ylabel('Averaged wall thickness (cm)')
        # plt.legend(['LV', 'RV'])
        # plt.savefig(self.results_dir + 'wall_thicknesss_transient.png')
        # diastasis_time_idx = np.nonzero(self.simulation_data['time'] == self.simulation_dict['diastasis_t'])
        # self.simulation_biomarkers['mean_lv_thickness_diastasis'] = lv_wall_thickness_transient[diastasis_time_idx]
        # self.simulation_biomarkers['mean_rv_thickness_diastasis'] = rv_wall_thickness_transient[diastasis_time_idx]
        # end_diastole_time_idx = np.nonzero(self.simulation_data['time'] == self.simulation_dict['end_diastole_t'])
        # self.simulation_biomarkers['mean_lv_thickness_end_diastole'] = lv_wall_thickness_transient[end_diastole_time_idx]
        # self.simulation_biomarkers['mean_rv_thickness_end_diastole'] = rv_wall_thickness_transient[end_diastole_time_idx]
        # end_systole_time_idx = np.nonzero(self.simulation_data['lv_phase'] == 3)[0][0] # Index at which phase first changes to 3 IVR.
        # self.simulation_biomarkers['mean_lv_thickness_end_systole'] = lv_wall_thickness_transient[end_systole_time_idx]
        # self.simulation_biomarkers['mean_rv_thickness_end_systole'] = rv_wall_thickness_transient[end_systole_time_idx]
        #
        #
        # print('Fibre stretch ratio transient')
        # lv_lambda_transients = np.mean(self.simulation_data['LAMBD'][lv_nodes, 0, :], axis=0)
        # rv_lambda_transients = np.mean(self.simulation_data['LAMBD'][rv_nodes, 0, :], axis=0)
        # np.savetxt(self.results_dir+'lv_fibre_stretch_ratio_transient.txt', lv_lambda_transients, delimiter=',')
        # np.savetxt(self.results_dir + 'rv_fibre_stretch_ratio_transient.txt', rv_lambda_transients, delimiter=',')
        # plt.figure()
        # plt.plot(self.simulation_data['time'], lv_lambda_transients, self.simulation_data['time'], rv_lambda_transients)
        # plt.xlabel('Time (s)')
        # plt.ylabel('Fibre stretch ratio')
        # plt.legend(['LV', 'RV'])
        # plt.savefig(self.results_dir + 'fibre_stretch_ratio_transient.png')


    def evaluate_strains(self):
        print('Strain evaluations in ventricular coordinates')

def evaluate_lat(time, vm, percentage, time_window):
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
