import numpy as np
import os
import json
# from mpi4py import MPI #needs install
# from multiprocessing import Queue
import pymp, multiprocessing

from meshstructure import MeshStructure
from myformat import Fields


class PostProcessing(MeshStructure):
    def __init__(self, name, geometric_data_dir, simulation_json_file, alyacsv_dir, results_dir, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.alyacsv_dir = alyacsv_dir
        self.results_dir = results_dir
        self.simulation_biomarkers = {}
        self.simulation_data = {}
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.post_nodefield = Fields(name, field_type='postnodefield', verbose=verbose)
        self.post_elementfield = Fields(name, field_type='postelementfield', verbose=verbose)
        print('Reading postprocessing csv fields')
        self.post_nodefield.read_csv_to_attributes(input_dir=results_dir, field_type='postnodefield')
        # self.post_elementfield.read_csv_to_attributes(input_dir=results_dir, field_type='postelementfield')
        # Evaluate various biomarkers
        self.read_csv_fields()
        self.evaluate_ep_maps()
        self.evaluate_ventricular_cvs()

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

    def read_csv_fields(self):
        print('Reading Vm')
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
            temp = pymp.shared.array((self.geometry.number_of_nodes, time.shape[0]), dtype=int)
            threadsNum = multiprocessing.cpu_count()
            with pymp.Parallel(min(threadsNum, time_index.shape[0])) as p1:
                for time_i in p1.range(time_index.shape[0]):
                    index = '{:06d}'.format(time_index[time_i])
                    filename = self.alyacsv_dir + name + '.' + field_type + '.' + field_name + '-' + index + '.csv'
                    temp[:, time_i] = np.loadtxt(filename, delimiter=',').astype(float)
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
