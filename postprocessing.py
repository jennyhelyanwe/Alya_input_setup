import numpy as np
import os
import json
# from mpi4py import MPI #needs install
# from multiprocessing import Queue

from meshstructure import MeshStructure


class PostProcessing:
    def __init__(self, name, geometric_data_dir, simulation_json_file, alyacsv_dir, results_dir, verbose):
        # super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.alyacsv_dir = alyacsv_dir
        self.results_dir = results_dir
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.read_csv()


    def read_csv(self):
        name = self.simulation_dict['name']
        self.time = np.loadtxt(self.alyacsv_dir+'timeset_1.csv', delimiter=',')
        self.time_index = np.loadtxt(self.alyacsv_dir+'timeindex_1.csv', delimiter=',').astype(int)
        assert self.time.shape == self.time_index.shape
        filenames = np.array([f for f in os.listdir(self.alyacsv_dir) if
                              os.path.isfile(os.path.join(self.alyacsv_dir, f)) and ('scalar' in f or 'vector' in f)])
        field_names_types = []
        for filename in filenames:
            field_names_types.append([filename.split('.')[2].split('-')[0], filename.split('.')[1]])
        field_names_types = np.unique(field_names_types, axis=0)
        self.data = {}
        for field_i in range(field_names_types.shape[0]):
            temp = []
            field_name = field_names_types[field_i, 0]
            field_type = field_names_types[field_i, 1]
            for time_i in range(self.time_index.shape[0]):
                index = '{:06d}'.format(self.time_index[time_i])
                filename= self.alyacsv_dir + name + '.' + field_type + '.' + field_name+'-'+index+'.csv'
                temp.append(list(np.loadtxt(filename, delimiter=',')))
            self.data[field_name] = np.array(temp)

    def evaluate_lat(self, percentage, time_window):
        vm = self.data['INTRA']
        time = self.time
        window_idx = np.nonzero(time > time_window[0] & time < time_window[1])[0]
        vm = vm[window_idx]
        time = time[window_idx] - time_window[0]  # Offset by beginning of time window.
        activation_map = np.zeros((vm.shape[0]))
        vm_range = np.amax(vm, axis=1) - np.amin(vm, axis=1)
        vm_threshold = vm_range * (1.0 - percentage / 100.0) + np.amin(vm, axis=1)
        for node_i in range(vm.shape[0]):  # Loop through every node in mesh
            local_vm = vm[node_i, :]
            index = np.nonzero(local_vm > vm_threshold[node_i])[0][0]  # np.searchsorted(local_vm, vm_threshold[node_i])
            activation_map[node_i] = time[index]
        return activation_map

    def evaluate_rt(self, percentage, time_window):
        vm = self.data['INTRA']
        time = self.time
        window_idx = np.nonzero(time > time_window[0] & time < time_window[1])[0]
        vm = vm[window_idx]
        time = time[window_idx] - time_window[0] # Offset by beginning of time window.
        repolarisation_map = np.ones((vm.shape[0]))
        vm_range = np.amax(vm, axis=1) - np.amin(vm, axis=1)
        vm_threshold = vm_range * (1.0 - percentage / 100.0) + np.amin(vm, axis=1)
        vm_max_idx = np.argmax(vm, axis=1)
        for node_i in range(vm.shape[0]):  # Loop through every node in mesh
            local_vm = vm[node_i, vm_max_idx[node_i]:]
            local_vm_fliped = np.flip(local_vm)
            fliped_index = np.searchsorted(local_vm_fliped, vm_threshold[node_i])
            index = local_vm.shape[0] - fliped_index - 1
            repolarisation_map[node_i] = time[index + vm_max_idx[node_i]]
        return repolarisation_map
