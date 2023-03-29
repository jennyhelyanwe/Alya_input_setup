from meshstructure import MeshStructure
from myformat import *


class FieldGeneration(MeshStructure):
    def __init__(self, name, geometric_data_dir, personalisation_data_dir, electrode_data_filename, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.personalisation_data_dir = personalisation_data_dir
        # Generate required fields for Alya simulations.
        print('Generate additional Alya fields')
        neighbours, unfolded_edges = evaluate_mesh_characteristics(self.geometry)
        self.node_fields.add_field(data=evaluate_celltype(number_of_nodes=self.geometry.number_of_nodes,
                                                          uvc_transmural=self.node_fields.dict['tm'],
                                                          endo_mid_divide=0.3, mid_epi_divide=0.7),
                                   data_name='cell-type', field_type='nodefield')
        activation_time = np.loadtxt(self.personalisation_data_dir + self.name + '_inferred_lat.csv')
        self.node_fields.add_field(data=activation_time, data_name='activation-time', field_type='nodefield')
        # Read in ionic scaling factors
        filenames = np.array([f for f in os.listdir(self.personalisation_data_dir) if
                              os.path.isfile(
                                  os.path.join(self.personalisation_data_dir, f)) and '.sf_' in f and 'ensi' not in f])
        for file_i in range(filenames.shape[0]):
            if verbose:
                print('Reading in ' + self.personalisation_data_dir + filenames[file_i])
            varname = filenames[file_i].split('.')[1]
            self.node_fields.add_field(data=load_txt(filename=self.personalisation_data_dir + filenames[file_i]),
                                       data_name=varname, field_type='nodefield')
        self.node_fields.add_field(data=evaluate_endocardial_activation_map(activation_times=activation_time,
                                                                            boundary_node_fields=self.boundary_node_fields),
                                   data_name='endocardial-activation-times', field_type='nodefield')
        self.node_fields.add_field(data=load_txt(electrode_data_filename), data_name='electrode_xyz', field_type='nodefield')

        print('Evaluate cavity landmark nodes')
        # LV cavity landmarks
        basal_ring_meta_idx = np.nonzero(self.node_fields.dict['ab'][self.node_fields.dict['ep-lvnodes'].astype(int)] == self.geometry.base)
        basal_ring = self.node_fields.dict['ep-lvnodes'][basal_ring_meta_idx].astype(int)
        rt_posterior = -np.pi/2.
        rt_anterior = np.pi/2.
        posterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_posterior))
        lv_posterior_node = basal_ring[posterior_meta_idx]
        anterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_anterior))
        lv_anterior_node = basal_ring[anterior_meta_idx]

        # RV cavity landmarks
        basal_ring_meta_idx = np.nonzero(self.node_fields.dict['ab'][self.node_fields.dict['ep-rvnodes'].astype(int)] == self.geometry.base)
        basal_ring = self.node_fields.dict['ep-rvnodes'][basal_ring_meta_idx].astype(int)
        posterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_posterior))
        rv_posterior_node = basal_ring[posterior_meta_idx]
        anterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_anterior))
        rv_anterior_node = basal_ring[anterior_meta_idx]
        self.node_fields.add_field(data=np.array([lv_posterior_node, lv_anterior_node]), data_name='lv-cavity-nodes',
                                   field_type='nodefield')
        self.node_fields.add_field(data=np.array([rv_posterior_node, rv_anterior_node]), data_name='rv-cavity-nodes',
                                   field_type='nodefield')
        self.save()


def evaluate_mesh_characteristics(geometry):
    unfolded_edges = np.concatenate((geometry.edges, np.flip(geometry.edges, axis=1))).astype(int)
    aux = [[] for i in range(0, geometry.number_of_nodes, 1)]
    for i in range(0, len(unfolded_edges)):
        aux[unfolded_edges[i, 0]].append(unfolded_edges[i, 1])
    neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
    return neighbours, unfolded_edges


def evaluate_dijkstra_endocardial_activation(number_of_nodes, number_of_faces, face_fields, root_node_locations):
    lv_endocardial_nodes = face_fields


def evaluate_endocardial_activation_map(activation_times, boundary_node_fields):
    print('Evaluating endocardial activation map')
    lv_activation_times = activation_times[boundary_node_fields.dict['ep-lvnodes'].astype(int)]
    rv_activation_times = activation_times[boundary_node_fields.dict['ep-rvnodes'].astype(int)]
    boundary_node_fields.add_field(data=np.concatenate(
        (boundary_node_fields.dict['ep-lvnodes'].astype(int), boundary_node_fields.dict['ep-rvnodes'].astype(int))),
                                   data_name='endocardial-nodes', field_type='nodefield')
    return np.concatenate((lv_activation_times, rv_activation_times))


def evaluate_celltype(number_of_nodes, uvc_transmural, endo_mid_divide, mid_epi_divide):
    print('Evaluating nodal cell type')
    celltype = [1] * number_of_nodes
    for node_i in range(number_of_nodes):
        if uvc_transmural[node_i] <= endo_mid_divide:
            celltype[node_i] = 1  # Endocardial cell type in Alya.
        if uvc_transmural[node_i] > endo_mid_divide and uvc_transmural[
            node_i] < mid_epi_divide:  # It is possible to have no mid-myocardial cells.
            celltype[node_i] = 2  # Midmyocardial cell type in Alya.
        if uvc_transmural[node_i] >= mid_epi_divide:
            celltype[node_i] = 3  # Epicardial cell type in Alya.
    return np.array(celltype)


def evaluate_ab_Gks_scaling(number_of_nodes, uvc_longitudinal, max_sf, min_sf):
    print('Evaluating apex to base gradient in GKs')
    gks_scaling = [0] * number_of_nodes
    for node_i in range(number_of_nodes):
        z_scale = uvc_longitudinal[node_i] * 2.0 - 1.0
        gks_scaling[node_i] = pow(0.2, z_scale)
        if (uvc_longitudinal[node_i] < 0):
            gks_scaling[node_i] = 0
    return np.array(gks_scaling)


def normalise_vector(vector):
    m = np.linalg.norm(vector)
    if m > 0:
        return vector / m
    else:
        return vector


def evaluate_hybrid_rodero_fibres(geometry, node_fields, element_fields, neighbours):
    transmural_vector = node_fields.dict['transmural-vector']
    fibres_nodes = node_fields.dict['fibres']
    number_of_nodes = geometry.number_of_nodes
    ortho = np.zeros([number_of_nodes, 9])
    for i in range(0, number_of_nodes):
        f = fibres_nodes[i, :]
        f = f / np.linalg.norm(f)
        s = transmural_vector[i, :]
        n = np.cross(f, s)
        n = normalise_vector(n)
        ortho[i, :] = list(f) + list(s) + list(n)
    plug_fibres(geometry=geometry, element_fields=element_fields, neighbours=neighbours, ortho=ortho)
    # Sanity check:
    zero_fvector_nodes = np.where(~ortho[:, 0:3].any(axis=1))[0]  # Vectors all zero
    zero_svector_nodes = np.where(~ortho[:, 3:6].any(axis=1))[0]
    zero_nvector_nodes = np.where(~ortho[:, 6:9].any(axis=1))[0]
    nan_vector_nodes = np.where(np.isnan(ortho).any(axis=1))[0]
    print('Ortho calculation sanity check: ')
    print('Number of zero fibre vectors: ', len(zero_fvector_nodes))
    print('Number of zero sheet vectors: ', len(zero_fvector_nodes))
    print('Number of zero normal vectors: ', len(zero_fvector_nodes))
    print('Number of vectors with NaNs: ', len(nan_vector_nodes))
    return ortho


def plug_fibres(geometry, element_fields, neighbours, ortho):
    print('Evaluating valvular plug fibre, sheet, and normal vectors...')
    # Identify all plug elements
    lvrvbase_elems = element_fields.dict['tv-element']
    elems = geometry.tetrahedrons
    nodes = geometry.nodes_xyz
    nnodes = geometry.number_of_nodes
    for k in range(0, 4):
        # Find all elements and nodes in current valve plug
        plug_elems = np.nonzero(lvrvbase_elems == k + 7)[0]  # Plug element is labelled as 7 to 10.
        plug_nodes = np.unique(elems[plug_elems])
        plug_nodes_coords = nodes[plug_nodes - 1, :]

        non_plug_elems = elems[np.nonzero(lvrvbase_elems <= 2)[0]]
        non_plug_nodes = np.unique(non_plug_elems)
        border_elems = np.nonzero(np.any(np.isin(non_plug_elems, plug_nodes), axis=1))
        border_nodes = np.unique(non_plug_elems[border_elems, :])
        ring_nodes = np.array(list(set(border_nodes).difference(set(plug_nodes))), dtype=int)
        ring_nodes_coords = nodes[ring_nodes - 1, :]
        non_plug_nodes = np.array(list(set(non_plug_nodes).difference(set(plug_nodes))), dtype=int)

        # Ordering all plug nodes according to how far they are from the nearest ventricular node
        dists = np.zeros(len(plug_nodes))
        for i in range(0, len(plug_nodes)):
            dists[i] = np.amin(np.linalg.norm(ring_nodes_coords - plug_nodes_coords[i, :], axis=1))
        sorted_plug_nodes = plug_nodes[np.argsort(dists)]
        mask = np.zeros(nnodes)
        mask[non_plug_nodes - 1] = 1
        for i in range(0, len(sorted_plug_nodes)):
            all_closest_nodes = neighbours[sorted_plug_nodes[i] - 1]
            closest_nodes = all_closest_nodes[np.where(mask[all_closest_nodes] > 0)]
            if (len(closest_nodes) == 0):
                print('No viable neighbours found for node number: ', sorted_plug_nodes[i],
                      ' plug number: ' + str(k + 1))
                print('finding next closest nodes...')
                temp = []
                for i in range(0, len(all_closest_nodes)):
                    temp = temp + list(neighbours[all_closest_nodes[i]])
                all_closest_nodes = np.unique(np.array(temp))
                closest_nodes = all_closest_nodes[np.where(mask[all_closest_nodes] > 0)]
            weight = 0.4
            fibres = np.concatenate(
                (ortho[closest_nodes, :], weight * np.array([[1, 0, 0, 0, 1, 0, 0, 0, 1]], dtype=float)))
            new_fibres = np.mean(fibres, axis=0)
            for j in [0, 3, 6]:
                new_fibres[j:j + 3] = normalise_vector(new_fibres[j:j + 3])
            ortho[sorted_plug_nodes[i] - 1, :] = new_fibres
            mask[sorted_plug_nodes[i] - 1] = 1
