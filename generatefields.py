from meshstructure import MeshStructure
from meshpreprocessing import map_indexes
from myformat import *
import pandas as pd
import vtk
from vtk.util import numpy_support as VN
import pymp

class FieldGeneration(MeshStructure):
    def __init__(self, name, geometric_data_dir, personalisation_data_dir, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.personalisation_data_dir = personalisation_data_dir
        self.geometric_data_dir = geometric_data_dir
        self.verbose = verbose

    def generate_celltype(self, endo_mid_divide=0.3, mid_epi_divide=0.7, read_celltype_filename=None, save=False):
        # Generate required fields for Alya simulations.
        if read_celltype_filename:
            print('Reading celltype info from ', read_celltype_filename)
            celltypes = pd.read_csv(read_celltype_filename)['celltype'].values
            celltypes_alya = np.zeros(celltypes.shape).astype(int)
            for i in range(celltypes.shape[0]):
                if celltypes[i] == 'endo':
                    celltypes_alya[i] = 1
                elif celltypes[i] == 'mid':
                    celltypes_alya[i] = 2
                elif celltypes[i] == 'epi':
                    celltypes_alya[i] = 3
            self.node_fields.add_field(data=celltypes_alya, data_name='cell-type', field_type='nodefield')
        else:
            print('Generate additional Alya fields with endo-mid divde ', endo_mid_divide, ' mid-epi divide' , mid_epi_divide)
            self.node_fields.add_field(data=evaluate_celltype(number_of_nodes=self.geometry.number_of_nodes,
                                                              uvc_transmural=self.node_fields.dict['tm'],
                                                              endo_mid_divide=endo_mid_divide, mid_epi_divide=mid_epi_divide),
                                       data_name='cell-type', field_type='nodefield')
        if save:
            self.save()

    def generate_stimuli(self, lat_filename=None, save=False):
        if lat_filename:
            activation_time = np.loadtxt(lat_filename)
        else:
            activation_time = np.loadtxt(self.personalisation_data_dir + self.name + '_nodefield_inferred-lat.csv')
        self.node_fields.add_field(data=activation_time, data_name='activation-time', field_type='nodefield')

        endocardial_activation_times = evaluate_endocardial_activation_map(activation_times=activation_time,
                                                                            boundary_node_fields=self.boundary_node_fields)
        endocardial_activation_times = endocardial_activation_times - np.amin(endocardial_activation_times) # Remove offset
        self.boundary_node_fields.add_field(data=endocardial_activation_times,
                                   data_name='endocardial-activation-times', field_type='boundarynodefield')
        stimulus_node_field = np.zeros(self.geometry.number_of_nodes)
        stimulus_node_field[self.boundary_node_fields.dict['endocardial-nodes']] = endocardial_activation_times
        self.node_fields.add_field(data=stimulus_node_field, data_name='stimulus', field_type='nodefield')
        if save:
            self.save()

    def generate_ionic_scaling_factors(self, read_biomarker_filename=None, save=False):
        if self.verbose:
            print('Generating ionic scaling factor fields...')
        if read_biomarker_filename:
            data = pd.read_csv(read_biomarker_filename)
            keys = pd.read_csv(read_biomarker_filename).keys()
            for key in keys:
                if 'sf_' in key:
                    self.node_fields.add_field(data=data[key].values, data_name=key, field_type='nodefield')
            # varname = read_biomarker_filename.split('.')[1]
            # self.node_fields.add_field(data=load_txt(filename=read_biomarker_filename), data_name=varname, field_type='nodefield')
        else:
            # Read in ionic scaling factors
            filenames = np.array([f for f in os.listdir(self.personalisation_data_dir) if
                                  os.path.isfile(
                                      os.path.join(self.personalisation_data_dir, f)) and '.sf_' in f and 'ensi' not in f])
            for file_i in range(filenames.shape[0]):
                if self.verbose:
                    print('Reading in ' + self.personalisation_data_dir + filenames[file_i])
                varname = filenames[file_i].split('.')[1]
                self.node_fields.add_field(data=load_txt(filename=self.personalisation_data_dir + filenames[file_i]),
                                           data_name=varname, field_type='nodefield')
        if save:
            self.save()

    def generate_orthogonal_local_bases(self, save):
        local_c = pymp.shared.array((self.geometry.number_of_nodes, 3))
        local_l = pymp.shared.array((self.geometry.number_of_nodes, 3))
        local_r = pymp.shared.array((self.geometry.number_of_nodes, 3))
        ln = pymp.shared.array(self.node_fields.dict['longitudinal-vector'].shape)
        tn = pymp.shared.array(self.node_fields.dict['transmural-vector'].shape)
        ln[:,:] = self.node_fields.dict['longitudinal-vector']
        tn[:, :] = self.node_fields.dict['transmural-vector']
        threadsNum = multiprocessing.cpu_count()
        with pymp.Parallel(min(threadsNum, self.geometry.number_of_nodes)) as p1:
            for node_i in p1.range(self.geometry.number_of_nodes):
                l = ln[node_i, :]
                t = tn[node_i, :]
                c = np.cross(l, t)
                r = np.cross(c, l)
                local_c[node_i, :] = c/np.linalg.norm(c)
                local_l[node_i, :] = l/np.linalg.norm(l)
                local_r[node_i, :] = r/np.linalg.norm(r)
        self.node_fields.add_field(data=local_l, data_name='local_l', field_type='nodefield')
        self.node_fields.add_field(data=local_c, data_name='local_c', field_type='nodefield')
        self.node_fields.add_field(data=local_r, data_name='local_r', field_type='nodefield')
        if save:
            self.save()


    def generate_short_long_axes_slices(self, save):
        num_slices = 3
        short_axis_slices = pymp.shared.array((self.geometry.number_of_nodes))
        long_axis_slices = pymp.shared.array((self.geometry.number_of_nodes))
        ab = pymp.shared.array(self.node_fields.dict['ab'].shape)
        rt = pymp.shared.array(self.node_fields.dict['rt'].shape)
        tm = pymp.shared.array(self.node_fields.dict['tm'].shape)
        tv = pymp.shared.array(self.node_fields.dict['tv'].shape)
        ab[:] = self.node_fields.dict['ab']
        rt[:] = self.node_fields.dict['rt']
        tm[:] = self.node_fields.dict['tm']
        tv[:] = self.node_fields.dict['tv']
        threadsNum = multiprocessing.cpu_count()
        with pymp.Parallel(min(threadsNum, self.geometry.number_of_nodes)) as p1:
            for node_i in p1.range(self.geometry.number_of_nodes):
                if (ab[node_i] > 0.27) & (ab[node_i] < 0.3) & (tv[node_i] == -1):
                    short_axis_slices[node_i] = 1 # Apex short-axis
                elif (ab[node_i] > 0.47) & (ab[node_i] < 0.5) & (tv[node_i] == -1):
                    short_axis_slices[node_i] = 2 # Mid short-axis
                elif (ab[node_i] > 0.67) & (ab[node_i] < 0.7) & (tv[node_i] == -1):
                    short_axis_slices[node_i] = 3 # Base short-axis
                else:
                    short_axis_slices[node_i] = 0

                if (rt[node_i] > 1.75) & (rt[node_i] < 1.95) & (tv[node_i] == -1):
                    long_axis_slices[node_i] = 1 # 2-chamber long-axis
                elif (rt[node_i] > -1.6) & (rt[node_i] < -1.4)  & (tv[node_i] == -1):
                    long_axis_slices[node_i] = 1 # 2-chamber long-axis
                elif (rt[node_i] > -2.5) & (rt[node_i] < -2.3)  & (tv[node_i] == -1):
                    long_axis_slices[node_i] = 2 # 4-chamber long-axis
                elif (rt[node_i] > 0.9) & (rt[node_i] < 1.1)  & (tv[node_i] == -1):
                    long_axis_slices[node_i] = 2 # 4-chamber long-axis
                elif (rt[node_i] > 2.8) & (rt[node_i] < 3.2)  & (tv[node_i] == -1):
                    long_axis_slices[node_i] = 3  # 3-chamber long-axis
                elif (rt[node_i] > -0.1) & (rt[node_i] < 0.3)  & (tv[node_i] == -1):
                    long_axis_slices[node_i] = 3 # 3-chamber long-axis
                else:
                    long_axis_slices[node_i] = 0
        self.node_fields.add_field(data=long_axis_slices, data_name='long-axis-slices', field_type='nodefield')
        self.node_fields.add_field(data=short_axis_slices, data_name='short-axis-slices', field_type='nodefield')
        if save:
            self.save()

    def generate_electrode_locations(self, electrode_data_filename, save=False):
        self.node_fields.add_field(data=load_txt(electrode_data_filename), data_name='electrode_xyz', field_type='nodefield')
        if save:
            self.save()

    def generate_cavity_landmarks(self):
        print('Evaluate cavity landmark nodes')
        # LV cavity landmarks
        basal_ring_meta_idx = np.nonzero(self.node_fields.dict['ab'][self.boundary_node_fields.dict['ep-lvnodes'].astype(int)] == self.geometry.base)
        basal_ring = self.boundary_node_fields.dict['ep-lvnodes'][basal_ring_meta_idx].astype(int)
        rt_posterior = -np.pi/2.
        rt_anterior = np.pi/2.
        posterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_posterior))
        lv_posterior_node = basal_ring[posterior_meta_idx]
        anterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_anterior))
        lv_anterior_node = basal_ring[anterior_meta_idx]

        # RV cavity landmarks
        basal_ring_meta_idx = np.nonzero(self.node_fields.dict['ab'][self.boundary_node_fields.dict['ep-rvnodes'].astype(int)] == self.geometry.base)
        basal_ring = self.boundary_node_fields.dict['ep-rvnodes'][basal_ring_meta_idx].astype(int)
        posterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_posterior))
        rv_posterior_node = basal_ring[posterior_meta_idx]
        anterior_meta_idx = np.argmin(abs(self.node_fields.dict['rt'][basal_ring] - rt_anterior))
        rv_anterior_node = basal_ring[anterior_meta_idx]
        self.node_fields.add_field(data=np.array([lv_posterior_node, lv_anterior_node]), data_name='lv-cavity-nodes',
                                   field_type='nodefield')
        self.node_fields.add_field(data=np.array([rv_posterior_node, rv_anterior_node]), data_name='rv-cavity-nodes',
                                   field_type='nodefield')
        self.save()

    def generate_prestress_field(self):
        # Prestress field
        prestress_field = np.zeros(self.element_fields.dict['tv-element'].shape[0]).astype(int)
        for element_i in range(self.element_fields.dict['tv-element'].shape[0]):
            if (self.element_fields.dict['tv-element'][element_i] == 7) or \
                    (self.element_fields.dict['tv-element'][element_i] == 9) or \
                    (self.element_fields.dict['tv-element'][element_i] == 1):
                prestress_field[element_i] = 1
            elif self.element_fields.dict['tv-element'][element_i] == 8 or \
                    (self.element_fields.dict['tv-element'][element_i] == 10) or \
                    (self.element_fields.dict['tv-element'][element_i] == 2):
                prestress_field[element_i] = 2
        self.element_fields.add_field(data=prestress_field, data_name='prestress', field_type='elementfield')
        self.save()

    def read_fibre_fields(self, fibre_field_filename, sheet_field_filename, normal_field_filename):
        self.node_fields.add_field(data=pd.read_csv(fibre_field_filename, header=None).values, data_name='fibre', field_type='nodefield')
        self.node_fields.add_field(data=pd.read_csv(sheet_field_filename, header=None).values, data_name='sheet',
                                   field_type='nodefield')
        self.node_fields.add_field(data=pd.read_csv(normal_field_filename, header=None).values, data_name='normal',
                                   field_type='nodefield')
        self.save()

    def read_doste_fibre_fields_vtk(self, fibre_vtk_filename, sheet_vtk_filename, normal_vtk_filename):
        if self.verbose:
            print('Reading in Doste fibres, ensuring node correspondence using nearest node mapping...')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(fibre_vtk_filename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        doste_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())
        # Nearest node mapping
        map = map_indexes(points_to_map_xyz=self.geometry.nodes_xyz, reference_points_xyz=doste_xyz)

        # Read fibres
        if self.verbose:
            print('Reading in fibre field...')
        fibre = VN.vtk_to_numpy(data.GetPointData().GetArray('Fibers'))
        # Check for zero vectors:
        n_zero_vectors = len(np.where(~fibre.any(axis=1))[0])
        n_nan_vectors = len(np.where(np.isnan(fibre).any(axis=1))[0])
        print('Number of zero fibre vectors: ', n_zero_vectors)
        print('Number of NaN fibre vectors: ', n_nan_vectors)
        self.node_fields.add_field(data=fibre[map], data_name='fibre', field_type='nodefield')

        # Read sheets
        if self.verbose:
            print('Reading in sheet field...')
        reader.SetFileName(sheet_vtk_filename)
        reader.Update()
        data = reader.GetOutput()
        sheet = VN.vtk_to_numpy(data.GetPointData().GetArray('Fibers'))
        # Check for zero vectors:
        n_zero_vectors = len(np.where(~sheet.any(axis=1))[0])
        n_nan_vectors = len(np.where(np.isnan(sheet).any(axis=1))[0])
        print('Number of zero sheet vectors: ', n_zero_vectors)
        print('Number of NaN sheet vectors: ', n_nan_vectors)
        self.node_fields.add_field(data=sheet[map], data_name='sheet', field_type='nodefield')

        # Read normal
        if self.verbose:
            print('Reading in normal field...')
        reader.SetFileName(normal_vtk_filename)
        reader.Update()
        data = reader.GetOutput()
        normal = VN.vtk_to_numpy(data.GetPointData().GetArray('Fibers'))
        # Check for zero vectors:
        n_zero_vectors = len(np.where(~normal.any(axis=1))[0])
        n_nan_vectors = len(np.where(np.isnan(normal).any(axis=1))[0])
        print('Number of zero normal vectors: ', n_zero_vectors)
        print('Number of NaN normal vectors: ', n_nan_vectors)
        self.node_fields.add_field(data=normal[map], data_name='normal', field_type='nodefield')
        self.save()

    def generate_infarct_borderzone(self):
        # Using UVC to define infarct and border zone geometry.

        # Evaluate ab, rt, and tm at tetrahedron centres
        # if not 'ab' in list(self.element_fields.dict.keys()):
        #     if self.verbose:
        #         print('Mapping UVC to element field...')
        #     map = map_indexes(points_to_map_xyz=self.geometry.tetrahedron_centres, reference_points_xyz=self.geometry.nodes_xyz)
        #     ab_tetra_centres = self.node_fields.dict['ab'][map]
        #     rt_tetra_centres = self.node_fields.dict['rt'][map]
        #     tm_tetra_centres = self.node_fields.dict['tm'][map]
        #     tv_tetra_centres = self.node_fields.dict['tv'][map]
        #     self.element_fields.add_field(data=ab_tetra_centres, data_name='ab_tetra', field_type='elementfield')
        #     self.element_fields.add_field(data=rt_tetra_centres, data_name='rt_tetra', field_type='elementfield')
        #     self.element_fields.add_field(data=tm_tetra_centres, data_name='tm_tetra', field_type='elementfield')
        #     self.element_fields.add_field(data=tv_tetra_centres, data_name='tv_tetra', field_type='elementfield')
        #     self.save()
        # quit()
        if self.verbose:
            print('Delineating infarct and borderzones...')
        materials = self.materials.dict['tetra']
        materials_shared = pymp.shared.array((materials.shape[0]), dtype=int)
        ab_shared = pymp.shared.array((materials.shape[0]), dtype=float)
        rt_shared = pymp.shared.array((materials.shape[0]), dtype=float)
        tm_shared = pymp.shared.array((materials.shape[0]), dtype=float)
        materials_shared[:] = self.materials.dict['tetra'][:]
        ab_shared[:] = self.element_fields.dict['ab_tetra'][:]
        rt_shared[:] = self.element_fields.dict['rt_tetra'][:]
        tm_shared[:] = self.element_fields.dict['tm_tetra'][:]
        threadsNum = multiprocessing.cpu_count()
        infarct_rt_range = [1.2, 2.1]
        infarct_ab_range = [0.1, 0.5]
        infarct_ab_centre = np.abs((infarct_ab_range[1] + infarct_ab_range[0]) / 2.)
        infarct_rt_centre = np.abs((infarct_rt_range[1] + infarct_rt_range[0]) / 2.)
        infarct_a = infarct_ab_range[1] - infarct_ab_centre
        infarct_b = infarct_rt_range[1] - infarct_rt_centre
        bz_rt_range = [1.0, 2.4]
        bz_ab_range = [0.0, 0.6]
        bz_ab_centre = np.abs((bz_ab_range[1] + bz_ab_range[0]) / 2.)
        bz_rt_centre = np.abs((bz_rt_range[1] + bz_rt_range[0]) / 2.)
        bz_a = bz_ab_range[1] - bz_ab_centre
        bz_b = bz_rt_range[1] - bz_rt_centre
        with pymp.Parallel(min(threadsNum, self.geometry.tetrahedron_centres.shape[0])) as p1:
            for i in p1.range(self.geometry.tetrahedron_centres.shape[0]):
        # if True:
        #     for i in range(self.geometry.tetrahedron_centres.shape[0]):
                bz_equation = ((ab_shared[i] - bz_ab_centre) / bz_a) ** 2 + ((rt_shared[i] - bz_rt_centre) / bz_b) ** 2 - 1
                if (bz_equation < 0):
                    materials_shared[i] = 4
                infarct_equation = ((ab_shared[i] - infarct_ab_centre) / infarct_a) ** 2 + ((rt_shared[i] - infarct_rt_centre) / infarct_b) ** 2 - 1
                if (infarct_equation < 0) & (tm_shared[i] < 0.75):
                    materials_shared[i] = 3 # Infarct
        self.materials.add_field(data=materials_shared, data_name='tetra_mi', field_type='material')
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
    if not 'endocardial-nodes' in boundary_node_fields.dict:
        lv_activation_times = activation_times[boundary_node_fields.dict['ep-lvnodes'].astype(int)]
        rv_activation_times = activation_times[boundary_node_fields.dict['ep-rvnodes'].astype(int)]
        boundary_node_fields.add_field(data=np.concatenate(
            (boundary_node_fields.dict['ep-lvnodes'].astype(int), boundary_node_fields.dict['ep-rvnodes'].astype(int))),
                                       data_name='endocardial-nodes', field_type='boundarynodefield')
        return np.concatenate((lv_activation_times, rv_activation_times))
    else:
        boundary_node_fields.dict['endocardial-nodes'] = boundary_node_fields.dict['endocardial-nodes'].astype(int)
        return activation_times[boundary_node_fields.dict['endocardial-nodes'].astype(int)]


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

