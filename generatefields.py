from meshstructure import MeshStructure
from meshpreprocessing import map_indexes
from myformat import *
import pandas as pd
import vtk
from vtk.util import numpy_support as VN
import pymp
from matplotlib import pyplot as plt
import numba
from dijkstra import Dijkstra

class FieldGeneration(MeshStructure):
    def __init__(self, name, geometric_data_dir, personalisation_data_dir, max_cores_used, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, max_cores_used=max_cores_used, verbose=verbose)
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
        threadsNum = np.amin((multiprocessing.cpu_count(), self.max_cores_used))
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
        threadsNum = np.amin((multiprocessing.cpu_count(), self.max_cores_used))
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

    def generate_cavity_landmarks_three_nodes(self):
        print('Evaluate cavity landmark nodes (old CC version, 3 nodes needed)')
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

    def read_cavity_landmarks(self, save=False):
        lv_posterior_node = 5078
        lv_anterior_node = 5151
        self.node_fields.add_field(data=np.array([lv_posterior_node, lv_anterior_node]), data_name='lv-cavity-nodes',
                                   field_type='nodefield')
        if save:
            self.save()

    def generate_prestress_field(self, save=False):
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
        if save:
            self.save()

    def generate_prestress_field_cardiax(self, save=False):
        prestress_field = np.zeros(self.geometry.number_of_elements).astype(int)
        prestress_field[:] = 1
        self.element_fields.add_field(data=prestress_field, data_name='prestress', field_type='elementfield')
        if save:
            self.save()

    def read_fibre_fields(self, fibre_field_filename, sheet_field_filename, normal_field_filename, save=False):
        self.node_fields.add_field(data=pd.read_csv(fibre_field_filename, header=None).values, data_name='fibre', field_type='nodefield')
        self.node_fields.add_field(data=pd.read_csv(sheet_field_filename, header=None).values, data_name='sheet',
                                   field_type='nodefield')
        self.node_fields.add_field(data=pd.read_csv(normal_field_filename, header=None).values, data_name='normal',
                                   field_type='nodefield')
        if save:
            self.save()

    def read_cardiax_ellipsoid_fibre_fields_vtk(self, fibre_vtk_filename, save=False):
        if self.verbose:
            print('Reading in fibres, ensuring node correspondence using nearest node mapping...')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(fibre_vtk_filename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        fibre_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())
        # Nearest node mapping
        map = map_indexes(points_to_map_xyz=self.geometry.nodes_xyz, reference_points_xyz=fibre_xyz)

        # Read fibres
        if self.verbose:
            print('Reading in fibre field...')
        fibre = VN.vtk_to_numpy(data.GetPointData().GetArray('f0'))
        fibre = fibre[map]

        # Randomly generate sheet and normal vectors
        sheet = np.zeros(fibre.shape)
        normal = np.zeros(fibre.shape)
        for i in range(sheet.shape[0]):
            f = fibre[i, :]
            random = [0, 1, 0]
            if (np.dot(f, random) < 0.00001):
                random = [1, 0, 0]
            s = np.cross(random, f)
            s = s/np.linalg.norm(s)
            n = np.cross(f, s)
            n = n/np.linalg.norm(n)
            sheet[i, :] = s
            normal[i, :] = n

        # Check for zero vectors:
        n_zero_vectors = len(np.where(~fibre.any(axis=1))[0])
        n_nan_vectors = len(np.where(np.isnan(fibre).any(axis=1))[0])
        print('Number of zero fibre vectors: ', n_zero_vectors)
        print('Number of NaN fibre vectors: ', n_nan_vectors)
        n_zero_vectors = len(np.where(~sheet.any(axis=1))[0])
        n_nan_vectors = len(np.where(np.isnan(sheet).any(axis=1))[0])
        print('Number of zero fibre vectors: ', n_zero_vectors)
        print('Number of NaN fibre vectors: ', n_nan_vectors)
        n_zero_vectors = len(np.where(~normal.any(axis=1))[0])
        n_nan_vectors = len(np.where(np.isnan(normal).any(axis=1))[0])
        print('Number of zero normal vectors: ', n_zero_vectors)
        print('Number of NaN normal vectors: ', n_nan_vectors)
        self.node_fields.add_field(data=fibre, data_name='fibre', field_type='nodefield')
        self.node_fields.add_field(data=sheet, data_name='sheet', field_type='nodefield')
        self.node_fields.add_field(data=normal, data_name='normal', field_type='nodefield')
        if save:
            self.save()

    def map_doste_nodes_to_rodero_nodes(self, fibre_vtk_filename, map_filename):
        if os.path.exists(map_filename):
            map = np.loadtxt(map_filename).astype(int)
        else:
            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(fibre_vtk_filename)
            reader.ReadAllVectorsOn()
            reader.ReadAllScalarsOn()
            reader.Update()
            data = reader.GetOutput()
            doste_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())
            # Nearest node mapping
            map = map_indexes(points_to_map_xyz=doste_xyz, reference_points_xyz=self.geometry.nodes_xyz)
            print(map.shape)
            print(self.geometry.nodes_xyz.shape)
            np.savetxt(map_filename, map)
        return map

    def read_doste_fibre_fields_vtk(self, fibre_vtk_filename, sheet_vtk_filename, normal_vtk_filename, map, save=False):
        if self.verbose:
            print('Reading in Doste fibres, ensuring node correspondence using nearest node mapping...')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(fibre_vtk_filename)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        doste_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())

        # Read fibres
        if self.verbose:
            print('Reading in fibre field...')
        fibre = VN.vtk_to_numpy(data.GetPointData().GetArray('Fibers'))
        print('VTK fibre: ', fibre)
        # Check for zero vectors:
        n_zero_vectors = len(np.where(~fibre.any(axis=1))[0])
        n_nan_vectors = len(np.where(np.isnan(fibre).any(axis=1))[0])
        print('Number of zero fibre vectors: ', n_zero_vectors)
        print('Number of NaN fibre vectors: ', n_nan_vectors)
        if 'fibre' in self.node_fields.dict.keys():
            self.node_fields.dict['fibre'] = fibre[map,:]
        else:
            self.node_fields.add_field(data=fibre[map,:], data_name='fibre', field_type='nodefield')

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
        if 'sheet' in self.node_fields.dict.keys():
            self.node_fields.dict['sheet'] = sheet[map,:]
        else:
            self.node_fields.add_field(data=sheet[map,:], data_name='sheet', field_type='nodefield')

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
        if 'normal' in self.node_fields.dict.keys():
            self.node_fields.dict['normal'] = normal[map,:]
        else:
            self.node_fields.add_field(data=normal[map,:], data_name='normal', field_type='nodefield')
        if save:
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

        if self.verbose:
            print('Delineating infarct and borderzones...')
        materials = self.materials.dict['tetra']
        materials_shared = pymp.shared.array(materials.shape, dtype=int)
        cell_centres = pymp.shared.array(self.geometry.tetrahedron_centres.shape, dtype=float)
        tm_shared = pymp.shared.array((self.node_fields.dict['tm'].shape[0]), dtype=float)
        materials_shared[:] = materials[:]
        cell_centres[:,:] = self.geometry.tetrahedron_centres[:, :]
        tm_shared[:] = self.node_fields.dict['tm'][:]
        threadsNum = np.amin((multiprocessing.cpu_count(), 10))
        infarct_centres = np.array([[-1.91, -2.977, 6.066], [-1.198, -2.075, 3.98], [-2.40, -1.36, 2.88], [-1.58, -3.91, 4.48],
                           [-1.56, 0.187, 1.65]])
        infarct_radii = [2, 2, 2, 1, 1]
        bz_width = 0.5
        with (pymp.Parallel(min(threadsNum, cell_centres.shape[0]))) as p1:
            for conf_i in p1.range(0, cell_centres.shape[0]):
                for i_centre in range(infarct_centres.shape[0]):
                    infarct_centre = infarct_centres[i_centre,:]
                    infarct_radius = infarct_radii[i_centre]
                    bz_radius = infarct_radius + bz_width
                    infarct_equation = np.sqrt((cell_centres[conf_i,0]-infarct_centre[0])**2+
                                               (cell_centres[conf_i,1]-infarct_centre[1])**2+
                                               (cell_centres[conf_i,2]-infarct_centre[2])**2)
                    if (infarct_equation <= infarct_radius) & (tm_shared[self.geometry.tetrahedrons[conf_i,0]] < 0.75):
                        materials_shared[conf_i] = 3
                    elif (infarct_equation<bz_radius) & (not (materials_shared[conf_i]==3)):
                        materials_shared[conf_i] = 4
        self.materials.add_field(data=materials_shared, data_name='tetra_mi', field_type='material')
        self.save()
        # infarct_centres = [[1.0, 1.0]] # rt and ab
        # infarct_rt_range = [1.2, 2.1]
        # infarct_ab_range = [0.1, 0.5]
        # infarct_ab_centre = np.abs((infarct_ab_range[1] + infarct_ab_range[0]) / 2.)
        # infarct_rt_centre = np.abs((infarct_rt_range[1] + infarct_rt_range[0]) / 2.)
        # infarct_a = infarct_ab_range[1] - infarct_ab_centre
        # infarct_b = infarct_rt_range[1] - infarct_rt_centre
        # bz_rt_range = [1.0, 2.4]
        # bz_ab_range = [0.0, 0.6]
        # bz_ab_centre = np.abs((bz_ab_range[1] + bz_ab_range[0]) / 2.)
        # bz_rt_centre = np.abs((bz_rt_range[1] + bz_rt_range[0]) / 2.)
        # bz_a = bz_ab_range[1] - bz_ab_centre
        # bz_b = bz_rt_range[1] - bz_rt_centre
        # Define centres and ranges for ab and rt coordinates

    def generate_infarct_borderzone_Wallman(self, visualise=False):
        print('Generating infarct and borderzone using Wallman code...')
        dijkstra = Dijkstra()
        # Re-implementation of MATLAB code from Wallman, based on algorithm as described in https://doi.org/10.1371/journal.pone.0149342 and
        # https://ieeexplore.ieee.org/document/7043132) and selects elements as infarcted
        # Select region for scar and BZ
        # LAD:
        scar_candidates = []
        LAD_ab = 0.4
        LAD_tm = 1
        LAD_rt = 1.61
        LAD_tv = -1 # LV
        TOL = 0.01
        lv_select = np.where((self.node_fields.dict['tv'] == LAD_tv))[0]
        ab_select = lv_select[np.where(abs(self.node_fields.dict['ab'][lv_select] - LAD_ab) < TOL)[0]]
        tm_select = ab_select[np.where((self.node_fields.dict['tm'][ab_select] == LAD_tm))[0]]
        rt_select = tm_select[np.where(abs(self.node_fields.dict['rt'][tm_select] - LAD_rt) < TOL)[0]]
        scar_candidates.append(rt_select[0])

        LAD_ab = 0.6
        LAD_tm = 1
        LAD_rt = 1.61
        LAD_tv = -1  # LV
        TOL = 0.01
        lv_select = np.where((self.node_fields.dict['tv'] == LAD_tv))[0]
        ab_select = lv_select[np.where(abs(self.node_fields.dict['ab'][lv_select] - LAD_ab) < TOL)[0]]
        tm_select = ab_select[np.where((self.node_fields.dict['tm'][ab_select] == LAD_tm))[0]]
        rt_select = tm_select[np.where(abs(self.node_fields.dict['rt'][tm_select] - LAD_rt) < TOL)[0]]
        scar_candidates.append(rt_select[0])

        LAD_ab = 0.4
        LAD_tm = 1
        LAD_rt = 1.05
        LAD_tv = -1  # LV
        TOL = 0.01
        lv_select = np.where((self.node_fields.dict['tv'] == LAD_tv))[0]
        ab_select = lv_select[np.where(abs(self.node_fields.dict['ab'][lv_select] - LAD_ab) < TOL)[0]]
        tm_select = ab_select[np.where((self.node_fields.dict['tm'][ab_select] == LAD_tm))[0]]
        rt_select = tm_select[np.where(abs(self.node_fields.dict['rt'][tm_select] - LAD_rt) < TOL)[0]]
        scar_candidates.append(rt_select[0])
        seed_times = len(scar_candidates) * [0]

        LAD_ab = 0.5
        LAD_tm = 1
        LAD_rt = 1.08
        LAD_tv = -1  # LV
        TOL = 0.01
        lv_select = np.where((self.node_fields.dict['tv'] == LAD_tv))[0]
        ab_select = lv_select[np.where(abs(self.node_fields.dict['ab'][lv_select] - LAD_ab) < TOL)[0]]
        tm_select = ab_select[np.where((self.node_fields.dict['tm'][ab_select] == LAD_tm))[0]]
        rt_select = tm_select[np.where(abs(self.node_fields.dict['rt'][tm_select] - LAD_rt) < TOL)[0]]
        scar_candidates.append(rt_select[0])

        seed_times = len(scar_candidates) * [0]
        print(scar_candidates)
        # Generate core and border zones
        print ('Generate scar region...')
        long_axis = self.node_fields.dict['longitudinal-vector'].mean(axis=0)
        long_axis = long_axis/np.linalg.norm(long_axis)
        projection = np.matmul(self.geometry.nodes_xyz, long_axis)
        heart_length = np.max(projection) - np.min(projection)
        cutoff = (1.2 - self.node_fields.dict['ab']) * (heart_length) * 0.5
        cutoff[np.where(self.node_fields.dict['ab'] == -10)[0]] = np.nan

        approx_djikstra_max_path_len = 100
        neighbour, nodes_xyz, unfoldedEdges, edgeVEC = dijkstra.prepare_for_dijkstra(edge=self.geometry.edges,
                                                                                     node_xyz=self.geometry.nodes_xyz)
        num_iterations = 2
        for j in range(num_iterations):
            distance = dijkstra.iso_eikonal(sub_node_coordinates=nodes_xyz, source_indexes=scar_candidates,
                                       seed_times=seed_times, sub_edge_indexes=self.geometry.edges,
                                       sub_edgeVEC=edgeVEC, sub_neighbour=neighbour, sub_unfoldedEdges=unfoldedEdges)
        # distance = dijkstra.sorted_dijkstra(source_indexes=np.asarray(scar_candidates, dtype=int), seed_times=seed_times,
        #                                   dijkstra_nodes_xyz=self.geometry.nodes_xyz,
        #                                   dijkstra_unfoldedEdges=unfoldedEdges,
        #                                   dijkstra_edgeVEC=edgeVEC,
        #                                   dijkstra_neighbours=neighbour,
        #                                   approx_dijkstra_max_path_len=approx_djikstra_max_path_len)
        scar = np.where(distance < cutoff)[0]
        expr = np.zeros(self.geometry.number_of_nodes)
        expr[scar] = 1

        print ('Generate core and border zones...')
        core_ext_noise = 6
        core_post_noise = 0
        ncores = 250 # Dense scar cores
        core = scar[np.random.uniform(0, len(scar), ncores).astype(int)]
        seed_times = abs(np.random.normal(0, core_ext_noise, len(core)))
        num_iterations = 2
        for j in range(num_iterations):
            core_d = dijkstra.iso_eikonal(sub_node_coordinates=nodes_xyz, source_indexes=np.asarray(core, dtype=int),
                                         seed_times=seed_times, sub_edge_indexes=self.geometry.edges,
                                         sub_edgeVEC=edgeVEC, sub_neighbour=neighbour, sub_unfoldedEdges=unfoldedEdges)
        # core_d = dijkstra.sorted_dijkstra(source_indexes=np.asarray(core, dtype=int), seed_times=seed_times,
        #                                   dijkstra_nodes_xyz=self.geometry.nodes_xyz,
        #                                   dijkstra_unfoldedEdges=unfoldedEdges,
        #                                   dijkstra_edgeVEC=edgeVEC,
        #                                   dijkstra_neighbours=neighbour,
        #                                   approx_dijkstra_max_path_len=approx_djikstra_max_path_len)
        core_d = core_d/np.amax(core_d) + np.random.normal(0, core_post_noise, len(core_d))
        cz_threshold = 0.13
        bz_threshold = 0.18
        cz_nodes = np.where(core_d < cz_threshold)[0]
        bz_nodes = np.where(core_d < bz_threshold)[0]
        cz_bz_regions = np.zeros(self.geometry.number_of_nodes)
        cz_bz_regions[bz_nodes] = 2
        cz_bz_regions[cz_nodes] = 1
        if visualise:
            print('Visualising...')
            def scatter_visualise(ax, xyz, field, title):
                p = ax.scatter(xyz[:, 0], xyz[:, 1], xyz[:, 2], c=field, marker='o', s=1)
                ax.set_title(title)
                return p
            node_downsample_skip = int(self.geometry.number_of_nodes / 10000.)
            node_xyz = self.geometry.nodes_xyz[::node_downsample_skip, :]
            d_downsample = distance[::node_downsample_skip]

            fig = plt.figure(figsize=(16, 12))
            ax1 = fig.add_subplot(131, projection='3d')
            p = scatter_visualise(ax1, node_xyz,d_downsample,'Dijkstra solution')
            ax1.scatter(self.geometry.nodes_xyz[scar_candidates, 0], self.geometry.nodes_xyz[scar_candidates, 1], self.geometry.nodes_xyz[scar_candidates, 2], marker='o', s=3, c='r')
            # ax1.quiver(0,0,0, long_axis[0], long_axis[1], long_axis[2], length=1, color='r' )
            # ax1.quiver(0, 0, 0, rvlv_axis[0], rvlv_axis[1], rvlv_axis[2], length=1, color='b')
            # ax1.quiver(0, 0, 0, slanted_axis[0], slanted_axis[1], slanted_axis[2], length=1, color='m')
            fig.colorbar(p)

            ax2 = fig.add_subplot(132, projection='3d')
            rvlv_downsample = self.node_fields.dict['rvlv'][::node_downsample_skip]
            scar_region = expr[::node_downsample_skip]
            # slanted_coordinate_downsample = slanted_coordinate[::node_downsample_skip]
            cutoff_downsample = cutoff[::node_downsample_skip]
            p = scatter_visualise(ax2, node_xyz, scar_region, 'Scar region')
            fig.colorbar(p)

            ax3 = fig.add_subplot(133, projection='3d')
            cz_bz_region_downsample = cz_bz_regions[::node_downsample_skip]
            core_d_downsample = core_d[::node_downsample_skip]
            p = scatter_visualise(ax3, node_xyz, cz_bz_region_downsample, 'Cores and Border Zones')
            fig.colorbar(p)
            plt.show()
        quit()

        # Save as new material fields.
        materials = []
        self.materials.add_field(data=materials, data_name='tetra_mi', field_type='material')
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




def dijkstra_specify_seed_node_times(nodes_xyz, edge, navigation_costs, seed_nodes, seed_times):
    # Adapted from Julia Camps' implementation of Dijkstra from https://github.com/juliacamps/Cardiac-Digital-Twin/blob/main/src/propagation_models.py
    adjacent_cost = numba.typed.List()
    unfolded_edge = np.concatenate((edge, np.flip(edge, axis=1))).astype(int)
    aux = [[] for i in range(0, np.shape(nodes_xyz)[0], 1)]
    for edge_i in range(0, len(unfolded_edge)):
        aux[unfolded_edge[edge_i, 0]].append(unfolded_edge[edge_i, 1])
    neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
    for i in range(0, nodes_xyz.shape[0], 1):
        not_nan_neighbours = neighbours[i][neighbours[i] != np.nan]
        adjacent_cost.append(np.concatenate((unfolded_edge[not_nan_neighbours][:, 1:2],
                                             np.expand_dims(navigation_costs[not_nan_neighbours % edge.shape[0]], -1)), axis=1))
    seed_nodes = np.array(seed_nodes).astype(np.int32)
    seed_times = np.array(seed_times)
    predicted_lat = np.zeros((nodes_xyz.shape[0],), np.float64)
    visited_nodes = np.zeros((nodes_xyz.shape[0],), dtype=np.bool_)
    temp_times = np.zeros((nodes_xyz.shape[0],), np.float64) + 1e6  # Initialise times to largely impossible values
    temp_times[seed_nodes] = seed_times
    time_sorting = np.array(np.argsort(seed_times))
    seed_nodes = seed_nodes[time_sorting]
    seed_times = seed_times[time_sorting]
    cumm_cost = seed_times[0]
    initial_root_nodes_indexes = seed_times <= cumm_cost
    initial_rootNodes = seed_nodes[initial_root_nodes_indexes]
    initial_rootActivationTimes = seed_times[initial_root_nodes_indexes]
    later_rootNodes = seed_nodes[np.logical_not(initial_root_nodes_indexes)]
    later_rootActivationTimes = seed_times[np.logical_not(initial_root_nodes_indexes)]
    visited_nodes[initial_rootNodes] = True  # Not simultaneous activation anymore
    predicted_lat[initial_rootNodes] = initial_rootActivationTimes  # Not simultaneous activation anymore
    next_nodes = (np.vstack([adjacent_cost[initial_rootNodes[rootNode_i]]
                             + np.array([0, initial_rootActivationTimes[rootNode_i]]) for rootNode_i in
                             range(initial_rootNodes.shape[0])])).tolist()  # Not simultaneous activation anymore
    for rootNode_i in range(later_rootNodes.shape[0]):
        next_nodes.append(np.array([later_rootNodes[rootNode_i], later_rootActivationTimes[rootNode_i]]))
    activeNode_i = seed_nodes[0]
    sortSecond = lambda x: x[1]
    next_nodes.sort(key=sortSecond, reverse=True)

    while visited_nodes[activeNode_i]:
        nextEdge = next_nodes.pop()
        activeNode_i = int(nextEdge[0])
    cumm_cost = nextEdge[1]
    if next_nodes:  # Check if the list is empty, which can happen while everything being Ok
        temp_times[(np.array(next_nodes)[:, 0]).astype(np.int32)] = np.array(next_nodes)[:,1]
    def insert_sorted(aList, newV):
        # Function to insert element
        ini_index = 0
        end_index = len(aList)
        index = int((end_index - ini_index) / 2)
        for i in range(0, len(aList), 1):
            if newV[1] < aList[index][1]:
                if end_index - ini_index <= 1 or index + 1 == end_index:
                    index = index + 1
                    break
                else:
                    ini_index = index + 1
                    index = int(index + (end_index - ini_index) / 2 + 1)
            elif newV[1] > aList[index][1]:
                if end_index - ini_index <= 1 or index == ini_index:
                    index = index  # Place before the current position
                    break
                else:
                    end_index = index
                    index = int(index - (end_index - ini_index) / 2)
            else:
                index = ini_index
                break
        aList.insert(index, newV)

    ## Run the whole algorithm
    for i in range(0, nodes_xyz.shape[0] - np.sum(visited_nodes), 1):
        visited_nodes[activeNode_i] = True
        predicted_lat[activeNode_i] = cumm_cost  # Instead of using cumCost, I could use the actual time cost for each node
        adjacents = (adjacent_cost[activeNode_i] + np.array([0, cumm_cost])).tolist()  # Instead of using cumCost, I could use the actual time cost for each node
        # If I use the actual costs, I only have to do it the first time and then it will just propagate, I will have to use decimals though, so no more uint type arrays.
        for adjacent_i in range(0, len(adjacents), 1):
            if (not visited_nodes[int(adjacents[adjacent_i][0])]
                    and (temp_times[int(adjacents[adjacent_i][0])] >
                         adjacents[adjacent_i][1])):
                insert_sorted(next_nodes, adjacents[adjacent_i])
                temp_times[int(adjacents[adjacent_i][0])] = adjacents[adjacent_i][1]
        while visited_nodes[activeNode_i] and len(next_nodes) > 0:
            nextEdge = next_nodes.pop()
            activeNode_i = int(nextEdge[0])
        cumm_cost = nextEdge[1]

    # Clean Memory
    adjacent_cost = None  # Clear Mem
    visited_nodes = None  # Clear Mem
    tempTimes = None  # Clear Mem
    next_nodes = None  # Clear Mem
    tempVisited = None  # Clear Mem
    navigationCosts = None  # Clear Mem
    return predicted_lat
