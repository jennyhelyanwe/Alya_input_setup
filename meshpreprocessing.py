import vtk
from vtk.util import numpy_support as VN

from meshstructure import MeshStructure
from myformat import *
import copy

class MeshPreprocessing(MeshStructure):
    def __init__(self, vtk_name, name, input_dir, geometric_data_dir, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        # Read and write geometry
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.vtk_name = vtk_name
        self.input_dir = input_dir
        self.read_geometry_from_vtk_rodero()
        self.lv_endocardium = 2
        self.rv_endocardium = 3
        self.epicardium = 1
        self.valve_plug = 4
        # Generate fields
        self.generate_boundary_data_rodero()
        self.generate_fibre_sheet_normal()
        self.generate_additional_vc()

        self.save()
        self.check_fields_for_qrs_inference()
        self.check_fields_for_twave_inference()

    def check_fields_for_qrs_inference(self):
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_xyz.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_tetra.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_edges.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_ab.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_tm.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_rt.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_tv.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_material_tetra.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_fibre-sheet-normal.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_ep-lvnodes.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_ep-rvnodes.csv')
        print('Everything set for QRS inference')

    def check_fields_for_twave_inference(self):
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_xyz.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_tetra.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_edges.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_ab.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_tm.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_rt.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_tv.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_aprt.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_nodefield_rvlv.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_material_tetra.csv')
        print('Everything set for T wave inference')

    def read_geometry_from_vtk_rodero(self):
        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '.vtk')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        self.geometry.number_of_elements = data.GetNumberOfCells()
        self.geometry.number_of_nodes = data.GetNumberOfPoints()
        self.geometry.tetrahedrons = np.zeros((self.geometry.number_of_elements, 4)).astype(int)
        self.geometry.nodes_xyz = np.zeros((self.geometry.number_of_nodes, 3))
        for i in range(0, self.geometry.number_of_elements):
            self.geometry.tetrahedrons[i, 0] = int(data.GetCell(i).GetPointId(0))
            self.geometry.tetrahedrons[i, 1] = int(data.GetCell(i).GetPointId(1))
            self.geometry.tetrahedrons[i, 2] = int(data.GetCell(i).GetPointId(2))
            self.geometry.tetrahedrons[i, 3] = int(data.GetCell(i).GetPointId(3))
        nodes_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())  # Convert from mm to cm
        self.geometry.nodes_xyz = nodes_xyz
        self.geometry.edges = []
        for element_i in range(self.geometry.number_of_elements):
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 0], self.geometry.tetrahedrons[element_i, 1]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 1], self.geometry.tetrahedrons[element_i, 2]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 2], self.geometry.tetrahedrons[element_i, 3]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 3], self.geometry.tetrahedrons[element_i, 0]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 1], self.geometry.tetrahedrons[element_i, 3]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 0], self.geometry.tetrahedrons[element_i, 2]])
        self.geometry.edges = np.unique(np.sort(self.geometry.edges, axis=1), axis=0)
        self.geometry.tetrahedron_centres = (nodes_xyz[self.geometry.tetrahedrons[:, 0], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 1], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 2], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 3], :]) / 4.
        # self.element_fields.lvrv = VN.vtk_to_numpy(data.GetCellData().GetArray('ID'))
        self.element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('ID')), 'tv-element', 'elementfield')
        self.element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('fibres')), 'fibres-element',
                                      'elementfield')
        self.element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('sheets')), 'sheets-element',
                                      'elementfield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('V.dat')), 'tv', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RHO.dat')), 'tm', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('Z.dat')), 'ab', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('PHI.dat')), 'rt', 'nodefield')
        print('Convert from UVC to Cobiveco')
        self.convert_uvc_to_cobiveco()
        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_fibres_at_nodes.vtk')
        reader.SetFileName(self.input_dir + self.vtk_name + '_fibres_at_nodes.vtk')
        reader.Update()
        data = reader.GetOutput()
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('fibres')), 'fibres', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('sheets')), 'sheets', 'nodefield')

        def _normalise_vector(vector):
            m = np.linalg.norm(vector)
            if m > 0:
                return vector / m
            else:
                return vector

        print('Reading vtk file from: ' + self.input_dir + 'transmural_vectors.vtk')
        reader.SetFileName(self.input_dir + 'transmural_vectors.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        transmural_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(transmural_vector.shape[0]):
            transmural_vector[node_i, :] = _normalise_vector(transmural_vector[node_i, :])
        self.node_fields.add_field(transmural_vector, 'transmural-vector', 'nodefield')
        print('Reading vtk file from: ' + self.input_dir + 'longitudinal_vectors.vtk')
        reader.SetFileName(self.input_dir + 'longitudinal_vectors.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        longitudinal_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(longitudinal_vector.shape[0]):
            longitudinal_vector[node_i, :] = _normalise_vector(longitudinal_vector[node_i, :])
        self.node_fields.add_field(longitudinal_vector, 'longitudinal-vector', 'nodefield')
        print('Reading vtk file from: ' + self.input_dir + 'circumferential_vectors.vtk')
        reader.SetFileName(self.input_dir + 'circumferential_vectors.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        circumferential_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(circumferential_vector.shape[0]):
            circumferential_vector[node_i, :] = _normalise_vector(circumferential_vector[node_i, :])
        self.node_fields.add_field(circumferential_vector, 'circumferential-vector', 'nodefield')

        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()

        # Read faces
        self.geometry.number_of_triangles = data.GetNumberOfCells()
        self.geometry.triangles = np.zeros([self.geometry.number_of_triangles, 3])
        self.boundary_node_fields.add_field(data=VN.vtk_to_numpy(data.GetPointData().GetArray('Ids')).astype(int),
                                            data_name='surface-node-id', field_type='nodefield')
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(0)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(1)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)
        self.geometry.save_to_csv(self.geometric_data_dir)
        self.node_fields.save_to_csv(self.geometric_data_dir)
        self.element_fields.save_to_csv(self.geometric_data_dir)
        self.geometry.save_to_ensight(self.geometric_data_dir + 'ensight/')
        self.node_fields.save_to_ensight(self.geometric_data_dir + 'ensight/', casename=self.name + '_nodefield',
                                         geometry=self.geometry)
        self.element_fields.save_to_ensight(self.geometric_data_dir + 'ensight/', casename=self.name + '_elementfield',
                                            geometry=self.geometry)
        materials = np.zeros(self.element_fields.dict['tv-element'].shape[0]).astype(int)
        for element_i in range(self.element_fields.dict['tv-element'].shape[0]):
            if self.element_fields.dict['tv-element'][element_i] > 3:
                materials[element_i] = 2
            else:
                materials[element_i] = 1
        self.materials.add_field(data=materials, data_name='tetra', field_type='material')
        self.materials.save_to_csv(self.geometric_data_dir)
        self.materials.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/', casename=self.name,
                                       geometry=self.geometry)

    def convert_uvc_to_cobiveco(self):
        # Map input geometry Cobiveco coordinates to output-style UVC
        cobiveco = {}
        cobiveco['ab'] = copy.deepcopy(self.node_fields.dict['ab'])
        # Handle mapping of TV split separately:
        cobiveco['tv'] = map_ventricular_coordinates(input_ranges=[[-1, 1]], output_ranges=[[0, 1]],
                                                     input_coordinate=self.node_fields.dict['tv'])
        cobiveco_rv_septum_indices = []
        for node_i in range(cobiveco['ab'].shape[0]):
            if (self.node_fields.dict['rt'][node_i] > -1.66) & (self.node_fields.dict['rt'][node_i] < 1.19) & (self.node_fields.dict['tv'][node_i] == -1) & (
                    self.node_fields.dict['tm'][node_i] > 0.6):
                cobiveco['tv'][node_i] = 1
                cobiveco_rv_septum_indices.append(node_i)
        cobiveco_lv_indices = np.where(cobiveco['tv'] == 0)
        cobiveco_rv_indices = np.where(cobiveco['tv'] == 1)
        # Translation of RT is ventricle dependent
        cobiveco['rt'] = copy.deepcopy(self.node_fields.dict['rt'])
        cobiveco['rt'][cobiveco_lv_indices] = map_ventricular_coordinates(
            input_ranges=[[-np.pi, -1.3], [-1.3, 0], [0, np.pi, ]],
            # output_ranges=[[0.4, 0], [1, 0.8], [0.8, 0.4]],
            output_ranges=[[0.3, 0], [1, 0.8], [0.8, 0.3]],
            input_coordinate=self.node_fields.dict['rt'][cobiveco_lv_indices])
        # output_ranges=[[0.4, 0], [1, 0.8], [0.8, 0.4]],
        cobiveco['rt'][cobiveco_rv_indices] = map_ventricular_coordinates(
            input_ranges=[[0, 1.2], [-1.3, 0]],
            output_ranges=[[0.4, 0.75], [0, 0.4]],
            input_coordinate=self.node_fields.dict['rt'][cobiveco_rv_indices])  # RT
        cobiveco['rt'][cobiveco_rv_septum_indices] = map_ventricular_coordinates(
            input_ranges=[[-1.3, 1.19]],
            output_ranges=[[1, 0.6]], input_coordinate=self.node_fields.dict['rt'][cobiveco_rv_septum_indices])
        cobiveco['tm'] = map_ventricular_coordinates(input_ranges=[[0, 1]], output_ranges=[[1, 0]],
                                                     input_coordinate=self.node_fields.dict['tm'])  # TM
        for key in cobiveco.keys():
            self.node_fields.add_field(data=cobiveco[key], data_name='cobiveco_'+key, field_type='nodefield')


    def generate_boundary_data_rodero(self):
        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        elem_ids = VN.vtk_to_numpy(data.GetCellData().GetArray('Ids'))
        self.boundary_element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId')).astype(int),
                                               'boundary-label', 'elementfield')
        self.boundary_node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RegionId')).astype(int),
                                            'boundary-label', 'nodefield')

        for i in range(0, len(self.boundary_element_fields.dict['boundary-label'])):
            if self.boundary_element_fields.dict['boundary-label'][i] == 0:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.epicardium  # Epicardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == 1:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium  # LV endocardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == 2:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.rv_endocardium  # RV endocardial

        for i in range(0, len(self.boundary_node_fields.dict['boundary-label'])):
            if self.boundary_node_fields.dict['boundary-label'][i] == 0:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.epicardium  # Epicardial
            elif self.boundary_node_fields.dict['boundary-label'][i] == 1:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium  # LV endocardial
            elif self.boundary_node_fields.dict['boundary-label'][i] == 2:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.rv_endocardium  # RV endocardial
        # Save input global label
        input_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['boundary-label'] == self.geometry.lv_endocardium]] = self.geometry.lv_endocardium
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['boundary-label'] == self.geometry.rv_endocardium]] = self.geometry.rv_endocardium
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['boundary-label'] == self.geometry.epicardium]] = self.geometry.epicardium
        self.boundary_node_fields.add_field(input_node_label_global, 'input-boundary-label', 'nodefield')

        # Set up boundary label for EP and mechanics
        mechanical_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        ep_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        mechanical_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        ep_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        for i in range(0, len(self.geometry.triangles)):
            for j in range(0, 3):
                local_index = \
                np.nonzero(self.boundary_node_fields.dict['surface-node-id'] == self.geometry.triangles[i, j])[0][0]
                if (self.boundary_node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & (
                        self.node_fields.dict['tm'][local_index] == -10):
                    mechanical_node_boundary_label[local_index] = self.valve_plug
                    ep_node_boundary_label[local_index] = self.valve_plug
                    mechanical_element_boundary_label[i] = self.valve_plug
                else:
                    mechanical_element_boundary_label[i] = \
                    self.boundary_element_fields.dict['boundary-label'][local_index]
                    ep_element_boundary_label[local_index] = self.boundary_element_fields.dict['boundary-label'][
                        local_index]
                    mechanical_node_boundary_label[local_index] = self.boundary_node_fields.dict['boundary-label'][
                        local_index]
                    ep_node_boundary_label[local_index] = self.boundary_node_fields.dict['boundary-label'][local_index]
                if self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10:
                    ep_node_boundary_label[local_index] = self.valve_plug
                # Restrict mechanical epicardium to only 80% of apex-to-base axis.
                if (self.boundary_node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & \
                        (self.node_fields.dict['ab'][self.geometry.triangles[i, j]] > 0.8):
                    # if (self.        self.boundary_element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId')).astype(int),.dict['boundary_label'][i] == self.epicardium) & \
                    #         (self.node_fields.dict['ab'][self.geometry.triangles[i, j]] > 0.8):
                    mechanical_element_boundary_label[
                        i] = self.valve_plug  # Apply epicardial spring BC only to 80% of apex-to-base
                    mechanical_node_boundary_label[local_index] = self.valve_plug
        self.boundary_node_fields.add_field(mechanical_node_boundary_label, 'mechanical-node-boundary-label',
                                            'nodefield')
        self.boundary_node_fields.add_field(ep_node_boundary_label, 'ep-node-boundary-label', 'nodefield')
        self.boundary_element_fields.add_field(mechanical_element_boundary_label, 'mechanical-element-boundary-label',
                                               'elementfield')
        self.boundary_element_fields.add_field(ep_element_boundary_label, 'ep-element-boundary-label', 'elementfield')

        # Get LV and RV endocardium nodes
        ep_lvnodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['ep-node-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        ep_rvnodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['ep-node-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        ep_epinodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['ep-node-boundary-label'] == self.geometry.epicardium].astype(int)
        ep_valvenodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['ep-node-boundary-label'] == self.geometry.valve_plug].astype(int)
        ep_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        ep_node_label_global[ep_lvnodes] = self.geometry.lv_endocardium
        ep_node_label_global[ep_rvnodes] = self.geometry.rv_endocardium
        ep_node_label_global[ep_epinodes] = self.geometry.epicardium
        ep_node_label_global[ep_valvenodes] = self.geometry.valve_plug
        self.boundary_node_fields.add_field(ep_node_label_global, 'ep-node-label-global', 'nodefield')
        self.boundary_node_fields.add_field(ep_lvnodes, 'ep-lvnodes', 'nodefield')
        self.boundary_node_fields.add_field(ep_rvnodes, 'ep-rvnodes', 'nodefield')
        self.boundary_node_fields.add_field(ep_epinodes, 'ep-epinodes', 'nodefield')
        self.boundary_node_fields.add_field(ep_valvenodes, 'ep-valvenodes', 'nodefield')

        mechanical_lvnodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        mechanical_rvnodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        mechanical_epinodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.epicardium].astype(int)
        mechanical_valvenodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.valve_plug].astype(int)
        mechanical_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        mechanical_node_label_global[mechanical_lvnodes] = self.geometry.lv_endocardium
        mechanical_node_label_global[mechanical_rvnodes] = self.geometry.rv_endocardium
        mechanical_node_label_global[mechanical_epinodes] = self.geometry.epicardium
        mechanical_node_label_global[mechanical_valvenodes] = self.geometry.valve_plug
        self.boundary_node_fields.add_field(mechanical_node_label_global, 'mechanical-node-label-global', 'nodefield')
        self.boundary_node_fields.add_field(mechanical_lvnodes, 'mechanical-lvnodes', 'nodefield')
        self.boundary_node_fields.add_field(mechanical_rvnodes, 'mechanical-rvnodes', 'nodefield')
        self.boundary_node_fields.add_field(mechanical_epinodes, 'mechanical-epinodes', 'nodefield')
        self.boundary_node_fields.add_field(mechanical_valvenodes, 'mechanical-valvenodes', 'nodefield')

        # Get LV and RV endocardium faces for Eikonal solution.
        ep_lvfaces = self.geometry.triangles[
            self.boundary_element_fields.dict['ep-element-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        ep_rvfaces = self.geometry.triangles[
            self.boundary_element_fields.dict['ep-element-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        self.boundary_element_fields.add_field(ep_lvfaces, 'ep-lvfaces', 'elementfield')
        self.boundary_element_fields.add_field(ep_rvfaces, 'ep-rvfaces', 'elementfield')

        self.boundary_node_fields.save_to_csv(self.geometric_data_dir)
        self.boundary_element_fields.save_to_csv(self.geometric_data_dir)
        self.geometry.save_to_ensight(self.geometric_data_dir + 'ensight/')
        self.boundary_node_fields.save_to_ensight(self.geometric_data_dir + 'ensight/',
                                                  casename=self.name + '_boundary_node_fields', geometry=self.geometry)
        self.boundary_element_fields.save_to_ensight(self.geometric_data_dir + 'ensight/',
                                                     casename=self.name + '_boundary_element_fields',
                                                     geometry=self.geometry)

    def generate_fibre_sheet_normal(self):
        unfolded_edges = np.concatenate((self.geometry.edges, np.flip(self.geometry.edges, axis=1))).astype(int)
        aux = [[] for i in range(0, self.geometry.number_of_nodes, 1)]
        for i in range(0, len(unfolded_edges)):
            aux[unfolded_edges[i, 0]].append(unfolded_edges[i, 1])
        neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
        transmural_vector = self.node_fields.dict['transmural-vector']
        fibres_nodes = self.node_fields.dict['fibres']
        number_of_nodes = self.geometry.number_of_nodes
        ortho = np.zeros([number_of_nodes, 9])
        for i in range(0, number_of_nodes):
            f = fibres_nodes[i, :]
            f = f / np.linalg.norm(f)
            s = transmural_vector[i, :]
            n = np.cross(f, s)
            n = normalise_vector(n)
            ortho[i, :] = list(f) + list(s) + list(n)
        self.plug_fibres(neighbours=neighbours, ortho=ortho)
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
        self.node_fields.add_field(data=ortho, data_name='fibre-sheet-normal', field_type='nodefield')
        self.node_fields.add_field(data=ortho[:, 0:3], data_name='fibre', field_type='nodefield')
        self.node_fields.add_field(data=ortho[:, 3:6], data_name='sheet', field_type='nodefield')
        self.node_fields.add_field(data=ortho[:, 6:9], data_name='normal', field_type='nodefield')

    def plug_fibres(self, neighbours, ortho):
        print('Evaluating valvular plug fibre, sheet, and normal vectors...')
        # Identify all plug elements
        lvrvbase_elems = self.element_fields.dict['tv-element']
        elems = self.geometry.tetrahedrons
        nodes = self.geometry.nodes_xyz
        nnodes = self.geometry.number_of_nodes
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
                if len(closest_nodes) == 0:
                    print('No viable neighbours found for node number: ', sorted_plug_nodes[i],
                          ' plug number: ' + str(k + 1))
                    print('finding next closest nodes...')
                    temp = []
                    for node_i in range(0, len(all_closest_nodes)):
                        temp = temp + list(neighbours[all_closest_nodes[node_i]])
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

    def generate_additional_vc(self):
        rt = self.node_fields.dict['cobiveco_rt']
        tv = self.node_fields.dict['cobiveco_tv']
        aprt = np.zeros(rt.shape)
        for i in range(rt.shape[0]):
            if (rt[i] >= 1 / 8) & (rt[i] <= 5 / 8):
                aprt[i] = 2 * rt[i] - 1 / 4
            elif rt[i] >= 5 / 8:
                aprt[i] = -20 / 8 * rt[i] + 21 / 8  # (10 / 8 - rt[i]) / (1/2)
            elif rt[i] <= 1 / 8:
                aprt[i] = -rt[i] + 1 / 8  # (1/ 8 - rt[i]) / (1/2)
        # RVLV rotational coordinate - rvlv
        rvlv = np.zeros(rt.shape)
        for i in range(rt.shape[0]):
            if tv[i] == 1:
                rvlv[i] = 1  # Don't apply gradient on RV.
            elif rt[i] <= 3 / 8:
                rvlv[i] = -8 / 3 * rt[i] + 1
            elif (rt[i] >= 3 / 8) & (rt[i] <= 2 / 3):
                rvlv[i] = 24 / 7 * rt[i] - 9 / 7
            elif rt[i] >= 2 / 3:
                rvlv[i] = 1
        self.node_fields.add_field(data=aprt, data_name='cobiveco_aprt', field_type='nodefield')
        self.node_fields.add_field(data=rvlv, data_name='cobiveco_rvlv', field_type='nodefield')


def normalise_vector(vector):
    m = np.linalg.norm(vector)
    if m > 0:
        return vector / m
    else:
        return vector


def map_ventricular_coordinates(input_ranges, output_ranges, input_coordinate):
    input_ranges = np.array(input_ranges)
    output_ranges = np.array(output_ranges)
    assert input_ranges.shape[0] == output_ranges.shape[0], 'Input and output ranges need to pair up for conversion'
    output_coordinate = np.zeros(input_coordinate.shape)
    for segment_i in range(input_ranges.shape[0]):
        a = (output_ranges[segment_i, 1] - output_ranges[segment_i, 0]) / (
                    input_ranges[segment_i, 1] - input_ranges[segment_i, 0])
        b = output_ranges[segment_i, 0] - a * input_ranges[segment_i, 0]
        for coord_i in range(output_coordinate.shape[0]):
            if (input_coordinate[coord_i] >= input_ranges[segment_i, 0]) & (
                    input_coordinate[coord_i] <= input_ranges[segment_i, 1]):
                output_coordinate[coord_i] = a * input_coordinate[coord_i] + b
    return output_coordinate
