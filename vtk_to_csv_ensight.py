import vtk
import numpy as np
from myformat import *
from vtk.util import numpy_support as VN


class MeshPreprocessing:
    def __init__(self, vtk_name, output_name, input_dir, output_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        self.vtk_name = vtk_name
        self.output_name = output_name
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.geometry = Geometry(self.output_name)
        self.node_fields = Fields(self.output_name, field_type='nodefield')
        self.element_fields = Fields(self.output_name, field_type='elementfield')
        self.read_vtk()
        self.write_to_csv()
        self.write_to_ensight()

    def read_vtk(self):
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
        self.geometry.lv_endocardium = 3
        self.geometry.rv_endocardium = 2
        self.geometry.epicardium = 1
        self.geometry.valve_plug = 4
        for i in range(0, self.geometry.number_of_elements):
            self.geometry.tetrahedrons[i, 0] = int(data.GetCell(i).GetPointId(0))
            self.geometry.tetrahedrons[i, 1] = int(data.GetCell(i).GetPointId(1))
            self.geometry.tetrahedrons[i, 2] = int(data.GetCell(i).GetPointId(2))
            self.geometry.tetrahedrons[i, 3] = int(data.GetCell(i).GetPointId(3))
        nodes_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())  # Convert from mm to cm
        self.geometry.nodes_xyz = nodes_xyz
        self.geometry.edges = []
        for element_i in range(self.geometry.number_of_elements):
            self.geometry.edges.append([self.geometry.tetrahedrons[element_i, 0], self.geometry.tetrahedrons[element_i, 1]])
            self.geometry.edges.append([self.geometry.tetrahedrons[element_i, 1], self.geometry.tetrahedrons[element_i, 2]])
            self.geometry.edges.append([self.geometry.tetrahedrons[element_i, 2], self.geometry.tetrahedrons[element_i, 3]])
            self.geometry.edges.append([self.geometry.tetrahedrons[element_i, 3], self.geometry.tetrahedrons[element_i, 0]])
            self.geometry.edges.append([self.geometry.tetrahedrons[element_i, 1], self.geometry.tetrahedrons[element_i, 3]])
            self.geometry.edges.append([self.geometry.tetrahedrons[element_i, 0], self.geometry.tetrahedrons[element_i, 2]])
        self.geometry.edges = np.unique(np.sort(self.geometry.edges, axis=1), axis=0)
        self.geometry.tetrahedron_centres = (nodes_xyz[self.geometry.tetrahedrons[:, 0], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 1], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 2], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 3], :])/4.
        # self.element_fields.lvrv = VN.vtk_to_numpy(data.GetCellData().GetArray('ID'))
        self.element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('ID')), 'tv-element', 'elementfield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('fibres')), 'fibres', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('sheets')), 'sheets', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('V.dat')), 'tv', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RHO.dat')), 'tm', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('Z.dat')), 'ab', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('PHI.dat')), 'rt', 'nodefield')

        def _normalise_vector(vector):
            m = np.linalg.norm(vector)
            if m > 0:
                return vector / m
            else:
                return vector
        print('Reading vtk file from: '+ self.input_dir + 'transmural_vectors.vtk')
        reader.SetFileName(self.input_dir + 'transmural_vectors.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        transmural_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(transmural_vector.shape[0]):
            transmural_vector[node_i, :] = _normalise_vector(transmural_vector[node_i,:])
        self.node_fields.add_field(transmural_vector, 'transmural-vector', 'nodefield')
        print('Reading vtk file from: ' + self.input_dir + 'longitudinal_vectors.vtk')
        reader.SetFileName(self.input_dir + 'longitudinal_vectors.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        longitudinal_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(longitudinal_vector.shape[0]):
            longitudinal_vector[node_i, :] = _normalise_vector(longitudinal_vector[node_i,:])
        self.node_fields.add_field(longitudinal_vector, 'longitudinal-vector', 'nodefield')
        print('Reading vtk file from: ' + self.input_dir + 'circumferential_vectors.vtk')
        reader.SetFileName(self.input_dir + 'circumferential_vectors.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        circumferential_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(circumferential_vector.shape[0]):
            circumferential_vector[node_i, :] = _normalise_vector(circumferential_vector[node_i,:])
        self.node_fields.add_field(circumferential_vector, 'circumferential-vector', 'nodefield')

        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        # elem_ids = VN.vtk_to_numpy(data.GetCellData().GetArray('Ids'))
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('Ids')).astype(int), 'surface-node-id', 'nodefield')
        self.element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId')).astype(int), 'boundary-label', 'elementfield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RegionId')).astype(int), 'boundary-label', 'nodefield')

        for i in range(0, len(self.element_fields.dict['boundary-label'])):
            if self.element_fields.dict['boundary-label'][i] == 0:
                self.element_fields.dict['boundary-label'][i] = self.geometry.epicardium # Epicardial
            elif self.element_fields.dict['boundary-label'][i] == 1:
                self.element_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium # LV endocardial
            elif self.element_fields.dict['boundary-label'][i] == 2:
                self.element_fields.dict['boundary-label'][i] = self.geometry.rv_endocardium # RV endocardial

        for i in range(0, len(self.node_fields.dict['boundary-label'])):
            if self.node_fields.dict['boundary-label'][i] == 0:
                self.node_fields.dict['boundary-label'][i] = self.geometry.epicardium # Epicardial
            elif self.node_fields.dict['boundary-label'][i] == 1:
                self.node_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium # LV endocardial
            elif self.node_fields.dict['boundary-label'][i] == 2:
                self.node_fields.dict['boundary-label'][i] = self.geometry.rv_endocardium # RV endocardial
        # Save input global label
        input_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        input_node_label_global[self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['boundary-label'] == self.geometry.lv_endocardium]] = self.geometry.lv_endocardium
        input_node_label_global[self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['boundary-label'] == self.geometry.rv_endocardium]] = self.geometry.rv_endocardium
        input_node_label_global[self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['boundary-label'] == self.geometry.epicardium]] = self.geometry.epicardium
        self.node_fields.add_field(input_node_label_global, 'input-boundary-label', 'nodefield')

        # Read faces
        self.geometry.number_of_triangles = data.GetNumberOfCells()
        self.geometry.triangles = np.zeros([self.geometry.number_of_triangles, 3])
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [int(self.node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(0)]),
                                             int(self.node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(1)]),
                                             int(self.node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)

        # Set up boundary label for EP and mechanics
        mechanical_element_boundary_label = np.zeros(self.element_fields.dict['boundary-label'].shape)
        ep_element_boundary_label = np.zeros(self.element_fields.dict['boundary-label'].shape)
        mechanical_node_boundary_label = np.zeros(self.node_fields.dict['boundary-label'].shape)
        ep_node_boundary_label = np.zeros(self.node_fields.dict['boundary-label'].shape)
        for i in range(0, len(self.geometry.triangles)):
            for j in range(0, 3):
                local_index = np.nonzero(self.node_fields.dict['surface-node-id'] == self.geometry.triangles[i, j])[
                    0][0]
                if (self.node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & (self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10):
                    mechanical_node_boundary_label[local_index] = self.geometry.valve_plug
                    ep_node_boundary_label[local_index] = self.geometry.valve_plug
                    mechanical_element_boundary_label[i] = self.geometry.valve_plug
                else:
                    mechanical_element_boundary_label[local_index] = self.element_fields.dict['boundary-label'][local_index]
                    ep_element_boundary_label[local_index] = self.element_fields.dict['boundary-label'][local_index]
                    mechanical_node_boundary_label[local_index] = self.node_fields.dict['boundary-label'][
                        local_index]
                    ep_node_boundary_label[local_index] = self.node_fields.dict['boundary-label'][local_index]
                if self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10:
                    ep_node_boundary_label[local_index] = self.geometry.valve_plug
                # Restrict mechanical epicardium to only 80% of apex-to-base axis.
                if (self.node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & \
                                (self.node_fields.dict['ab'][self.geometry.triangles[i, j]] > 0.8):
                # if (self.element_fields.dict['boundary_label'][i] == self.geometry.epicardium) & \
                #         (self.node_fields.dict['ab'][self.geometry.triangles[i, j]] > 0.8):
                    mechanical_element_boundary_label[i] = self.geometry.valve_plug  # Apply epicardial spring BC only to 80% of apex-to-base
                    mechanical_node_boundary_label[local_index] = self.geometry.valve_plug
        self.node_fields.add_field(mechanical_node_boundary_label, 'mechanical-node-boundary-label', 'nodefield')
        self.node_fields.add_field(ep_node_boundary_label, 'ep-node-boundary-label', 'nodefield')
        self.element_fields.add_field(mechanical_element_boundary_label, 'mechanical-element-boundary-label', 'elementfield')
        self.element_fields.add_field(ep_element_boundary_label, 'ep-element-boundary-label', 'elementfield')

        # Get LV and RV endocardium nodes
        ep_lvnodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['ep-node-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        ep_rvnodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['ep-node-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        ep_epinodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['ep-node-boundary-label'] == self.geometry.epicardium].astype(int)
        ep_valvenodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['ep-node-boundary-label'] == self.geometry.valve_plug].astype(int)
        ep_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        ep_node_label_global[ep_lvnodes] = self.geometry.lv_endocardium
        ep_node_label_global[ep_rvnodes] = self.geometry.rv_endocardium
        ep_node_label_global[ep_epinodes] = self.geometry.epicardium
        ep_node_label_global[ep_valvenodes] = self.geometry.valve_plug
        self.node_fields.add_field(ep_node_label_global, 'ep-node-label-global', 'nodefield')
        self.node_fields.add_field(ep_lvnodes, 'ep-lvnodes', 'nodefield')
        self.node_fields.add_field(ep_rvnodes, 'ep-rvnodes', 'nodefield')
        self.node_fields.add_field(ep_epinodes, 'ep-epinodes', 'nodefield')
        self.node_fields.add_field(ep_valvenodes, 'ep-valvenodes', 'nodefield')

        mechanical_lvnodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['mechanical-node-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        mechanical_rvnodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['mechanical-node-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        mechanical_epinodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['mechanical-node-boundary-label'] == self.geometry.epicardium].astype(int)
        mechanical_valvenodes = self.node_fields.dict['surface-node-id'][
            self.node_fields.dict['mechanical-node-boundary-label'] == self.geometry.valve_plug].astype(int)
        mechanical_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        mechanical_node_label_global[mechanical_lvnodes] = self.geometry.lv_endocardium
        mechanical_node_label_global[mechanical_rvnodes] = self.geometry.rv_endocardium
        mechanical_node_label_global[mechanical_epinodes] = self.geometry.epicardium
        mechanical_node_label_global[mechanical_valvenodes] = self.geometry.valve_plug
        self.node_fields.add_field(mechanical_node_label_global, 'mechanical-node-label-global', 'nodefield')
        self.node_fields.add_field(mechanical_lvnodes, 'mechanical-lvnodes', 'nodefield')
        self.node_fields.add_field(mechanical_rvnodes, 'mechanical-rvnodes', 'nodefield')
        self.node_fields.add_field(mechanical_epinodes, 'mechanical-epinodes', 'nodefield')
        self.node_fields.add_field(mechanical_valvenodes, 'mechanical-valvenodes', 'nodefield')

        # Get LV and RV endocardium faces for Eikonal solution.
        ep_lvfaces = self.geometry.triangles[self.element_fields.dict['ep-element-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        ep_rvfaces = self.geometry.triangles[self.element_fields.dict['ep-element-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        self.element_fields.add_field(ep_lvfaces, 'ep-lvfaces', 'elementfield')
        self.element_fields.add_field(ep_rvfaces, 'ep-rvfaces', 'elementfield')


    def write_to_csv(self):
        print('Writing .csv files to '+self.output_dir)
        self.geometry.save_to_csv(self.output_dir)
        self.node_fields.save_to_csv(self.output_dir)
        self.element_fields.save_to_csv(self.output_dir)

    def write_to_ensight(self):
        print('Writing Ensight files to ' + self.output_dir)
        self.geometry.save_to_ensight(self.output_dir)
        self.node_fields.save_to_ensight(self.output_dir, casename=self.output_name+'_nodefield', geometry=self.geometry)
        self.element_fields.save_to_ensight(self.output_dir, casename=self.output_name+'_elementfield', geometry=self.geometry)