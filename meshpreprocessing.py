import vtk
from vtk.util import numpy_support as VN
import pymp
from meshstructure import MeshStructure
from myformat import *
import copy
from xml.dom import minidom
# import pymesh
import json
from matplotlib import pyplot as plt
import pandas as pd

class MeshPreprocessing(MeshStructure):
    def __init__(self, vtk_name, name, input_dir, geometric_data_dir, max_cores_used, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, max_cores_used=max_cores_used, verbose=verbose)
        # Read and write geometry
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.vtk_name = vtk_name
        self.input_dir = input_dir
        # self.read_geometry_from_vtk_rodero()
        self.lv_endocardium = 3 # default values for Rodero meshes
        self.rv_endocardium = 2
        self.epicardium = 1
        self.valve_plug = 4
        # Generate fields
        # self.generate_boundary_data_rodero()
        # self.generate_fibre_sheet_normal()

        # self.save()
        # self.check_fields_for_qrs_inference()
        # self.check_fields_for_twave_inference()

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
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_boundarynodefield_ep-lvnodes.csv')
        assert os.path.exists(self.geometric_data_dir + self.geometry.name + '_boundarynodefield_ep-rvnodes.csv')
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

    def read_geometry_from_vtu_hannah_smith(self, vtu_name, save=False):
        print('Reading VTU files from ' + self.input_dir + vtu_name + '.vtu')
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(self.input_dir + vtu_name + '.vtu')
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
        unfolded_edges = np.concatenate((self.geometry.edges, np.flip(self.geometry.edges, axis=1))).astype(int)
        aux = [[] for i in range(0, self.geometry.number_of_nodes, 1)]
        for i in range(0, len(unfolded_edges)):
            aux[unfolded_edges[i, 0]].append(unfolded_edges[i, 1])
        neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
        self.neighbours = neighbours
        self.geometry.tetrahedron_centres = (nodes_xyz[self.geometry.tetrahedrons[:, 0], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 1], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 2], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 3], :]) / 4.
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('tv')), 'tv', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('tm')), 'tm', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('ab')), 'ab', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('rt')), 'rt', 'nodefield')
        def _normalise_vector(vector):
            m = np.linalg.norm(vector)
            if m > 0:
                return vector / m
            else:
                return vector

        print('Reading vtk file from: ' + self.input_dir + 'transmural_vectors.vtk')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(self.input_dir + 'transmural_vectors.vtk')
        reader.Update()
        data = reader.GetOutput()
        transmural_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(transmural_vector.shape[0]):
            transmural_vector[node_i, :] = _normalise_vector(transmural_vector[node_i, :])
        self.node_fields.add_field(transmural_vector, 'transmural-vector', 'nodefield')
        print('Reading vtk file from: ' + self.input_dir + 'longitudinal_vectors.vtk')
        reader.SetFileName(self.input_dir + 'longitudinal_vectors.vtk')
        reader.Update()
        data = reader.GetOutput()
        longitudinal_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(longitudinal_vector.shape[0]):
            longitudinal_vector[node_i, :] = _normalise_vector(longitudinal_vector[node_i, :])
        self.node_fields.add_field(longitudinal_vector, 'longitudinal-vector', 'nodefield')
        print('Reading vtk file from: ' + self.input_dir + 'circumferential_vectors.vtk')
        reader.SetFileName(self.input_dir + 'circumferential_vectors.vtk')
        reader.Update()
        data = reader.GetOutput()
        circumferential_vector = VN.vtk_to_numpy(data.GetPointData().GetArray('ScalarGradient'))
        for node_i in range(circumferential_vector.shape[0]):
            circumferential_vector[node_i, :] = _normalise_vector(circumferential_vector[node_i, :])
        self.node_fields.add_field(circumferential_vector, 'circumferential-vector', 'nodefield')
        if save:
            self.save()

    def add_cap_to_mesh(self):
        # This function only goes as far as generating new surface data points for the endo and epi surfaces
        # These points are then used to generate surface meshes using meshlab
        # Then they are passed to Ruben Doste's pipeline to generate UVCs and fibre orientations.
        vtk_names = ['epicardium.vtk', 'lv_endocardium.vtk', 'rv_endocardium.vtk']
        edge_names = ['lv_endo_edge_node_start_vectors.csv',
                  'rv_endo_edge_node_start_vectors.csv',
                  'epicardium_edge_node_start_vectors.csv']
        save_names = ['lvendo_cap_and_surface.txt',
                  'rvendo_cap_and_surface.txt',
                  'epicardium_cap_and_surface.txt']
        lid_rise_percentages = [0.15, 0.15, 0.2]

        all_start_vectors = []
        for edge_name in edge_names:
            data = pd.read_csv(self.geometric_data_dir + edge_name)
            edge_node_idx = data['PointIds'].values
            start_vectors = np.zeros((edge_node_idx.shape[0], 3))
            start_vectors[:, 0] = data['ScalarGradient:0'].values
            start_vectors[:, 1] = data['ScalarGradient:1'].values
            start_vectors[:, 2] = data['ScalarGradient:2'].values
            all_start_vectors.append(start_vectors)
        all_start_vectors = np.concatenate(all_start_vectors, axis=0)
        normal_vector = np.average(all_start_vectors, axis=0)
        for vtk_name, edge_name, save_name, lid_rise_percentage in zip(vtk_names, edge_names, save_names,
                                                                       lid_rise_percentages):
            print(' Reading ', self.geometric_data_dir + vtk_name, '...')
            reader = vtk.vtkUnstructuredGridReader()
            reader.SetFileName(self.geometric_data_dir + vtk_name)
            reader.ReadAllVectorsOn()
            reader.ReadAllScalarsOn()
            reader.Update()
            data = reader.GetOutput()
            nodes_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())  # Convert from mm to cm
            fig = plt.figure(figsize=(8, 8))
            ax = fig.add_subplot(projection='3d')
            ax.scatter(nodes_xyz[:, 0], nodes_xyz[:, 1], nodes_xyz[:, 2], marker='*', c='b')
            node_ids = VN.vtk_to_numpy(data.GetPointData().GetArray('PointIds'))
            existing_lv_length = np.amax(np.abs(np.sum(normal_vector * nodes_xyz, axis=1, keepdims=True)))
            lid_rise = existing_lv_length * lid_rise_percentage
            data = pd.read_csv(self.geometric_data_dir + edge_name)
            edge_node_idx = data['PointIds'].values
            start_vectors = np.zeros((edge_node_idx.shape[0], 3))
            start_vectors[:, 0] = data['ScalarGradient:0'].values
            start_vectors[:, 1] = data['ScalarGradient:1'].values
            start_vectors[:, 2] = data['ScalarGradient:2'].values
            # edge_node_idx = np.loadtxt('DTI003/epi_edge_node_idx.csv', delimiter=',')[:, 0].astype(int)
            edge_nodes = np.zeros((edge_node_idx.shape[0], 3))
            for i, idx in enumerate(edge_node_idx):
                print(nodes_xyz[np.where(node_ids == idx)[0], :])
                print(edge_nodes[i, :])
                edge_nodes[i, :] = nodes_xyz[np.where(node_ids == idx)[0], :]
            # ax.scatter(edge_nodes[:, 0], edge_nodes[:, 1], edge_nodes[:, 2], marker='*', c='r')
            # start_vectors = np.zeros((edge_nodes.shape[0], 3))
            # start_vectors[:, 0] = normal_vector[0]
            # start_vectors[:, 1] = normal_vector[1]
            # start_vectors[:, 2] = normal_vector[2]
            print('Generating cap data for ', edge_nodes.shape[0], ' edge nodes...')
            cap_data = generate_cap_surface_data(edge_nodes=edge_nodes, start_vectors=start_vectors,
                                                 longitudinal_vector=normal_vector, lid_rise=lid_rise,
                                                 flat_percentage=0.9, resolution_along_radius=0.05)
            cap_and_surface_data = np.concatenate((cap_data, nodes_xyz), axis=0)
            ax.scatter(cap_data[:, 0], cap_data[:, 1], cap_data[:, 2], marker='*', c='m')
            plt.show()
            np.savetxt(self.geometric_data_dir + save_name, cap_and_surface_data)
            reader = None

    def read_geometry_from_xml_vtk_cardiax_ellipsoid(self, xml_name, vtk_name, save=False):
        # Fenics XML element numbering is compatible with Alya. CardiaX number will give negative Jacobians.
        print('Reading xml file from: ' + self.input_dir + xml_name + '.xml')
        file = minidom.parse(self.input_dir + xml_name + '.xml')
        self.geometry.number_of_elements = int(file.getElementsByTagName('cells')[0].getAttribute('size'))
        tets = file.getElementsByTagName('cells')[0].getElementsByTagName('tetrahedron')
        self.geometry.tetrahedrons = np.zeros((self.geometry.number_of_elements, 4)).astype(int)
        for element_i in range(self.geometry.number_of_elements):
            self.geometry.tetrahedrons[element_i, 0] = int(tets[element_i].getAttribute('v0'))
            self.geometry.tetrahedrons[element_i, 1] = int(tets[element_i].getAttribute('v1'))
            self.geometry.tetrahedrons[element_i, 2] = int(tets[element_i].getAttribute('v2'))
            self.geometry.tetrahedrons[element_i, 3] = int(tets[element_i].getAttribute('v3'))
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
        unfolded_edges = np.concatenate((self.geometry.edges, np.flip(self.geometry.edges, axis=1))).astype(int)
        aux = [[] for i in range(0, self.geometry.number_of_nodes, 1)]
        for edge_i in range(0, len(unfolded_edges)):
            aux[unfolded_edges[edge_i, 0]].append(unfolded_edges[edge_i, 1])
        neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
        self.neighbours = neighbours
        # Read nodes/vertices
        vertices = file.getElementsByTagName('vertices')[0].getElementsByTagName('vertex')
        self.geometry.number_of_nodes = int(file.getElementsByTagName('vertices')[0].getAttribute('size'))
        nodes_xyz = np.zeros((self.geometry.number_of_nodes, 3))
        for node_i in range(self.geometry.number_of_nodes):
            nodes_xyz[node_i, 0] = float(vertices[node_i].getAttribute('x'))
            nodes_xyz[node_i, 1] = float(vertices[node_i].getAttribute('y'))
            nodes_xyz[node_i, 2] = float(vertices[node_i].getAttribute('z'))
        self.geometry.nodes_xyz = nodes_xyz
        self.geometry.tetrahedron_centres = (nodes_xyz[self.geometry.tetrahedrons[:, 0], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 1], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 2], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 3], :]) / 4.

        print('Reading vtk file from: ' + self.input_dir + vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.input_dir + vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()

        # Read faces from VTK file.
        self.geometry.number_of_triangles = data.GetNumberOfCells()
        self.geometry.triangles = np.zeros([self.geometry.number_of_triangles, 3])
        self.boundary_node_fields.add_field(data=VN.vtk_to_numpy(data.GetPointData().GetArray('PointIds')).astype(int),
                                            data_name='surface-node-id', field_type='boundarynodefield')
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(0)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(1)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)
        materials = np.zeros(self.geometry.number_of_elements).astype(int)
        materials[:] = 1
        self.materials.add_field(data=materials, data_name='tetra', field_type='material')
        if save:
            self.save()

    def read_geometry_from_vtk_cardiax_ellipsoid(self, save=False):
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
            self.geometry.tetrahedrons[i, 1  ] = int(data.GetCell(i).GetPointId(0))
            self.geometry.tetrahedrons[i, 2] = int(data.GetCell(i).GetPointId(1))
            self.geometry.tetrahedrons[i, 0] = int(data.GetCell(i).GetPointId(2))
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
        unfolded_edges = np.concatenate((self.geometry.edges, np.flip(self.geometry.edges, axis=1))).astype(int)
        aux = [[] for i in range(0, self.geometry.number_of_nodes, 1)]
        for i in range(0, len(unfolded_edges)):
            aux[unfolded_edges[i, 0]].append(unfolded_edges[i, 1])
        neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
        self.neighbours = neighbours
        self.geometry.tetrahedron_centres = (nodes_xyz[self.geometry.tetrahedrons[:, 0], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 1], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 2], :] +
                                             nodes_xyz[self.geometry.tetrahedrons[:, 3], :]) / 4.

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
        self.boundary_node_fields.add_field(data=VN.vtk_to_numpy(data.GetPointData().GetArray('PointIds')).astype(int),
                                            data_name='surface-node-id', field_type='boundarynodefield')
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(0)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(1)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)
        materials = np.zeros(self.geometry.number_of_elements).astype(int)
        materials[:] = 1
        self.materials.add_field(data=materials, data_name='tetra', field_type='material')
        if save:
            self.save()

    def generate_boundary_data_cardiax_ellipsoid(self, save=False):
        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        elem_ids = VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId'))
        self.boundary_element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId')).astype(int),
                                               'boundary-label', 'boundaryelementfield')
        self.boundary_node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RegionId')).astype(int),
                                            'boundary-label', 'boundarynodefield')
        for i in range(0, len(self.boundary_element_fields.dict['boundary-label'])):
            if self.boundary_element_fields.dict['boundary-label'][i] == 1:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.epicardium  # Epicardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == 2:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium  # LV endocardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == 0:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.lid  # Basal truncated surface

        for i in range(0, len(self.boundary_node_fields.dict['boundary-label'])):
            if self.boundary_node_fields.dict['boundary-label'][i] == 1:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.epicardium  # Epicardial
            elif self.boundary_node_fields.dict['boundary-label'][i] == 2:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium  # LV endocardial
            elif self.boundary_node_fields.dict['boundary-label'][i] == 0:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.lid  # Basal truncated surface

        # Save input global label
        input_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict[
                'boundary-label'] == self.geometry.lv_endocardium]] = self.geometry.lv_endocardium
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict[
                'boundary-label'] == self.geometry.epicardium]] = self.geometry.epicardium
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict[
                'boundary-label'] == self.geometry.lid]] = self.geometry.lid
        self.boundary_node_fields.add_field(input_node_label_global, 'input-boundary-label',
                                            'boundarynodefield')
        # Set up boundary label for mechanics
        mechanical_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        mechanical_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        for i in range(self.geometry.triangles.shape[0]):
            for j in range(self.geometry.triangles.shape[1]):
                local_index = \
                    np.nonzero(self.boundary_node_fields.dict['surface-node-id'] == self.geometry.triangles[i, j])[0][0]
                # Epicardial plug surface
                mechanical_element_boundary_label[i] = \
                    self.boundary_element_fields.dict['boundary-label'][i]
                mechanical_node_boundary_label[local_index] = self.boundary_node_fields.dict['boundary-label'][
                    local_index]
        self.boundary_node_fields.add_field(mechanical_node_boundary_label, 'mechanical-node-boundary-label',
                                            'boundarynodefield')
        self.boundary_element_fields.add_field(mechanical_element_boundary_label, 'mechanical-element-boundary-label',
                                               'boundaryelementfield')

        # Get LV and RV endocardium nodes
        mechanical_lvnodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.lv_endocardium].astype(
            int)
        mechanical_epinodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.epicardium].astype(int)
        mechanical_lidnodes = self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.lid].astype(int)
        mechanical_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        mechanical_node_label_global[mechanical_lvnodes] = self.geometry.lv_endocardium
        mechanical_node_label_global[mechanical_epinodes] = self.geometry.epicardium
        mechanical_node_label_global[mechanical_lidnodes] = self.geometry.lid
        self.boundary_node_fields.add_field(mechanical_node_label_global, 'mechanical-node-label-global',
                                            'boundarynodefield')
        self.boundary_node_fields.add_field(mechanical_lvnodes, 'mechanical-lvnodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(mechanical_epinodes, 'mechanical-epinodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(mechanical_lidnodes, 'mechanical-lidnodes', 'boundarynodefield')
        if save:
            self.save()

    def generate_boundary_data_UKB(self, save=False):
        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        self.boundary_element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId')).astype(int),
                                               'boundary-label', 'boundaryelementfield')
        self.boundary_node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RegionId')).astype(int),
                                            'boundary-label', 'boundarynodefield')
        # Read faces
        self.geometry.number_of_triangles = data.GetNumberOfCells()
        self.geometry.triangles = np.zeros([self.geometry.number_of_triangles, 3])
        self.boundary_node_fields.add_field(data=VN.vtk_to_numpy(data.GetPointData().GetArray('PointIds')).astype(int),
                                            data_name='surface-node-id', field_type='boundarynodefield')
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(0)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(1)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)
        print('Assigning boundary labels to each surface element')
        for i in range(0, len(self.boundary_element_fields.dict['boundary-label'])):
            if self.boundary_element_fields.dict['boundary-label'][i] == 0:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.epicardium  # Epicardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == 1:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium  # LV endocardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == 2:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.rv_endocardium  # RV endocardial
        print('Assigning boundary labels to each node')
        for i in range(0, len(self.boundary_node_fields.dict['boundary-label'])):
            if self.boundary_node_fields.dict['boundary-label'][i] == 0:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.epicardium  # Epicardial
            elif self.boundary_node_fields.dict['boundary-label'][i] == 1:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium  # LV endocardial
            elif self.boundary_node_fields.dict['boundary-label'][i] == 2:
                self.boundary_node_fields.dict['boundary-label'][i] = self.geometry.rv_endocardium  # RV endocardial
        # Save input global label
        print('Saving boundary labels as boundary node fields')
        input_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict[
                'boundary-label'] == self.geometry.lv_endocardium]] = self.geometry.lv_endocardium
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict[
                'boundary-label'] == self.geometry.rv_endocardium]] = self.geometry.rv_endocardium
        input_node_label_global[self.boundary_node_fields.dict['surface-node-id'][
            self.boundary_node_fields.dict['boundary-label'] == self.geometry.epicardium]] = self.geometry.epicardium
        self.boundary_node_fields.add_field(input_node_label_global, 'input-boundary-label', 'boundarynodefield')

        # Set up boundary label for EP and mechanics
        print('Set up boundary labels for mechanics...')
        # mechanical_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        # mechanical_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        # epi_nodes = np.where(self.boundary_node_fields.dict['boundary-label'] == self.geometry.epicardium)[0]
        # pericardial_nodes = epi_nodes[np.where(self.node_fields.dict['ab'][epi_nodes] > self.geometry.pericardial_ab_extent)[0]]
        # mechanical_node_boundary_label[:] = self.boundary_node_fields.dict['boundary-label'][:]
        # mechanical_node_boundary_label[pericardial_nodes] = self.valve_plug
        # self.boundary_node_fields.add_field(mechanical_node_boundary_label, 'mechanical-node-boundary-label',
        #                                     'boundarynodefield')
        #
        # for i in self.geometry.number_of_triangles:
        #

        mechanical_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        # ep_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        mechanical_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        # ep_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        for i in range(self.geometry.triangles.shape[0]):
            for j in range(self.geometry.triangles.shape[1]):
                local_index = \
                    np.nonzero(self.boundary_node_fields.dict['surface-node-id'] == self.geometry.triangles[i, j])[
                        0][0]
                # Epicardial plug surface
                if (self.boundary_node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & (
                        self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10):
                    mechanical_node_boundary_label[local_index] = self.valve_plug
                    # ep_node_boundary_label[local_index] = self.valve_plug
                    mechanical_element_boundary_label[i] = self.valve_plug
                else:
                    mechanical_element_boundary_label[i] = \
                        self.boundary_element_fields.dict['boundary-label'][i]
                    # ep_element_boundary_label[local_index] = self.boundary_element_fields.dict['boundary-label'][i]
                    mechanical_node_boundary_label[local_index] = self.boundary_node_fields.dict['boundary-label'][
                        local_index]
                    # ep_node_boundary_label[local_index] = self.boundary_node_fields.dict['boundary-label'][
                    #     local_index]
                # Restrict mechanical epicardium to only 80% of apex-to-base axis.
                if (self.boundary_node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & \
                        (self.node_fields.dict['ab'][
                             self.geometry.triangles[i, j]] > self.geometry.pericardial_ab_extent):
                    mechanical_element_boundary_label[
                        i] = self.valve_plug  # Apply epicardial spring BC only to 80% of apex-to-base
                    mechanical_node_boundary_label[local_index] = self.valve_plug
        self.boundary_node_fields.add_field(mechanical_node_boundary_label, 'mechanical-node-boundary-label',
                                            'boundarynodefield')
        # self.boundary_node_fields.add_field(ep_node_boundary_label, 'ep-node-boundary-label', 'boundarynodefield')
        self.boundary_element_fields.add_field(mechanical_element_boundary_label,
                                               'mechanical-element-boundary-label',
                                               'boundaryelementfield')
        # self.boundary_element_fields.add_field(ep_element_boundary_label, 'ep-element-boundary-label',
        #                                        'boundaryelementfield')

    # # Get LV and RV endocardium nodes
        # print('Get endocardial nodes for activation layer')
        # ep_lvnodes = self.boundary_node_fields.dict['surface-node-id'][
        #     self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        # ep_rvnodes = self.boundary_node_fields.dict['surface-node-id'][
        #     self.boundary_node_fields.dict['mechanical-node-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        # ep_node_label_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        # ep_node_label_global[ep_lvnodes] = self.geometry.lv_endocardium
        # ep_node_label_global[ep_rvnodes] = self.geometry.rv_endocardium
        # self.boundary_node_fields.add_field(ep_node_label_global, 'ep-node-label-global', 'boundarynodefield')
        # self.boundary_node_fields.add_field(ep_lvnodes, 'ep-lvnodes', 'boundarynodefield')
        # self.boundary_node_fields.add_field(ep_rvnodes, 'ep-rvnodes', 'boundarynodefield')

        print('Saving boundary node and element fields and writing to ensight at ', self.geometric_data_dir + 'ensight/')
        self.boundary_node_fields.save_to_csv(self.geometric_data_dir)
        self.boundary_element_fields.save_to_csv(self.geometric_data_dir)
        self.geometry.save_to_ensight(self.geometric_data_dir + 'ensight/')
        self.boundary_node_fields.save_to_ensight(self.geometric_data_dir + 'ensight/',
                                                  casename=self.name + '_boundary_node_fields', geometry=self.geometry)
        self.boundary_element_fields.save_to_ensight(self.geometric_data_dir + 'ensight/',
                                                     casename=self.name + '_boundary_element_fields',
                                                     geometry=self.geometry)
        if save:
            self.save()

    def read_geometry_from_vtk_rodero(self, save=False):
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
        unfolded_edges = np.concatenate((self.geometry.edges, np.flip(self.geometry.edges, axis=1))).astype(int)
        aux = [[] for i in range(0, self.geometry.number_of_nodes, 1)]
        for i in range(0, len(unfolded_edges)):
            aux[unfolded_edges[i, 0]].append(unfolded_edges[i, 1])
        neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
        self.neighbours = neighbours
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
        # print('Map UVC to element fields using nearest node')
        # map = map_indexes(points_to_map_xyz=self.geometry.tetrahedron_centres, reference_points_xyz=self.geometry.nodes_xyz)
        reader.SetFileName(self.input_dir + self.vtk_name + '_element_fields.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        self.element_fields.add_field(data=VN.vtk_to_numpy(data.GetCellData().GetArray('V.dat')), data_name='tv', field_type='elementfield')
        self.element_fields.add_field(data=VN.vtk_to_numpy(data.GetCellData().GetArray('RHO.dat')), data_name='tm', field_type='elementfield')
        self.element_fields.add_field(data=VN.vtk_to_numpy(data.GetCellData().GetArray('Z.dat')), data_name='ab', field_type='elementfield')
        self.element_fields.add_field(data=VN.vtk_to_numpy(data.GetCellData().GetArray('PHI.dat')), data_name='rt', field_type='elementfield')
        # print('Getting ab, rt, tm, and tv at tetrahedron centres')
        # ab_tetra_centres = self.node_fields.dict['ab'][map]
        # rt_tetra_centres = self.node_fields.dict['rt'][map]
        # tm_tetra_centres = self.node_fields.dict['tm'][map]
        # tv_tetra_centres = self.node_fields.dict['tv'][map]
        # self.element_fields.add_field(data=ab_tetra_centres, data_name='ab', field_type='elementfield')
        # self.element_fields.add_field(data=rt_tetra_centres, data_name='rt', field_type='elementfield')
        # self.element_fields.add_field(data=tm_tetra_centres, data_name='tm', field_type='elementfield')
        # self.element_fields.add_field(data=tv_tetra_centres, data_name='tv', field_type='elementfield')

        # print('Convert from UVC to Cobiveco')
        # self.convert_uvc_to_cobiveco()
        # print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_fibres_at_nodes.vtk')
        # reader.SetFileName(self.input_dir + self.vtk_name + '_fibres_at_nodes.vtk')
        # reader.Update()
        # data = reader.GetOutput()
        # self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('fibres')), 'fibres', 'nodefield')
        # self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('sheets')), 'sheets', 'nodefield')

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
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()

        # Read faces
        self.geometry.number_of_triangles = data.GetNumberOfCells()
        self.geometry.triangles = np.zeros([self.geometry.number_of_triangles, 3])
        self.boundary_node_fields.add_field(data=VN.vtk_to_numpy(data.GetPointData().GetArray('Ids')).astype(int),
                                            data_name='surface-node-id', field_type='boundarynodefield')
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(0)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(1)]),
                int(self.boundary_node_fields.dict['surface-node-id'][data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)
        materials = np.zeros(self.element_fields.dict['tv-element'].shape[0]).astype(int)
        for element_i in range(self.element_fields.dict['tv-element'].shape[0]):
            if self.element_fields.dict['tv-element'][element_i] > 3:
                materials[element_i] = 2
            else:
                materials[element_i] = 1
        self.materials.add_field(data=materials, data_name='tetra', field_type='material')
        if save:
            self.save()

    def convert_uvc_to_cobiveco(self):
        # Map input geometry Cobiveco coordinates to output-style UVC
        cobiveco = {}
        cobiveco['ab'] = copy.deepcopy(self.node_fields.dict['ab'])
        cobiveco['rt'] = copy.deepcopy(self.node_fields.dict['rt'])
        # Handle mapping of TV split separately:
        cobiveco['ab'] = map_ventricular_coordinates(input_ranges=[[0, 1]],
                                                     output_ranges=[[0, 1.5]],
                                                     input_coordinate=self.node_fields.dict['ab'])
        cobiveco['tv'] = map_ventricular_coordinates(input_ranges=[[self.geometry.lv, self.geometry.rv]], output_ranges=[[0, 1]],
                                                     input_coordinate=self.node_fields.dict['tv'])
        cobiveco['tm'] = copy.deepcopy(self.node_fields.dict['tm'])
        rvlv = copy.deepcopy(self.node_fields.dict['rt'])

        # Identify insertion points
        prev_ab = 0.0
        has_done_last_iteration = False
        for current_ab in np.arange(0, 1.01, 0.01):
        # for ab_section in [0.6]:
            mid_ventricular_nodes = np.nonzero((self.node_fields.dict['ab'] >= prev_ab) & (self.node_fields.dict['ab'] <= current_ab))[0].astype(int)
            rv_mid_ventricular_nodes_meta_idx = np.nonzero(self.node_fields.dict['tv'][mid_ventricular_nodes] == self.geometry.rv)[0].astype(int)
            rv_mid_ventricular_nodes = mid_ventricular_nodes[rv_mid_ventricular_nodes_meta_idx].astype(int)
            candidate_insertion_nodes = []
            for node_i in range(rv_mid_ventricular_nodes.shape[0]):
                neighbours = self.neighbours[rv_mid_ventricular_nodes[node_i]]
                meta_idx = np.nonzero(self.node_fields.dict['tv'][neighbours] == self.geometry.lv)[0]
                if len(meta_idx) > 0:
                    candidate_insertion_nodes.append(rv_mid_ventricular_nodes[node_i])
            candidate_insertion_nodes = np.unique(candidate_insertion_nodes).astype(int)
            if (len(candidate_insertion_nodes) == 0):
                continue
            else:
                negative_candidate_nodes = candidate_insertion_nodes[np.nonzero(self.node_fields.dict['rt'][candidate_insertion_nodes] < 0)[0]]
                positive_candidate_nodes = candidate_insertion_nodes[np.nonzero(self.node_fields.dict['rt'][candidate_insertion_nodes] > 0)[0]]
                if (len(positive_candidate_nodes) == 0) or (len(negative_candidate_nodes) == 0):
                    continue
                else:
                    sorted_negative_rt = np.sort(self.node_fields.dict['rt'][negative_candidate_nodes])
                    negative_median_rt = self.node_fields.dict['rt'][negative_candidate_nodes[int(len(sorted_negative_rt)/2)]]
                    negative_median_meta_idx = np.nonzero(self.node_fields.dict['rt'][negative_candidate_nodes] == negative_median_rt)[0][0]
                    posterior_insertion_node = negative_candidate_nodes[negative_median_meta_idx]
                    sorted_positive_rt = np.sort(self.node_fields.dict['rt'][positive_candidate_nodes])
                    positive_median_rt = self.node_fields.dict['rt'][positive_candidate_nodes[int(len(sorted_positive_rt) / 2)]]
                    positive_median_meta_idx = np.nonzero(self.node_fields.dict['rt'][positive_candidate_nodes] == positive_median_rt)[0][0]
                    anterior_insertion_node = positive_candidate_nodes[positive_median_meta_idx]
                    rt_anterior_insertion = self.node_fields.dict['rt'][anterior_insertion_node]
                    rt_posterior_insertion = self.node_fields.dict['rt'][posterior_insertion_node]
                    cobiveco_rv_septum_indices = []
                    for node_i in range(mid_ventricular_nodes.shape[0]):
                        if (self.node_fields.dict['rt'][mid_ventricular_nodes[node_i]] > rt_posterior_insertion) & \
                                (self.node_fields.dict['rt'][mid_ventricular_nodes[node_i]] < rt_anterior_insertion) & \
                                (self.node_fields.dict['tv'][mid_ventricular_nodes[node_i]] == self.geometry.lv) & \
                                (self.node_fields.dict['tm'][mid_ventricular_nodes[node_i]] > 0.6):
                            cobiveco['tv'][mid_ventricular_nodes[node_i]] = 1
                            cobiveco_rv_septum_indices.append(mid_ventricular_nodes[node_i])
                    rt_mid_septum = 0.0
                    rt_rv_lateral = 0.0
                    # rt_lv_lateral = 0.0 or 1.0
                    cobiveco_lv_indices = mid_ventricular_nodes[np.where(cobiveco['tv'][mid_ventricular_nodes] == 0)]
                    cobiveco_rv_indices = mid_ventricular_nodes[np.where(cobiveco['tv'][mid_ventricular_nodes] == 1)]
                    # Translation of RT is ventricle dependent
                    cobiveco_rt_lv_lateral = 0.4
                    cobiveco_rt_rv_lateral = 0.4
                    cobiveco_rt_mid_septum = 0.8
                    # cobiveco_rt_posterior_insertion = 0.0 or 1.0
                    cobiveco_rt_anterior_insertion = (1 - (cobiveco_rt_mid_septum - cobiveco_rt_lv_lateral)) * rt_anterior_insertion/np.pi
                    cobiveco['rt'][cobiveco_lv_indices] = map_ventricular_coordinates(
                        input_ranges=[[-np.pi, rt_posterior_insertion],
                                      [rt_posterior_insertion, rt_mid_septum],
                                      [rt_mid_septum, rt_anterior_insertion],
                                      [rt_anterior_insertion, np.pi]],
                        output_ranges=[[cobiveco_rt_lv_lateral, 0],
                                       [1, cobiveco_rt_mid_septum],
                                       [cobiveco_rt_mid_septum, cobiveco_rt_anterior_insertion],
                                       [cobiveco_rt_anterior_insertion, cobiveco_rt_lv_lateral]],
                        input_coordinate=self.node_fields.dict['rt'][cobiveco_lv_indices])
                    cobiveco['rt'][cobiveco_rv_indices] = map_ventricular_coordinates(
                        input_ranges=[[rt_rv_lateral, rt_anterior_insertion],
                                      [rt_posterior_insertion, rt_rv_lateral]],
                        output_ranges=[[cobiveco_rt_rv_lateral, cobiveco_rt_anterior_insertion],
                                       [0, cobiveco_rt_rv_lateral]],
                        input_coordinate=self.node_fields.dict['rt'][cobiveco_rv_indices])  # RT
                    cobiveco['rt'][cobiveco_rv_septum_indices] = map_ventricular_coordinates(
                        input_ranges=[[rt_posterior_insertion, rt_anterior_insertion]],
                        output_ranges=[[1, cobiveco_rt_anterior_insertion]],
                        input_coordinate=self.node_fields.dict['rt'][cobiveco_rv_septum_indices])
                    cobiveco['tm'][mid_ventricular_nodes] = map_ventricular_coordinates(input_ranges=[[0, 1]],
                                                                 output_ranges=[[1, 0]],
                                                                 input_coordinate=self.node_fields.dict['tm'][mid_ventricular_nodes])  # TM
                    rvlv_lv_lateral = 1.0
                    rvlv_rv_lateral = 0.0
                    rvlv_posterior_insertion = 0.0
                    rvlv_anterior_insertion = 0.0
                    rvlv[mid_ventricular_nodes] = map_ventricular_coordinates(
                        input_ranges=[[-np.pi, rt_posterior_insertion],
                                      [rt_anterior_insertion, np.pi]],
                        output_ranges=[[rvlv_lv_lateral, rvlv_posterior_insertion],
                                       [rvlv_anterior_insertion, rvlv_lv_lateral]],
                        input_coordinate=self.node_fields.dict['rt'][mid_ventricular_nodes])
                    rvlv[cobiveco_rv_indices] = rvlv_rv_lateral
                    prev_ab = current_ab
                    if current_ab == 1.01:
                        has_done_last_iteration = True
        if not has_done_last_iteration:
            current_ab = 1.01
            for prev_ab in np.arange(prev_ab, 0.0, -0.01):
                # for ab_section in [0.6]:
                mid_ventricular_nodes = \
                np.nonzero((self.node_fields.dict['ab'] >= prev_ab) & (self.node_fields.dict['ab'] <= current_ab))[
                    0].astype(int)
                rv_mid_ventricular_nodes_meta_idx = \
                np.nonzero(self.node_fields.dict['tv'][mid_ventricular_nodes] == self.geometry.rv)[0].astype(int)
                rv_mid_ventricular_nodes = mid_ventricular_nodes[rv_mid_ventricular_nodes_meta_idx].astype(int)
                candidate_insertion_nodes = []
                for node_i in range(rv_mid_ventricular_nodes.shape[0]):
                    neighbours = self.neighbours[rv_mid_ventricular_nodes[node_i]]
                    meta_idx = np.nonzero(self.node_fields.dict['tv'][neighbours] == self.geometry.lv)[0]
                    if len(meta_idx) > 0:
                        candidate_insertion_nodes.append(rv_mid_ventricular_nodes[node_i])
                candidate_insertion_nodes = np.unique(candidate_insertion_nodes).astype(int)
                if (len(candidate_insertion_nodes) == 0):
                    continue
                else:
                    negative_candidate_nodes = candidate_insertion_nodes[
                        np.nonzero(self.node_fields.dict['rt'][candidate_insertion_nodes] < 0)[0]]
                    positive_candidate_nodes = candidate_insertion_nodes[
                        np.nonzero(self.node_fields.dict['rt'][candidate_insertion_nodes] > 0)[0]]
                    if (len(positive_candidate_nodes) == 0) or (len(negative_candidate_nodes) == 0):
                        continue
                    else:
                        sorted_negative_rt = np.sort(self.node_fields.dict['rt'][negative_candidate_nodes])
                        negative_median_rt = self.node_fields.dict['rt'][
                            negative_candidate_nodes[int(len(sorted_negative_rt) / 2)]]
                        negative_median_meta_idx = \
                        np.nonzero(self.node_fields.dict['rt'][negative_candidate_nodes] == negative_median_rt)[0][0]
                        posterior_insertion_node = negative_candidate_nodes[negative_median_meta_idx]
                        sorted_positive_rt = np.sort(self.node_fields.dict['rt'][positive_candidate_nodes])
                        positive_median_rt = self.node_fields.dict['rt'][
                            positive_candidate_nodes[int(len(sorted_positive_rt) / 2)]]
                        positive_median_meta_idx = \
                        np.nonzero(self.node_fields.dict['rt'][positive_candidate_nodes] == positive_median_rt)[0][0]
                        anterior_insertion_node = positive_candidate_nodes[positive_median_meta_idx]
                        rt_anterior_insertion = self.node_fields.dict['rt'][anterior_insertion_node]
                        rt_posterior_insertion = self.node_fields.dict['rt'][posterior_insertion_node]
                        cobiveco_rv_septum_indices = []
                        for node_i in range(mid_ventricular_nodes.shape[0]):
                            if (self.node_fields.dict['rt'][mid_ventricular_nodes[node_i]] > rt_posterior_insertion) & \
                                    (self.node_fields.dict['rt'][
                                         mid_ventricular_nodes[node_i]] < rt_anterior_insertion) & \
                                    (self.node_fields.dict['tv'][mid_ventricular_nodes[node_i]] == self.geometry.lv) & \
                                    (self.node_fields.dict['tm'][mid_ventricular_nodes[node_i]] > 0.6):
                                cobiveco['tv'][mid_ventricular_nodes[node_i]] = 1
                                cobiveco_rv_septum_indices.append(mid_ventricular_nodes[node_i])
                        rt_mid_septum = 0.0
                        rt_rv_lateral = 0.0
                        # rt_lv_lateral = 0.0 or 1.0
                        cobiveco_lv_indices = mid_ventricular_nodes[
                            np.where(cobiveco['tv'][mid_ventricular_nodes] == 0)]
                        cobiveco_rv_indices = mid_ventricular_nodes[
                            np.where(cobiveco['tv'][mid_ventricular_nodes] == 1)]
                        # Translation of RT is ventricle dependent
                        cobiveco_rt_lv_lateral = 0.4
                        cobiveco_rt_rv_lateral = 0.4
                        cobiveco_rt_mid_septum = 0.8
                        # cobiveco_rt_posterior_insertion = 0.0 or 1.0
                        cobiveco_rt_anterior_insertion = (1 - (
                                    cobiveco_rt_mid_septum - cobiveco_rt_lv_lateral)) * rt_anterior_insertion / np.pi
                        cobiveco['rt'][cobiveco_lv_indices] = map_ventricular_coordinates(
                            input_ranges=[[-np.pi, rt_posterior_insertion],
                                          [rt_posterior_insertion, rt_mid_septum],
                                          [rt_mid_septum, rt_anterior_insertion],
                                          [rt_anterior_insertion, np.pi]],
                            output_ranges=[[cobiveco_rt_lv_lateral, 0],
                                           [1, cobiveco_rt_mid_septum],
                                           [cobiveco_rt_mid_septum, cobiveco_rt_anterior_insertion],
                                           [cobiveco_rt_anterior_insertion, cobiveco_rt_lv_lateral]],
                            input_coordinate=self.node_fields.dict['rt'][cobiveco_lv_indices])
                        cobiveco['rt'][cobiveco_rv_indices] = map_ventricular_coordinates(
                            input_ranges=[[rt_rv_lateral, rt_anterior_insertion],
                                          [rt_posterior_insertion, rt_rv_lateral]],
                            output_ranges=[[cobiveco_rt_rv_lateral, cobiveco_rt_anterior_insertion],
                                           [0, cobiveco_rt_rv_lateral]],
                            input_coordinate=self.node_fields.dict['rt'][cobiveco_rv_indices])  # RT
                        cobiveco['rt'][cobiveco_rv_septum_indices] = map_ventricular_coordinates(
                            input_ranges=[[rt_posterior_insertion, rt_anterior_insertion]],
                            output_ranges=[[1, cobiveco_rt_anterior_insertion]],
                            input_coordinate=self.node_fields.dict['rt'][cobiveco_rv_septum_indices])
                        cobiveco['tm'][mid_ventricular_nodes] = map_ventricular_coordinates(input_ranges=[[0, 1]],
                                                                                            output_ranges=[[1, 0]],
                                                                                            input_coordinate=
                                                                                            self.node_fields.dict['tm'][
                                                                                                mid_ventricular_nodes])  # TM
                        rvlv_lv_lateral = 1.0
                        rvlv_rv_lateral = 0.0
                        rvlv_posterior_insertion = 0.0
                        rvlv_anterior_insertion = 0.0
                        rvlv[mid_ventricular_nodes] = map_ventricular_coordinates(
                            input_ranges=[[-np.pi, rt_posterior_insertion],
                                          [rt_anterior_insertion, np.pi]],
                            output_ranges=[[rvlv_lv_lateral, rvlv_posterior_insertion],
                                           [rvlv_anterior_insertion, rvlv_lv_lateral]],
                            input_coordinate=self.node_fields.dict['rt'][mid_ventricular_nodes])
                        rvlv[cobiveco_rv_indices] = rvlv_rv_lateral
                        prev_ab = current_ab
                        break
        rt_posterior = -np.pi/2
        aprt_posterior = 0
        aprt_lateral = 0.5
        rt_anterior = np.pi/2
        aprt_anterior = 1
        rt_septum = 0.0
        aprt_septum = 0.5
        aprt = map_ventricular_coordinates(input_ranges=[[-np.pi, rt_posterior],
                                                         [rt_posterior, rt_septum],
                                                         [rt_septum, rt_anterior],
                                                         [rt_anterior, np.pi]],
                                           output_ranges=[[aprt_lateral, aprt_posterior],
                                                          [aprt_posterior, aprt_septum],
                                                          [aprt_septum, aprt_anterior],
                                                          [aprt_anterior, aprt_lateral]],
                                           input_coordinate=self.node_fields.dict['rt'])
        self.node_fields.add_field(data=aprt, data_name='aprt', field_type='nodefield')
        self.node_fields.add_field(data=aprt, data_name='cobiveco-aprt', field_type='nodefield')
        self.node_fields.add_field(data=rvlv, data_name='rvlv', field_type='nodefield')
        self.node_fields.add_field(data=rvlv, data_name='cobiveco-rvlv', field_type='nodefield')
        for key in cobiveco.keys():
            self.node_fields.add_field(data=cobiveco[key], data_name='cobiveco-'+key, field_type='nodefield')

    def generate_boundary_data_rodero(self, surface_label_json, save=False):
        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        elem_ids = VN.vtk_to_numpy(data.GetCellData().GetArray('Ids'))
        self.boundary_element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId')).astype(int),
                                               'boundary-label', 'boundaryelementfield')
        self.boundary_node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RegionId')).astype(int),
                                            'boundary-label', 'boundarynodefield')
        surface_labels = json.load(open(surface_label_json, 'r'))

        for i in range(0, len(self.boundary_element_fields.dict['boundary-label'])):
            if self.boundary_element_fields.dict['boundary-label'][i] == surface_labels['epi']:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.epicardium  # Epicardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == surface_labels['lv']:
                self.boundary_element_fields.dict['boundary-label'][i] = self.geometry.lv_endocardium  # LV endocardial
            elif self.boundary_element_fields.dict['boundary-label'][i] == surface_labels['rv']:
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
        self.boundary_node_fields.add_field(input_node_label_global, 'input-boundary-label', 'boundarynodefield')

        # Set up boundary label for EP and mechanics
        mechanical_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        ep_element_boundary_label = np.zeros(self.boundary_element_fields.dict['boundary-label'].shape)
        mechanical_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        ep_node_boundary_label = np.zeros(self.boundary_node_fields.dict['boundary-label'].shape)
        for i in range(self.geometry.triangles.shape[0]):
            for j in range(self.geometry.triangles.shape[1]):
                local_index = \
                np.nonzero(self.boundary_node_fields.dict['surface-node-id'] == self.geometry.triangles[i, j])[0][0]
                # Epicardial plug surface
                if (self.boundary_node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & (
                        self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10):
                    mechanical_node_boundary_label[local_index] = self.valve_plug
                    ep_node_boundary_label[local_index] = self.valve_plug
                    mechanical_element_boundary_label[i] = self.valve_plug
                else:
                    mechanical_element_boundary_label[i] = \
                    self.boundary_element_fields.dict['boundary-label'][i]
                    ep_element_boundary_label[local_index] = self.boundary_element_fields.dict['boundary-label'][i]
                    mechanical_node_boundary_label[local_index] = self.boundary_node_fields.dict['boundary-label'][
                        local_index]
                    ep_node_boundary_label[local_index] = self.boundary_node_fields.dict['boundary-label'][local_index]
                # Endocardial plug surface is isolated for EP simulations
                if self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10:
                    ep_node_boundary_label[local_index] = self.valve_plug
                # Restrict mechanical epicardium to only 80% of apex-to-base axis.
                if (self.boundary_node_fields.dict['boundary-label'][local_index] == self.geometry.epicardium) & \
                        (self.node_fields.dict['ab'][self.geometry.triangles[i, j]] > self.geometry.pericardial_ab_extent):
                    mechanical_element_boundary_label[i] = self.valve_plug  # Apply epicardial spring BC only to 80% of apex-to-base
                    mechanical_node_boundary_label[local_index] = self.valve_plug
        self.boundary_node_fields.add_field(mechanical_node_boundary_label, 'mechanical-node-boundary-label',
                                            'boundarynodefield')
        self.boundary_node_fields.add_field(ep_node_boundary_label, 'ep-node-boundary-label', 'boundarynodefield')
        self.boundary_element_fields.add_field(mechanical_element_boundary_label, 'mechanical-element-boundary-label',
                                               'boundaryelementfield')
        self.boundary_element_fields.add_field(ep_element_boundary_label, 'ep-element-boundary-label', 'boundaryelementfield')

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
        self.boundary_node_fields.add_field(ep_node_label_global, 'ep-node-label-global', 'boundarynodefield')
        self.boundary_node_fields.add_field(ep_lvnodes, 'ep-lvnodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(ep_rvnodes, 'ep-rvnodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(ep_epinodes, 'ep-epinodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(ep_valvenodes, 'ep-valvenodes', 'boundarynodefield')

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
        self.boundary_node_fields.add_field(mechanical_node_label_global, 'mechanical-node-label-global', 'boundarynodefield')
        self.boundary_node_fields.add_field(mechanical_lvnodes, 'mechanical-lvnodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(mechanical_rvnodes, 'mechanical-rvnodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(mechanical_epinodes, 'mechanical-epinodes', 'boundarynodefield')
        self.boundary_node_fields.add_field(mechanical_valvenodes, 'mechanical-valvenodes', 'boundarynodefield')

        # Get LV and RV endocardium faces for Eikonal solution.
        ep_lvfaces = self.geometry.triangles[
            self.boundary_element_fields.dict['ep-element-boundary-label'] == self.geometry.lv_endocardium].astype(int)
        ep_rvfaces = self.geometry.triangles[
            self.boundary_element_fields.dict['ep-element-boundary-label'] == self.geometry.rv_endocardium].astype(int)
        self.boundary_element_fields.add_field(ep_lvfaces, 'ep-lvfaces', 'boundaryelementfield')
        self.boundary_element_fields.add_field(ep_rvfaces, 'ep-rvfaces', 'boundaryelementfield')

        self.boundary_node_fields.save_to_csv(self.geometric_data_dir)
        self.boundary_element_fields.save_to_csv(self.geometric_data_dir)
        self.geometry.save_to_ensight(self.geometric_data_dir + 'ensight/')
        self.boundary_node_fields.save_to_ensight(self.geometric_data_dir + 'ensight/',
                                                  casename=self.name + '_boundary_node_fields', geometry=self.geometry)
        self.boundary_element_fields.save_to_ensight(self.geometric_data_dir + 'ensight/',
                                                     casename=self.name + '_boundary_element_fields',
                                                     geometry=self.geometry)
        if save:
            self.save()

    def generate_fibre_sheet_normal(self):
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
        self.plug_fibres(neighbours=self.neighbours, ortho=ortho)
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

    def generate_lower_resolution_mesh(self):
        # To run pymesh on ARCHER2, you need to install it first:
        # pip install pymesh2
        # this has a dependency on GMP, which you'll need to install like this (https://github.com/PyMesh/PyMesh/issues/293):
        # PREFIX=/my/target/install/location  # Change this to install GMP wherever you want.
        #
        # curl -fLO https://gmplib.org/download/gmp/gmp-6.2.1.tar.xz
        # tar -xf gmp-6.2.1.tar.xz
        # cd gmp-6.2.1
        # ./configure --prefix=$PREFIX
        # make all
        # make install
        # cp gmpxx.h $PREFIX/include
        # export GMP_INC=$PREFIX/include
        # export GMP_LIB=$PREFIX/lib
        mesh = pymesh.form_mesh(self.geometry.nodes_xyz, self.geometry.triangles, self.geometry.tetrahedrons)

        new_surface, info = pymesh.collapse_short_edges(mesh, 0.2)
        new_volume = pymesh.tetrahedralize(new_surface, 0.2, engine='cgal')
        print(self.geometric_data_dir + self.name + '_coarse')
        quit()
        if not os.path.exists(self.geometric_data_dir + self.name + '_coarse'):
            os.mkdir(os.path.join(self.geometric_data_dir, self.name + '_coarse'))
        coarse_geometry = Geometry(name=self.name+'_coarse', max_cores_used=self.max_cores_used, verbose=self.verbose)
        coarse_geometry.nodes_xyz = new_volume.vertices
        coarse_geometry.tetrahedrons = new_volume.voxels
        coarse_geometry.triangles = new_volume.faces
        coarse_geometry.save_to_ensight(output_dir=os.path.join(self.geometric_data_dir, self.name + '_coarse'))

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

def map_indexes(points_to_map_xyz, reference_points_xyz):
    mapped_indexes = pymp.shared.array((points_to_map_xyz.shape[0]), dtype=int)
    threadsNum = np.amin((int(multiprocessing.cpu_count()), 30))
    print('Mapping indices using ', str(threadsNum) , ' threads...')
    with pymp.Parallel(min(threadsNum, points_to_map_xyz.shape[0])) as p1:
        for conf_i in p1.range(points_to_map_xyz.shape[0]):
            mapped_indexes[conf_i] = np.argmin(
                np.linalg.norm(reference_points_xyz - points_to_map_xyz[conf_i, :], ord=2, axis=1)).astype(int)
    return mapped_indexes

def calculate_tetrahedral_mesh_volume(nodes, elements):
    D = nodes[elements[:, 3] , :]
    AD = nodes[elements[:, 0], :] - D
    BD = nodes[elements[:, 1], :] - D
    CD = nodes[elements[:, 2], :] - D
    D = None # Clear memory
    tV = np.reshape(np.abs(np.matmul(np.moveaxis(AD[:, :, np.newaxis], 1, -1), (np.cross(BD, CD)[:, :, np.newaxis]))),
                    elements.shape[0])
    AD = None
    BD = None
    CD = None
    return np.sum(tV) / 6.

def generate_cap_surface_data(edge_nodes, start_vectors, longitudinal_vector, lid_rise, flat_percentage, resolution_along_radius):
    # Create curves from each edge node to the middle
    middle_coord = np.mean(edge_nodes, axis=0)
    middle_lid_coord = middle_coord + longitudinal_vector*lid_rise
    # plot_3D_vector(ax, middle_lid_coord, longitudinal_vector, 'k')
    data = np.array([[0, 0, 0]])
    for edge_node, start_vector in zip(edge_nodes, start_vectors):
        # Create local axes
        x_vector = middle_coord - edge_node
        y_axis = np.cross(x_vector/np.linalg.norm(x_vector), longitudinal_vector)
        y_axis = y_axis/np.linalg.norm(y_axis)
        x_axis = np.cross(longitudinal_vector, y_axis)
        x_axis = x_axis/np.linalg.norm(x_axis)
        dx = resolution_along_radius * (x_vector)/np.linalg.norm(x_vector)
        # plot_3D_vector(ax, edge_node, dx, 'm')
        # plot_3D_vector(ax, edge_node, x_axis, 'g')
        # plot_3D_vector(ax, edge_node, y_axis, 'r')
        # plot_3D_vector(ax, edge_node, longitudinal_vector, 'b')
        n_points = np.floor(np.linalg.norm(x_vector)/ resolution_along_radius).astype(int)
        curve = np.zeros((n_points, 3))
        curve[0, :] = edge_node
        flat_x = flat_percentage * np.linalg.norm(x_vector)
        end_vector = x_vector / np.linalg.norm(x_vector)
        start_angle = np.arccos(np.dot(start_vector, end_vector))
        # plot_3D_vector(ax, edge_node, start_vector, 'k')
        current_coord = edge_node
        for i in range(1, n_points):
            angle = vector_angle(i * resolution_along_radius, start_angle, flat_x)
            basis_transformation = np.vstack((x_axis, y_axis, longitudinal_vector)).transpose()
            R = np.array([[np.cos(angle), 0, -np.sin(angle)], [0, 1, 0], [np.sin(angle), 0, np.cos(angle)]])
            new_vector = np.matmul(R, np.matmul(np.linalg.inv(basis_transformation),end_vector))
            current_coord = current_coord + dx
            new_vector_global = np.matmul(basis_transformation, new_vector)
            new_vector_global = new_vector_global/np.linalg.norm(new_vector_global)
            # plot_3D_vector(ax, current_coord, new_vector_global, 'b')
            # plot_3D_vector(ax, current_coord, dx, 'm')
            # scaling = resolution_along_radius / (np.dot(new_vector_global, dx/np.linalg.norm(dx)))
            curve[i, :] = curve[i-1,:] + new_vector_global * resolution_along_radius
            # ax.scatter(curve[0, 0], curve[0, 1], curve[0, 2], marker='o', c='m')
            # ax.scatter(curve[i, 0], curve[i, 1], curve[i, 2], marker='o', c='g')
            # plt.show()
            # quit()
        # lid_scale = lid_rise/np.dot(longitudinal_vector, curve[-1,:] - middle_coord)
        # x_scale = np.linalg.norm(x_vector)/np.dot(x_axis, np.linalg.norm(curve[-1,:] - edge_node))
        lid_scale = lid_rise / (np.dot(longitudinal_vector, curve[-1,:] - edge_node))
        x_scale =  np.linalg.norm(x_vector)/np.dot(x_axis, curve[-1,:] - edge_node)
        for i in range(0, n_points):
            curve_x = np.dot((curve[i,:]-edge_node), x_axis)
            curve_y = np.dot((curve[i,:]-edge_node), y_axis)
            curve_z = np.dot((curve[i,:]-edge_node), longitudinal_vector)
            curve[i, :] = (curve_x* x_axis*x_scale + curve_y*y_axis + curve_z * longitudinal_vector * lid_scale) + edge_node
        curve = curve[1:, :]
        # ax.scatter(curve[:, 0], curve[:, 1], curve[:, 2], marker='*', c='m')
        # plt.show()
        # quit()
        data = np.concatenate((data,curve))
    data = data[1:,:]
    return data


def vector_angle(x, start_angle, flat_x):
    if x <= flat_x:
        return start_angle/(flat_x)**2 * (x-flat_x)**2 # This is a quadratic function that begins at start angle and goes to zero at end_x
    else:
        return 0

def plot_3D_vector(ax, coord, vector, color):
    ax.quiver(coord[0], coord[1], coord[2], vector[0], vector[1], vector[2], color=color)

