import vtk
import numpy as np
from myformat import *
from vtk.util import numpy_support as VN


class VTKtoCSVEnsight:
    def __init__(self, vtk_name, output_name, input_dir, output_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
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
            self.geometry.edges.append([self.geometry.tetrahedrons[element_i,0], self.geometry.tetrahedrons[element_i,1]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 1], self.geometry.tetrahedrons[element_i, 2]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 2], self.geometry.tetrahedrons[element_i, 3]])
            self.geometry.edges.append(
                [self.geometry.tetrahedrons[element_i, 3], self.geometry.tetrahedrons[element_i, 0]])
        self.geometry.edges = np.unique(np.sort(self.geometry.edges,axis=1),axis=0)
        self.geometry.tetrahedron_centres = (nodes_xyz[self.geometry.tetrahedrons[:,0],:] +
                                             nodes_xyz[self.geometry.tetrahedrons[:,1],:] +
                                             nodes_xyz[self.geometry.tetrahedrons[:,2],:] +
                                             nodes_xyz[self.geometry.tetrahedrons[:,3],:])/4.
        # self.element_fields.lvrv = VN.vtk_to_numpy(data.GetCellData().GetArray('ID'))
        self.element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('ID')), 'tv_element', 'elementfield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('fibres')), 'fibres', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('sheets')), 'sheets', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('V.dat')), 'tv', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RHO.dat')), 'tm', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('Z.dat')), 'ab', 'nodefield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('PHI.dat')), 'rt', 'nodefield')

        # self.node_fields.fibres = VN.vtk_to_numpy(data.GetCellData().GetArray('fibres'))
        # self.node_fields.sheets  = VN.vtk_to_numpy(data.GetCellData().GetArray('sheets'))
        # # Rodero mesh does not have normal vectors built in
        # self.node_fields.lvrv = list(VN.vtk_to_numpy(data.GetPointData().GetArray('V.dat')))
        # self.node_fields.tm = list(VN.vtk_to_numpy(data.GetPointData().GetArray('RHO.dat')))
        # self.node_fields.ab = list(VN.vtk_to_numpy(data.GetPointData().GetArray('Z.dat')))
        # self.node_fields.rt = list(VN.vtk_to_numpy(data.GetPointData().GetArray('PHI.dat')))

        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        # elem_ids = VN.vtk_to_numpy(data.GetCellData().GetArray('Ids'))
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('Ids')).astype(int), 'surface_node_ids', 'nodefield')
        self.element_fields.add_field(VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId')).astype(int), 'boundary_labels', 'elementfield')
        self.node_fields.add_field(VN.vtk_to_numpy(data.GetPointData().GetArray('RegionId')).astype(int), 'boundary_labels', 'nodefield')

        for i in range(0, len(self.element_fields.dict['boundary_labels'])):
            if self.element_fields.dict['boundary_labels'][i] == 0:
                self.element_fields.dict['boundary_labels'][i] = self.geometry.epicardium # Epicardial
            elif self.element_fields.dict['boundary_labels'][i] == 1:
                self.element_fields.dict['boundary_labels'][i] = self.geometry.lv_endocardium # LV endocardial
            elif self.element_fields.dict['boundary_labels'][i] == 2:
                self.element_fields.dict['boundary_labels'][i] = self.geometry.rv_endocardium # RV endocardial

        for i in range(0, len(self.node_fields.dict['boundary_labels'])):
            if self.node_fields.dict['boundary_labels'][i] == 0:
                self.node_fields.dict['boundary_labels'][i] = self.geometry.epicardium # Epicardial
            elif self.node_fields.dict['boundary_labels'][i] == 1:
                self.node_fields.dict['boundary_labels'][i] = self.geometry.lv_endocardium # LV endocardial
            elif self.node_fields.dict['boundary_labels'][i] == 2:
                self.node_fields.dict['boundary_labels'][i] = self.geometry.rv_endocardium # RV endocardial
        # Save input global labels
        input_node_labels_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        input_node_labels_global[self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['boundary_labels'] == self.geometry.lv_endocardium]] = self.geometry.lv_endocardium
        input_node_labels_global[self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['boundary_labels'] == self.geometry.rv_endocardium]] = self.geometry.rv_endocardium
        input_node_labels_global[self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['boundary_labels'] == self.geometry.epicardium]] = self.geometry.epicardium
        self.node_fields.add_field(input_node_labels_global, 'input_boundary_labels', 'nodefield')

        # Read faces
        self.geometry.number_of_triangles = data.GetNumberOfCells()
        self.geometry.triangles = np.zeros([self.geometry.number_of_triangles, 3])
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [int(self.node_fields.dict['surface_node_ids'][data.GetCell(i).GetPointId(0)]),
                                             int(self.node_fields.dict['surface_node_ids'][data.GetCell(i).GetPointId(1)]),
                                             int(self.node_fields.dict['surface_node_ids'][data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)

        # Set up boundary labels for EP and mechanics
        mechanical_element_boundary_labels = np.zeros(self.element_fields.dict['boundary_labels'].shape)
        ep_element_boundary_labels = np.zeros(self.element_fields.dict['boundary_labels'].shape)
        mechanical_node_boundary_labels = np.zeros(self.node_fields.dict['boundary_labels'].shape)
        ep_node_boundary_labels = np.zeros(self.node_fields.dict['boundary_labels'].shape)
        for i in range(0, len(self.geometry.triangles)):
            for j in range(0, 3):
                local_index = np.nonzero(self.node_fields.dict['surface_node_ids'] == self.geometry.triangles[i, j])[
                    0][0]
                if (self.node_fields.dict['boundary_labels'][local_index] == self.geometry.epicardium) & (self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10):
                    mechanical_node_boundary_labels[local_index] = self.geometry.valve_plug
                    ep_node_boundary_labels[local_index] = self.geometry.valve_plug
                    mechanical_element_boundary_labels[i] = self.geometry.valve_plug
                else:
                    mechanical_element_boundary_labels[local_index] = self.element_fields.dict['boundary_labels'][local_index]
                    ep_element_boundary_labels[local_index] = self.element_fields.dict['boundary_labels'][local_index]
                    mechanical_node_boundary_labels[local_index] = self.node_fields.dict['boundary_labels'][
                        local_index]
                    ep_node_boundary_labels[local_index] = self.node_fields.dict['boundary_labels'][local_index]
                if self.node_fields.dict['tm'][self.geometry.triangles[i, j]] == -10:
                    ep_node_boundary_labels[local_index] = self.geometry.valve_plug
                # Restrict mechanical epicardium to only 80% of apex-to-base axis.
                if (self.node_fields.dict['boundary_labels'][local_index] == self.geometry.epicardium) & \
                                (self.node_fields.dict['ab'][self.geometry.triangles[i, j]] > 0.8):
                # if (self.element_fields.dict['boundary_labels'][i] == self.geometry.epicardium) & \
                #         (self.node_fields.dict['ab'][self.geometry.triangles[i, j]] > 0.8):
                    mechanical_element_boundary_labels[i] = self.geometry.valve_plug  # Apply epicardial spring BC only to 80% of apex-to-base
                    mechanical_node_boundary_labels[local_index] = self.geometry.valve_plug
        self.node_fields.add_field(mechanical_node_boundary_labels, 'mechanical_node_boundary_labels', 'nodefield')
        self.node_fields.add_field(ep_node_boundary_labels, 'ep_node_boundary_labels', 'nodefield')
        self.element_fields.add_field(mechanical_element_boundary_labels, 'mechanical_element_boundary_labels', 'elementfield')
        self.element_fields.add_field(ep_element_boundary_labels, 'ep_element_boundary_labels', 'elementfield')

        # Get LV and RV endocardium nodes
        ep_lvnodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['ep_node_boundary_labels'] == self.geometry.lv_endocardium].astype(int)
        ep_rvnodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['ep_node_boundary_labels'] == self.geometry.rv_endocardium].astype(int)
        ep_epinodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['ep_node_boundary_labels'] == self.geometry.epicardium].astype(int)
        ep_valvenodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['ep_node_boundary_labels'] == self.geometry.valve_plug].astype(int)
        ep_node_labels_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        ep_node_labels_global[ep_lvnodes] = self.geometry.lv_endocardium
        ep_node_labels_global[ep_rvnodes] = self.geometry.rv_endocardium
        ep_node_labels_global[ep_epinodes] = self.geometry.epicardium
        ep_node_labels_global[ep_valvenodes] = self.geometry.valve_plug
        self.node_fields.add_field(ep_node_labels_global, 'ep_node_labels_global', 'nodefield')
        self.node_fields.add_field(ep_lvnodes, 'ep_lvnodes', 'nodefield')
        self.node_fields.add_field(ep_rvnodes, 'ep_rvnodes', 'nodefield')
        self.node_fields.add_field(ep_epinodes, 'ep_epinodes', 'nodefield')
        self.node_fields.add_field(ep_valvenodes, 'ep_valvenodes', 'nodefield')

        mechanical_lvnodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['mechanical_node_boundary_labels'] == self.geometry.lv_endocardium].astype(int)
        mechanical_rvnodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['mechanical_node_boundary_labels'] == self.geometry.rv_endocardium].astype(int)
        mechanical_epinodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['mechanical_node_boundary_labels'] == self.geometry.epicardium].astype(int)
        mechanical_valvenodes = self.node_fields.dict['surface_node_ids'][
            self.node_fields.dict['mechanical_node_boundary_labels'] == self.geometry.valve_plug].astype(int)
        mechanical_node_labels_global = np.zeros(self.geometry.number_of_nodes).astype(int)
        mechanical_node_labels_global[mechanical_lvnodes] = self.geometry.lv_endocardium
        mechanical_node_labels_global[mechanical_rvnodes] = self.geometry.rv_endocardium
        mechanical_node_labels_global[mechanical_epinodes] = self.geometry.epicardium
        mechanical_node_labels_global[mechanical_valvenodes] = self.geometry.valve_plug
        self.node_fields.add_field(mechanical_node_labels_global, 'mechanical_node_labels_global', 'nodefield')
        self.node_fields.add_field(mechanical_lvnodes, 'mechanical_lvnodes', 'nodefield')
        self.node_fields.add_field(mechanical_rvnodes, 'mechanical_rvnodes', 'nodefield')
        self.node_fields.add_field(mechanical_epinodes, 'mechanical_epinodes', 'nodefield')
        self.node_fields.add_field(mechanical_valvenodes, 'mechanical_valvenodes', 'nodefield')

        # Get LV and RV endocardium faces for Eikonal solution.
        ep_lvfaces = self.geometry.triangles[self.element_fields.dict['ep_element_boundary_labels'] == self.geometry.lv_endocardium].astype(int)
        ep_rvfaces = self.geometry.triangles[self.element_fields.dict['ep_element_boundary_labels'] == self.geometry.rv_endocardium].astype(int)
        self.element_fields.add_field(ep_lvfaces, 'ep_lvfaces', 'elementfield')
        self.element_fields.add_field(ep_rvfaces, 'ep_rvfaces', 'elementfield')


    def write_to_csv(self):
        print('Writing .csv files to '+self.output_dir)
        self.geometry.save_to_csv(self.output_dir)
        self.node_fields.save_to_csv(self.output_dir)
        self.element_fields.save_to_csv(self.output_dir)

    def write_to_ensight(self):
        print('Writing Ensight files to ' + self.output_dir)
        self.geometry.save_to_ensight(self.output_dir)
        self.node_fields.save_to_ensight(self.output_dir, casename=self.output_name+'_nodefields', geometry=self.geometry)
        self.element_fields.save_to_ensight(self.output_dir, casename=self.output_name+'_elementfields', geometry=self.geometry)