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
        self.element_fields.lvrv = VN.vtk_to_numpy(data.GetCellData().GetArray('tags'))
        self.node_fields.fibres = VN.vtk_to_numpy(data.GetCellData().GetArray('fiber'))
        self.node_fields.sheets  = VN.vtk_to_numpy(data.GetCellData().GetArray('sheet'))
        # Rodero mesh does not have normal vectors built in
        self.node_fields.lvrv = list(VN.vtk_to_numpy(data.GetPointData().GetArray('uvc_intraventricular')))
        self.node_fields.tm = list(VN.vtk_to_numpy(data.GetPointData().GetArray('uvc_transmural')))
        self.node_fields.ab = list(VN.vtk_to_numpy(data.GetPointData().GetArray('uvc_longitudinal')))
        self.node_fields.rt = list(VN.vtk_to_numpy(data.GetPointData().GetArray('uvc_rotational')))

        print('Reading vtk file from: ' + self.input_dir + self.vtk_name + '_surfaces_connectivity.vtk')
        reader = vtk.vtkPolyDataReader()
        reader.SetFileName(self.input_dir + self.vtk_name + '_surfaces_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        # elem_ids = VN.vtk_to_numpy(data.GetCellData().GetArray('Ids'))
        node_ids = VN.vtk_to_numpy(data.GetPointData().GetArray('Ids'))
        self.element_fields.boundaries = VN.vtk_to_numpy(data.GetCellData().GetArray('RegionId'))
        self.node_fields.boundaries = VN.vtk_to_numpy(data.GetPointData().GetArray('RegionId'))

        for i in range(0, len(self.element_fields.boundaries)):
            if self.element_fields.boundaries[i] == 0:
                self.element_fields.boundaries[i] = 1
            elif self.element_fields.boundaries[i] == 1:
                self.element_fields.boundaries[i] = 3
            elif self.element_fields.boundaries[i] == 2:
                self.element_fields.boundaries[i] = 2

        self.node_fields.global_boundaries = np.full(self.geometry.number_of_nodes, np.inf)
        for i in range(0, len(self.node_fields.boundaries)):
            self.node_fields.global_boundaries[node_ids[i]] = self.node_fields.boundaries[i]
        # Read faces
        self.geometry.number_of_triangles = data.GetNumberOfCells()
        self.geometry.triangles = np.zeros([self.geometry.number_of_triangles, 3])
        for i in range(0, self.geometry.number_of_triangles):
            self.geometry.triangles[i, :] = [int(node_ids[data.GetCell(i).GetPointId(0)]),
                                             int(node_ids[data.GetCell(i).GetPointId(1)]),
                                             int(node_ids[data.GetCell(i).GetPointId(2)])]
        self.geometry.triangles = self.geometry.triangles.astype(int)
        for i in range(0, len(self.geometry.triangles)):
            for j in range(0, 3):
                if self.node_fields.tm[self.geometry.triangles[i, j]] == -10:
                    if self.element_fields.boundaries[i] == 1:
                        self.element_fields.boundaries[i] = 4
                if (self.element_fields.boundaries[i] == 1) & \
                        (self.node_fields.ab[self.geometry.triangles[i, j]] > 0.8):
                    self.element_fields.boundaries[i] = 4  # Apply epicardial spring BC only to 80% of apex-to-base

    def write_to_csv(self):
        print('Writing .csv files to '+self.output_dir)
        self.geometry.save_to_csv(self.output_dir)
        self.node_fields.save_to_csv(self.output_dir)
        self.element_fields.save_to_csv(self.output_dir)

    def write_to_ensight(self):
        print('Writing Ensight files to ' + self.output_dir)
        self.geometry.save_to_ensight(self.output_dir)
        self.node_fields.save_to_ensight(self.output_dir, casename=self.output_name+'_nodefields')
        self.element_fields.save_to_ensight(self.output_dir, casename=self.output_name+'_elementfields')