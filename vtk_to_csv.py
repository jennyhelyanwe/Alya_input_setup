import vtk
import numpy as np
import pandas as pd
from myformat import *
from vtk.util import numpy_support as VN


class VTKtoCSV:
    def __init__(self, name, input_dir, output_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        self.name = name
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.geometry = Geometry(self.name)
        self.node_fields = NodeFields(self.name)
        self.element_fields = ElementFields(self.name)
        self.surface_fields = SurfaceFields(self.name)
        self.read_vtk()
        self.write_to_csv()

    def read_vtk(self):
        print('Reading vtk file from: ' + self.input_dir + self.name + '.vtk')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(self.input_dir + self.name + '.vtk')
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
        self.element_fields.lvrv = VN.vtk_to_numpy(data.GetCellData().GetArray('ID'))
        self.node_fields.fibres = VN.vtk_to_numpy(data.GetCellData().GetArray('fibres'))
        self.node_fields.sheets  = VN.vtk_to_numpy(data.GetCellData().GetArray('sheets'))
        # Rodero mesh does not have normal vectors built in
        self.node_fields.lvrv = list(VN.vtk_to_numpy(data.GetPointData().GetArray('V.dat')))
        self.node_fields.tm = list(VN.vtk_to_numpy(data.GetPointData().GetArray('RHO.dat')))
        self.node_fields.ab = list(VN.vtk_to_numpy(data.GetPointData().GetArray('Z.dat')))
        self.node_fields.rt = list(VN.vtk_to_numpy(data.GetPointData().GetArray('PHI.dat')))

        reader = vtkPolyDataReader()
        reader.SetFileName(self.input_dir + self.name + '_surface_connectivity.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        elem_ids = VN.vtk_to_numpy(data.GetCellData().GetArray('Ids'))
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