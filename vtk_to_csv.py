import vtk
import numpy as np
import pandas as pd
from myformat import *
from vtk.util import numpy_support as VN

class vtk_to_csv():
    def __init__(self, name, input_dir, output_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        self.name = name
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.geometry = geometry(self.name)
        self.node_fields = node_fields(self.name)
        self.element_fields = element_fields(self.name)
        self.surface_fields = surface_fields(self.name)

        self.read_vtk()
        self.write_to_csv()

    def read_vtk(self):
        print('Reading vtk file from: '+ self.input_dir + self.name + '.vtk')
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(self.input_dir + self.name + '.vtk')
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        self.geometry.number_of_elements = data.GetNumberOfCells()
        self.geometry.number_of_nodes = data.GetNumberOfPoints()
        self.geometry.tetrahedrons_1 = np.zeros((self.geometry.number_of_elements,)).astype(int)
        self.geometry.tetrahedrons_2 = np.zeros((self.geometry.number_of_elements,)).astype(int)
        self.geometry.tetrahedrons_3 = np.zeros((self.geometry.number_of_elements,)).astype(int)
        self.geometry.tetrahedrons_4 = np.zeros((self.geometry.number_of_elements,)).astype(int)
        self.geometry.nodes_xyz_1 = np.zeros((self.geometry.number_of_nodes,))
        self.geometry.nodes_xyz_2 = np.zeros((self.geometry.number_of_nodes,))
        self.geometry.nodes_xyz_3 = np.zeros((self.geometry.number_of_nodes,))
        for i in range(0, self.geometry.number_of_elements):
            self.geometry.tetrahedrons_1[i] = int(data.GetCell(i).GetPointId(0))
            self.geometry.tetrahedrons_2[i] = int(data.GetCell(i).GetPointId(1))
            self.geometry.tetrahedrons_3[i] = int(data.GetCell(i).GetPointId(2))
            self.geometry.tetrahedrons_4[i] = int(data.GetCell(i).GetPointId(3))
        self.geometry.tetrahedrons_1 = list(self.geometry.tetrahedrons_1)
        self.geometry.tetrahedrons_2 = list(self.geometry.tetrahedrons_2)
        self.geometry.tetrahedrons_3 = list(self.geometry.tetrahedrons_3)
        self.geometry.tetrahedrons_4 = list(self.geometry.tetrahedrons_4)
        nodes_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())  # Convert from mm to cm
        self.geometry.nodes_xyz_1 = list(nodes_xyz[:, 0])
        self.geometry.nodes_xyz_2 = list(nodes_xyz[:, 1])
        self.geometry.nodes_xyz_3 = list(nodes_xyz[:, 2])
        self.element_fields.lvrv = VN.vtk_to_numpy(data.GetCellData().GetArray('ID'))
        fibres = VN.vtk_to_numpy(data.GetCellData().GetArray('fibres'))
        self.node_fields.fibres_1 = list(fibres[:, 0])
        self.node_fields.fibres_2 = list(fibres[:, 1])
        self.node_fields.fibres_3 = list(fibres[:, 2])
        sheets = VN.vtk_to_numpy(data.GetCellData().GetArray('sheets'))
        self.node_fields.sheets_1 = list(sheets[:, 0])
        self.node_fields.sheets_2 = list(sheets[:, 1])
        self.node_fields.sheets_3 = list(sheets[:, 2])
        # Rodero mesh does not have normal vectors built in
        self.node_fields.lvrv = list(VN.vtk_to_numpy(data.GetPointData().GetArray('V.dat')))
        self.node_fields.tm = list(VN.vtk_to_numpy(data.GetPointData().GetArray('RHO.dat')))
        self.node_fields.ab = list(VN.vtk_to_numpy(data.GetPointData().GetArray('Z.dat')))
        self.node_fields.rt = list(VN.vtk_to_numpy(data.GetPointData().GetArray('PHI.dat')))

    def write_to_csv(self):
        print('Writing .csv files to '+self.output_dir)
        self.geometry.save_to_csv(self.output_dir)
        self.node_fields.save_to_csv(self.output_dir)
        self.element_fields.save_to_csv(self.output_dir)