import vtk
import numpy as np
from vtk.util import numpy_support as VN

input_dir = '/users/jenang/RoderoNidererMeshHealthy/01/'
vtk_name = '01'
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName(input_dir + vtk_name + '.vtk')
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()
data = reader.GetOutput()

# Threshold based on

number_of_elements = data.GetNumberOfCells()
number_of_nodes = data.GetNumberOfPoints()
tetrahedrons = np.zeros((number_of_elements, 4)).astype(int)
nodes_xyz = np.zeros((number_of_nodes, 3))

