import os, sys
import numpy as np
import vtk
from vtk.util import numpy_support as VN


class I():
    def __init__(self, name, input_dir, input_format, save_dir, refresh ):
        self.name = name
        self.input_dir = input_dir
        self.input_format = input_format
        self.save_dir = save_dir
        self.refresh = refresh

    def _read_mesh(self):
        if (not self.refresh) & (os.path.exists(self.save_dir+'nodes.npy')):
            print('# SAVED: Reading in mesh '+self.name+' with wall axes')
            self.nodes = np.load(self.save_dir+'nodes.npy')
            self.elems = np.load(self.save_dir+'elems.npy')
            self.longitudinal_coordinate = np.load(self.save_dir+'longitudinal_coordinate.npy')
            self.transmural_coordinate = np.load(self.save_dir+'transmural_coordinate.npy')
            self.rotational_coordinate = np.load(self.save_dir+'rotational_coordinate.npy')
            self.longitudinal_vector = np.load(self.save_dir+'longitudinal_vector.npy')
            self.transmural_vector = np.load(self.save_dir+'transmural_vector.npy')
            self.rotational_vector = np.load(self.save_dir+'rotational_vector.npy')
            self.nnodes = len(self.nodes)
            self.nelems = len(self.elems)
        else:
            print('# NEW: Reading in mesh '+self.name)
            if self.input_format == 'vtk':
                reader = vtk.vtkUnstructuredGridReader()
                reader.SetFileName(self.input_dir+self.name+'.vtk')
                reader.ReadAllVectorsOn()
                reader.ReadAllScalarsOn()
                reader.Update()
                self.data = reader.GetOutput()
                self.nelems = self.data.GetNumberOfCells()
                self.nnodes = self.data.GetNumberOfPoints()
                self.elems = np.zeros((self.nelems, 4)).astype(int)
                for i in range(0, self.nelems):
                    for j in range(0, 4):
                        self.elems[i, j] = int(data.GetCell(i).GetPointId(j))+1
                self.nodes = VN.vtk_to_numpy(data.GetPoints().GetData()) # Convert from mm to cm
                np.ssave(self.save_dir+'elems', self.elems)
                np.save(self.save_dir+'nodes', self.nodes)
                self.lvrv_elems = VN.vtk_to_numpy(self.data.GetCellData().GetArray('ID'))
                # self.fibres = VN.vtk_to_numpy(self.data.GetCellData().GetArray('fibres'))
                # self.sheets = VN.vtk_to_numpy(self.data.GetCellData().GetArray('sheets'))
                self.lvrv_nodes = VN.vtk_to_numpy(self.data.GetPointData().GetArray('V.dat'))
                self.transmural_coordinate = VN.vtk_to_numpy(self.data.GetPointData().GetArray('RHO.dat'))
                self.longitudinal_coordinate = VN.vtk_to_numpy(self.data.GetPointData().GetArray('Z.dat'))
                self.rotational_coordinate = VN.vtk_to_numpy(self.data.GetPointData().GetArray('PHI.dat'))
            elif self.input_format == 'cobiveco':
                reader = vtk.vtkXMLUnstructuredGridReader()
                reader.SetFileName(self.input_dir+self.name+'.vtu')
                reader.Update()
                self.data = reader.GetOutput()
                self.nelems = self.data.GetNumberOfCells()
                self.nnodes = self.data.GetNumberOfPoints()
                self.elems = np.zeros((self.nelems, 4)).astype(int)
                for i in range(0, self.nelems):
                    for j in range(0, 4):
                        self.elems[i, j] = int(self.data.GetCell(i).GetPointId(j))+1
                self.nodes = VN.vtk_to_numpy(self.data.GetPoints().GetData()) # Convert from mm to cm
                np.save(self.save_dir+'elems', self.elems)
                np.save(self.save_dir+'nodes', self.nodes)
                self.longitudinal_coordinate = VN.vtk_to_numpy(self.data.GetPointData().GetArray('ab'))
                self.transmural_coordinate = VN.vtk_to_numpy(self.data.GetPointData().GetArray('tm'))
                self.rotational_coordinate = VN.vtk_to_numpy(self.data.GetPointData().GetArray('rt'))
            else:
                sys.exit('Input format '+self.mesh_input_format+' not recognised. ')
            grad = vtk.vtkGradientFilter()
            grad.SetInputData(self.data)
            grad.SetInputArrayToProcess(0,0,0,0,'ab')
            grad.SetResultArrayName('ab_grad')
            grad.Update()
            self.longitudinal_vector = VN.vtk_to_numpy(grad.GetOutput().GetPointData().GetArray('ab_grad'))
            grad.SetInputArrayToProcess(0,0,0,0,'tm')
            grad.SetResultArrayName('tm_grad')
            grad.Update()
            self.transmural_vector = VN.vtk_to_numpy(grad.GetOutput().GetPointData().GetArray('tm_grad'))
            grad.SetInputArrayToProcess(0,0,0,0,'rt')
            grad.SetResultArrayName('rt_grad')
            grad.Update()
            self.rotational_vector = VN.vtk_to_numpy(grad.GetOutput().GetPointData().GetArray('rt_grad'))
            np.save(self.save_dir+'longitudinal_coordinate',self.longitudinal_coordinate)
            np.save(self.save_dir+'transmural_coordinate',self.transmural_coordinate)
            np.save(self.save_dir+'rotational_coordinate',self.rotational_coordinate)
            np.save(self.save_dir+'longitudinal_vector',self.longitudinal_vector)
            np.save(self.save_dir+'transmural_vector',self.transmural_vector)
            np.save(self.save_dir+'rotational_vector',self.rotational_vector)

        print('Number of nodes: '+str(self.nnodes)+', Number of elements: '+str(self.nelems))


    def _read_surfaces(self):

        if ((self.refresh == False) & (os.path.exists(self.save_dir+'faces.npy'))):
            print('# SAVED: Reading surfaces and labels.')
            self.faces = np.load(self.save_dir+'faces.npy')
            self.faces_label = np.load(self.save_dir+'faces_label.npy')
            self.facesnodes_label = np.load(self.save_dir+'facesnodes_label.npy')
        else:
            print('# NEW: Reading surfaces and labels.')
            if self.input_format == 'cobiveco':
                surf = vtk.vtkXMLPolyDataReader()
                surf.SetFileName(self.input_dir+self.name+'.vtp')
                surf.Update()
                data = surf.GetOutput()
                self.surface_nodes = VN.vtk_to_numpy(data.GetPoints().GetData())
                # Find correspondence between face nodes and mesh nodes
                surface_idx_to_mesh_idx = np.zeros((len(self.surface_nodes),)).astype(int)
                mins = np.zeros((len(self.surface_nodes),))
                for i in range(0, len(self.surface_nodes)):
                    ndist = np.linalg.norm((self.nodes - self.surface_nodes[i,:]), axis=1)
                    surface_idx_to_mesh_idx[i] = np.where(ndist == np.amin(ndist))[0][0]
                    mins[i] = np.amin(ndist)
                cobiveco_faces_node_label = VN.vtk_to_numpy(data.GetPointData().GetArray('class'))
                alya_faces_node_label = np.zeros((len(cobiveco_faces_node_label),)).astype(int)
                for i in range(0, len(cobiveco_faces_node_label)):
                    # Cobiveco vs Alya surface classes are as follows:
                    # 1 - lid  |  epi
                    # 2 - epi  |  LV endo
                    # 3 - LV endo | RV endo
                    # 4 - RV endo | lid
                    if cobiveco_faces_node_label[i] == 1:
                        alya_faces_node_label[i] = 4
                    elif cobiveco_faces_node_label[i] == 2:
                        alya_faces_node_label[i] = 1
                    elif cobiveco_faces_node_label[i] == 3:
                        alya_faces_node_label[i] = 2
                    elif cobiveco_faces_node_label[i] == 4:
                        alya_faces_node_label[i] = 3
                self.nfaces = data.GetNumberOfCells()
                self.faces = np.zeros((self.nfaces, 3)).astype(int)
                self.faces_label = np.zeros((self.nfaces,)).astype(int)
                for i in range(0, self.nfaces):
                    node_idxs = []
                    for j in range(0, 3):
                        node_idx = data.GetCell(i).GetPointId(j)
                        self.faces[i, j] = int(surface_idx_to_mesh_idx[node_idx])+1
                        node_idxs.append(data.GetCell(i).GetPointId(j))
                    label = int(np.floor(np.mean(alya_faces_node_label[node_idxs].astype(float))))
                self.facesnodes_label = np.zeros((self.nnodes,)).astype(int)
                self.facesnodes_label[surface_idx_to_mesh_idx] = alya_faces_node_label
                np.save(self.save_dir+'faces', self.faces)
                np.save(self.save_dir+'faces_label', self.faces_label)
                np.save(self.save_dir+'facesnodes_label', self.facesnodes_label)


    def _read_endocardial_activation(self):
        print('Reading in endocardial activation is currently not implemented...')
