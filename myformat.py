# import pandas as pd
import os
import numpy as np


class Geometry:
    def __init__(self, name):
        # Geometry initialisations
        self.name = name
        self.number_of_nodes = 0
        self.number_of_elements = 0
        self.number_of_triangles = 0
        self.nodes_xyz = None
        self.tetrahedrons = None
        self.triangles = None
        self.tetrahedron_centres = None
        self.edges = None

    def _savetxt(self, filename, var):
        if not var is None:
            np.savetxt(filename, var, delimiter=',')

    def _loadtxt(self, filename):
        if os.path.exists(filename):
            return np.loadtxt(filename, delimiter=',')
        else:
            return None

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        self._savetxt(filename=output_dir + self.name + '_nodes_xyz.csv', var=self.nodes_xyz)
        self._savetxt(filename=output_dir + self.name + '_tetrahedrons.csv', var=self.tetrahedrons)
        self._savetxt(filename=output_dir + self.name + '_triangles.csv', var=self.triangles)
        self._savetxt(filename=output_dir + self.name + '_tetrahedron_centres.csv', var=self.tetrahedron_centres)
        self._savetxt(filename=output_dir + self.name + '_edges.csv', var=self.edges)

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.nodes_xyz = self._loadtxt(filename=input_dir+self.name+'_nodes_xyz.csv')
        self.tetrahedrons = self._loadtxt(filename=input_dir + self.name + '_tetrahedrons.csv')
        self.triangles = self._loadtxt(filename=input_dir + self.name + '_triangles.csv')
        self.tetrahedron_centres = self._loadtxt(filename=input_dir + self.name + '_tetrahedron_centres.csv')
        self.edges = self._loadtxt(filename=input_dir + self.name + '_edges.csv')


class NodeFields:
    def __init__(self, name):
        self.name = name
        # Field initialisations
        self.fibre = None
        self.sheet = None
        self.normal = None
        self.lvrv = None
        self.tm = None
        self.ab = None
        self.rt = None
        self.celltype = None
        self.ab_Gks_sf = None
        self.stimulus = None
        self.boundaries = None

    def _savetxt(self, filename, var):
        if not var is None:
            np.savetxt(filename, var, delimiter=',')

    def _loadtxt(self, filename):
        if os.path.exists(filename):
            return np.loadtxt(filename, delimiter=',')
        else:
            return None

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_fibre.csv', var=self.fibre)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_sheet.csv', var=self.sheet)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_normal.csv', var=self.normal)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_lvrv.csv', var=self.lvrv)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_tm.csv', var=self.tm)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_ab.csv', var=self.ab)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_rt.csv', var=self.rt)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_celltype.csv', var=self.celltype)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_ab_Gks_sf.csv', var=self.ab_Gks_sf)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_stimulus.csv', var=self.stimulus)
        self._savetxt(filename=output_dir + self.name + '_nodefield' + '_boundaries.csv', var=self.boundaries)

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.fibre = self._loadtxt(filename=input_dir+self.name + '_nodefield' + '_fibre.csv')
        self.sheet = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_sheet.csv')
        self.normal = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_normal.csv')
        self.lvrv = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_lvrv.csv')
        self.tm = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_tm.csv')
        self.ab = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_ab.csv')
        self.rt = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_rt.csv')
        self.celltype = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_celltype.csv')
        self.ab_Gks_sf = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_ab_Gks_sf.csv')
        self.stimulus = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_stimulus.csv')
        self.boundaries = self._loadtxt(filename=input_dir + self.name + '_nodefield' + '_boundaries.csv')


class ElementFields:
    def __init__(self, name):
        self.name = name
        self.materials = None
        self.lvrv = None
        self.boundaries = None

    def _savetxt(self, filename, var):
        if not var is None:
            np.savetxt(filename, var, delimiter=',')

    def _loadtxt(self, filename):
        if os.path.exists(filename):
            return np.loadtxt(filename, delimiter=',')
        else:
            return None

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'

        self._savetxt(filename=output_dir + self.name + '_elementfield' + '_materials.csv', var=self.materials)
        self._savetxt(filename=output_dir + self.name + '_elementfield' + '_lvrv.csv', var=self.lvrv)
        self._savetxt(filename=output_dir + self.name + '_elementfield' + '_boundaries.csv', var=self.boundaries)

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.materials = self._loadtxt(filename=input_dir+self.name + '_elementfield' + '_materials.csv')
        self.lvrv = self._loadtxt(filename=input_dir + self.name + '_elementfield' + '_lvrv.csv')
        self.boundaries = self._loadtxt(filename=input_dir + self.name + '_elementfield' + '_boundaries.csv')


if __name__ == '__main__':
    wd = os.getcwd()
    geometry = Geometry('test')
    geometry.save_to_csv(wd)
    node_fields = NodeFields('test')
    node_fields.save_to_csv(wd)
    element_fields = ElementFields('test')
    element_fields.save_to_csv(wd)
