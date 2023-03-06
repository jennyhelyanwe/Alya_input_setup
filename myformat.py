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

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        if self.nodes_xyz:
            np.savetxt(output_dir + self.name+'_nodes_xyz.csv', self.nodes_xyz, delimiter=',')
        if self.tetrahedrons:
            np.savetxt(output_dir + self.name+'_tetrahedrons.csv', self.tetrahedrons, delimiter=',')
        if self.triangles:
            np.savetxt(output_dir + self.name + '_triangles.csv', self.triangles, delimiter=',')

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.nodes_xyz = np.loadtxt(input_dir+self.name+'_nodes_xyz.csv', delimiter=',')
        self.tetrahedrons = np.loadtxt(input_dir + self.name + '_tetrahedrons.csv', delimiter=',')
        self.triangles = np.loadtxt(input_dir + self.name + '_triangles.csv', delimiter=',')


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

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        np.savetxt(output_dir + self.name + '_nodefield' + '_fibre.csv', self.fibre, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_sheet.csv', self.sheet, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_normal.csv', self.normal, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_lvrv.csv', self.lvrv, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_tm.csv', self.tm, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_ab.csv', self.ab, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_rt.csv', self.rt, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_celltype.csv', self.celltype, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_ab_Gks_sf.csv', self.ab_Gks_sf, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_stimulus.csv', self.stimulus, delimiter=',')
        np.savetxt(output_dir + self.name + '_nodefield' + '_boundaries.csv', self.boundaries, delimiter=',')

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.fibre = np.loadtxt(input_dir+self.name + '_nodefield' + '_fibre.csv', delimiter=',')
        self.sheet = np.loadtxt(input_dir + self.name + '_nodefield' + '_sheet.csv', delimiter=',')
        self.normal = np.loadtxt(input_dir + self.name + '_nodefield' + '_normal.csv', delimiter=',')
        self.lvrv = np.loadtxt(input_dir + self.name + '_nodefield' + '_lvrv.csv', delimiter=',')
        self.tm = np.loadtxt(input_dir + self.name + '_nodefield' + '_tm.csv', delimiter=',')
        self.ab = np.loadtxt(input_dir + self.name + '_nodefield' + '_ab.csv', delimiter=',')
        self.rt = np.loadtxt(input_dir + self.name + '_nodefield' + '_rt.csv', delimiter=',')
        self.celltype = np.loadtxt(input_dir + self.name + '_nodefield' + '_celltype.csv', delimiter=',')
        self.ab_Gks_sf = np.loadtxt(input_dir + self.name + '_nodefield' + '_ab_Gks_sf.csv', delimiter=',')
        self.stimulus = np.loadtxt(input_dir + self.name + '_nodefield' + '_stimulus.csv', delimiter=',')
        self.boundaries = np.loadtxt(input_dir + self.name + '_nodefield' + '_boundaries.csv', delimiter=',')


class ElementFields:
    def __init__(self, name):
        self.name = name
        self.materials = None
        self.lvrv = None
        self.boundaries = None

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        np.savetxt(output_dir + self.name + '_elementfield' + '_materials.csv', self.materials, delimiter=',')
        np.savetxt(output_dir + self.name + '_elementfield' + '_lvrv.csv', self.lvrv, delimiter=',')
        np.savetxt(output_dir + self.name + '_elementfield' + '_boundaries.csv', self.boundaries, delimiter=',')

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.materials = np.loadtxt(input_dir+self.name + '_elementfield' + '_materials.csv', delimiter=',')
        self.lvrv = np.loadtxt(input_dir + self.name + '_elementfield' + '_lvrv.csv', delimiter=',')
        self.boundaries = np.loadtxt(input_dir + self.name + '_elementfield' + '_boundaries.csv', delimiter=',')


if __name__ == '__main__':
    wd = os.getcwd()
    geometry = Geometry('test')
    geometry.save_to_csv(wd)
    node_fields = NodeFields('test')
    node_fields.save_to_csv(wd)
    element_fields = ElementFields('test')
    element_fields.save_to_csv(wd)
