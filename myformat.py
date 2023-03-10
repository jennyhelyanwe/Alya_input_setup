# import pandas as pd
import os
import numpy as np


def savetxt(filename, var):
    if not var is None:
        np.savetxt(filename, var, delimiter=',')

def loadtxt(filename):
    if os.path.exists(filename):
        return np.loadtxt(filename, delimiter=',')
    else:
        return None


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

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        savetxt(filename=output_dir + self.name + '_nodes_xyz.csv', var=self.nodes_xyz)
        savetxt(filename=output_dir + self.name + '_tetrahedrons.csv', var=self.tetrahedrons)
        savetxt(filename=output_dir + self.name + '_triangles.csv', var=self.triangles)
        savetxt(filename=output_dir + self.name + '_tetrahedron_centres.csv', var=self.tetrahedron_centres)
        savetxt(filename=output_dir + self.name + '_edges.csv', var=self.edges)

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.nodes_xyz = loadtxt(filename=input_dir+self.name+'_nodes_xyz.csv')
        self.tetrahedrons = loadtxt(filename=input_dir + self.name + '_tetrahedrons.csv')
        self.triangles = loadtxt(filename=input_dir + self.name + '_triangles.csv')
        self.tetrahedron_centres = loadtxt(filename=input_dir + self.name + '_tetrahedron_centres.csv')
        self.edges = loadtxt(filename=input_dir + self.name + '_edges.csv')


class Fields:
    def __init__(self, name, field_type='nodefield'):
        self.name = name
        self.field_type = field_type
        self.dict = {'celltype': np.zeros((10))} # Placeholder.
        self.number_of_fields = len(self.dict)

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        list_fields = list(self.dict.keys())
        for field_i in range(self.number_of_fields):
            varname = list_fields[field_i]
            savetxt(filename=output_dir + self.name + '_'+self.field_type+'_' + varname + '.csv', var=self.dict[varname])

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        filenames = np.array([f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f)) and 'nodefield' in f])
        for file_i in range(filenames.shape[0]):
            varname = filenames[file_i].split('_')[-1].split('.')[0]
            self.dict[varname] = loadtxt(filename=filenames[file_i])
        self.number_of_fields = len(self.dict)


class Materials:
    def __init__(self, name):
        self.name = name
        self.materials = None

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'

        savetxt(filename=output_dir + self.name + '_materials.csv', var=self.materials)

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.materials = loadtxt(filename=input_dir + self.name + '_materials.csv')


if __name__ == '__main__':
    wd = os.getcwd()
    geometry = Geometry('test')
    geometry.save_to_csv(wd)
    node_fields = Fields('test', field_type='nodefield')
    node_fields.save_to_csv(wd)
    element_fields = Fields('test', field_type='elementfield')
    element_fields.save_to_csv(wd)
    materials = Materials('test')
    materials.save_to_csv(wd)
