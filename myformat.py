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

def save_ensight_node(dir, name, field_name, var):
    print(var)
    print(var.shape)
    if var.shape[1] == 1:
        with open(dir + name + '.ensi.'+field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Scalar per-node variables file\npart\n         1\ncoordinates\n')
            for node_i in range(0, var.shape[0]):
                f.write(str(var[node_i]) + '\n')
    else:
        with open(dir + name + '.ensi.' + field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Vector per-node variables file\npart\n         1\ncoordinates\n')
            for c in range(var.shape[1]):
                for i in range(var.shape[0]):
                    f.write(str(var[i, c]) + '\n')

def save_ensight_element(dir, name, field_name, var):
    if var.shape[1] == 1:
        with open(dir + name + '.ensi.'+field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Scalar per-cell variables file\npart\n         1\ncoordinates\n')
            for elem_i in range(0, var.shape[0]):
                f.write(str(var[elem_i]) + '\n')
    else:
        with open(dir + name + '.ensi.' + field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Vector per-cell variables file\npart\n         1\ncoordinates\n')
            for c in range(var.shape[1]):
                for elem_i in range(var.shape[0]):
                    f.write(str(var[elem_i, c]) + '\n')

def save_ensight_geometry(dir, name, nodes_xyz, tetrahedrons):
    with open(dir+name+'.ensi.geo', 'w') as f:
        f.write(
            'Problem name:  '+str(name)+'\nGeometry file\nnode id given\nelement id given\npart\n\t1\nVolume Mesh\ncoordinates\n' + str(
                nodes_xyz.shape[0]) + '\n')
        for node_i in range(nodes_xyz.shape[0]):
            f.write(str(node_i + 1) + '\n')
        for c in [0, 1, 2]:
            for i in range(nodes_xyz.shape[0]):
                f.write(str(nodes_xyz[i, c]) + '\n')
        f.write('tetra4\n  ' + str(tetrahedrons.shape[0]) + '\n')
        for elem_i in range(tetrahedrons.shape[0]):
            f.write('  ' + str(elem_i + 1) + '\n')
        for i in range(tetrahedrons.shape[0]):
            f.write(str(tetrahedrons[i, 0]) + '\t' + str(tetrahedrons[i, 1]) + '\t' + str(tetrahedrons[i, 2]) + '\t' + str(
                tetrahedrons[i, 3]) + '\n')

def save_ensight_case(dir, name, field_names, field_dimensions, field_type):
    with open(dir+name+'.ensi.case', 'w') as f:
        f.write('#\n# Alya generated post-process files\n# Ensight Gold Format\n#\n# Problem name:\t'+str(name)+'\n#\n')
        f.write('FORMAT\ntype:\tensight gold\nGEOMETRY\nmodel:\t1\t'+name+'.ensi.geo\nVARIABLE\n')
        for field_i in range(len(field_names)):
            field_name = field_names[field_i]
            field_dimension = field_dimensions[field_i]
            if field_type == 'nodefield':
                ensight_field_type = 'node'
            elif field_type == 'elementfield':
                ensight_field_type = 'cell'
            if field_dimension == 1:
                f.write(
                    'scalar per ' + ensight_field_type + ':\t' + field_dimension + '\t' + field_name + '\t' + name + '.ensi.' + field_name + '\n')
            elif field_dimension > 1:
                f.write(
                    'vector per ' + ensight_field_type + ':\t' + field_dimension + '\t' + field_name + '\t' + name + '.ensi.' + field_name + '\n')


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


    def save_to_ensight(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        save_ensight_geometry(dir=output_dir, name=self.name, nodes_xyz=self.nodes_xyz, tetrahedrons=self.tetrahedrons)

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

    def add_field(self, data, data_name, field_type):
        assert field_type == self.field_type, 'Input field type: '+field_type+' does not matching existing field type: '+field_type
        self.dict[data_name] = data

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        list_fields = list(self.dict.keys())
        for field_i in range(self.number_of_fields):
            varname = list_fields[field_i]
            savetxt(filename=output_dir + self.name + '_'+self.field_type+'_' + varname + '.csv', var=self.dict[varname])


    def save_to_ensight(self, output_dir, casename):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        list_fields = list(self.dict.keys())
        field_dimensions = []
        for field_i in range(self.number_of_fields):
            varname = list_fields[field_i]
            if self.field_type == 'nodefield':
                save_ensight_node(dir=output_dir, name=self.name, field_name=varname, var=self.dict[varname])
                field_dimensions.append(self.dict[varname].shape[1])
            elif self.field_type == 'elementfield':
                save_ensight_element(dir=output_dir, name=self.name, field_name=varname, var=self.dict[varname])
                field_dimensions.append(self.dict[varname].shape[1])
            else:
                raise ValueError('save_to_ensight: Field type '+self.field_type+' not found.')
        save_ensight_case(dir=output_dir, name=casename, field_names=list_fields, field_dimensions=field_dimensions, field_type = self.field_type)

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
