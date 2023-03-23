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
    if len(var.shape) == 1:
        with open(dir + name + '.ensi.'+field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Scalar per-node variables file\npart\n         1\ncoordinates\n')
            for node_i in range(0, var.shape[0]):
                f.write(str(var[node_i]) + '\n')
    else:
        if var.shape[1] == 1:
            with open(dir + name + '.ensi.' + field_name, 'w') as f:
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
    if len(var.shape) == 1:
        with open(dir + name + '.ensi.'+field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Scalar per-cell variables file\npart\n         1\ncoordinates\n')
            for elem_i in range(0, var.shape[0]):
                f.write(str(var[elem_i]) + '\n')
    else:
        if var.shape[1] == 1:
            with open(dir + name + '.ensi.' + field_name, 'w') as f:
                f.write('Alya Ensight Gold --- Scalar per-cell variables file\npart\n         1\ncoordinates\n')
                for node_i in range(0, var.shape[0]):
                    f.write(str(var[node_i]) + '\n')
        with open(dir + name + '.ensi.' + field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Vector per-cell variables file\npart\n         1\ncoordinates\n')
            for c in range(var.shape[1]):
                for elem_i in range(var.shape[0]):
                    f.write(str(var[elem_i, c]) + '\n')

def save_ensight_geometry(dir, name, nodes_xyz, tetrahedrons):
    if np.amin(tetrahedrons) == 0:
        tetrahedrons = tetrahedrons + 1 # Ensight takes node indices starting from 1.
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

def save_ensight_case(dir, name, geometry_name, field_names, field_dimensions, field_type):
    with open(dir+name+'.ensi.case', 'w') as f:
        f.write('#\n# Alya generated post-process files\n# Ensight Gold Format\n#\n# Problem name:\t'+str(name)+'\n#\n')
        f.write('FORMAT\ntype:\tensight gold\nGEOMETRY\nmodel:\t1\t'+geometry_name+'.ensi.geo\nVARIABLE\n')
        for field_i in range(len(field_names)):
            field_name = field_names[field_i]
            field_dimension = field_dimensions[field_i]
            if field_type == 'nodefield':
                ensight_field_type = 'node'
            elif field_type == 'elementfield':
                ensight_field_type = 'cell'
            if field_dimension == 1:
                f.write(
                    'scalar per ' + ensight_field_type + ':\t' + str(field_dimension) + '\t' + field_name + '\t' + geometry_name + '.ensi.' + field_name + '\n')
            elif field_dimension > 1:
                f.write(
                    'vector per ' + ensight_field_type + ':\t' + str(field_dimension) + '\t' + field_name + '\t' + geometry_name + '.ensi.' + field_name + '\n')


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
        savetxt(filename=output_dir + self.name + '_xyz.csv', var=self.nodes_xyz)
        savetxt(filename=output_dir + self.name + '_tri.csv', var=self.tetrahedrons)
        savetxt(filename=output_dir + self.name + '_triangles.csv', var=self.triangles)
        savetxt(filename=output_dir + self.name + '_tetrahedron_centers.csv', var=self.tetrahedron_centres)
        savetxt(filename=output_dir + self.name + '_edges.csv', var=self.edges)


    def save_to_ensight(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        save_ensight_geometry(dir=output_dir, name=self.name, nodes_xyz=self.nodes_xyz, tetrahedrons=self.tetrahedrons)

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        self.nodes_xyz = loadtxt(filename=input_dir+self.name+'_xyz.csv')
        self.tetrahedrons = loadtxt(filename=input_dir + self.name + '_tri.csv').astype(int)
        if np.amin(self.tetrahedrons == 1):
            self.tetrahedrons = self.tetrahedrons - 1
        self.number_of_elements = self.tetrahedrons.shape[0]
        self.number_of_nodes = self.nodes_xyz.shape[0]
        self.triangles = loadtxt(filename=input_dir + self.name + '_triangles.csv')
        self.number_of_triangles = self.triangles.shape[0]
        self.tetrahedron_centres = loadtxt(filename=input_dir + self.name + '_tetrahedron_centres.csv')
        self.edges = loadtxt(filename=input_dir + self.name + '_edges.csv')


class Fields:
    def __init__(self, name, field_type='nodefield'):
        self.name = name
        self.field_type = field_type
        self.dict = {} # Placeholder.
        self.number_of_fields = len(self.dict)

    def add_field(self, data, data_name, field_type):
        assert field_type == self.field_type, 'Input field type: '+field_type+' does not matching existing field type: '+field_type
        self.dict[data_name] = data

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        list_fields = list(self.dict.keys())
        self.number_of_fields = len(list_fields)
        for field_i in range(len(list_fields)):
            varname = list_fields[field_i]
            savetxt(filename=output_dir + self.name + '_'+self.field_type+'_' + varname + '.csv', var=self.dict[varname])


    def save_to_ensight(self, output_dir, casename, geometry):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        list_fields = list(self.dict.keys())
        field_dimensions = []
        list_fields_output = []
        for field_i in range(self.number_of_fields):
            varname = list_fields[field_i]
            if (self.field_type == 'nodefield') & (self.dict[varname].shape[0] == geometry.number_of_nodes): # Don't write out fields that aren't meant for this geometry.
                save_ensight_node(dir=output_dir, name=self.name, field_name=varname, var=self.dict[varname])
                if len(self.dict[varname].shape) == 1:
                    field_dimensions.append(1)
                else:
                    field_dimensions.append(self.dict[varname].shape[1])
                list_fields_output.append(varname)
            elif (self.field_type == 'elementfield') & (self.dict[varname].shape[0] == geometry.number_of_elements):  # Don't write out fields that aren't meant for this geometry.
                save_ensight_element(dir=output_dir, name=self.name, field_name=varname, var=self.dict[varname])
                if len(self.dict[varname].shape) == 1:
                    field_dimensions.append(1)
                else:
                    field_dimensions.append(self.dict[varname].shape[1])
                list_fields_output.append(varname)
            # else:
            #     raise ValueError('save_to_ensight: Field type '+self.field_type+' not found.')
        save_ensight_case(dir=output_dir, name=casename, geometry_name=geometry.name, field_names=list_fields_output, field_dimensions=field_dimensions, field_type = self.field_type)

    def read_csv_to_attributes(self, input_dir, field_type):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        filenames = np.array([f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f)) and '_'+field_type+'_' in f])
        for file_i in range(filenames.shape[0]):
            print('Reading in '+input_dir+filenames[file_i])
            varname = filenames[file_i].split('_')[-1].split('.')[0]
            self.dict[varname] = loadtxt(filename=input_dir+ filenames[file_i])
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
