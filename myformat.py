import os
import numpy as np
import pymp, multiprocessing

class Geometry:
    def __init__(self, name, verbose):
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
        self.lv_endocardium = 3
        self.rv_endocardium = 2
        self.epicardium = 1
        self.valve_plug = 4
        self.ra_endocardium = 5
        self.la_endocardium = 6
        self.lv = -1
        self.rv = 1
        self.base = 1
        self.apex = 0
        self.tm_endo = 0
        self.tm_epi = 1
        self.pericardial_ab_extent = 0.75
        self.verbose = verbose

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        save_txt(filename=output_dir + self.name + '_xyz.csv', var=self.nodes_xyz)
        save_txt(filename=output_dir + self.name + '_tetra.csv', var=self.tetrahedrons)
        save_txt(filename=output_dir + self.name + '_triangles.csv', var=self.triangles)
        save_txt(filename=output_dir + self.name + '_tetrahedron_centers.csv', var=self.tetrahedron_centres)
        save_txt(filename=output_dir + self.name + '_edges.csv', var=self.edges)

    def save_to_ensight(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        save_ensight_geometry(directory=output_dir, name=self.name, nodes_xyz=self.nodes_xyz,
                              tetrahedrons=self.tetrahedrons)

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        if os.path.exists(input_dir + self.name + '_xyz.csv'):
            self.nodes_xyz = load_txt(filename=input_dir + self.name + '_xyz.csv')
            self.tetrahedrons = load_txt(filename=input_dir + self.name + '_tetra.csv').astype(int)
            if np.amin(self.tetrahedrons == 1):
                self.tetrahedrons = self.tetrahedrons - 1
            self.number_of_elements = self.tetrahedrons.shape[0]
            self.number_of_nodes = self.nodes_xyz.shape[0]
            self.triangles = load_txt(filename=input_dir + self.name + '_triangles.csv')
            self.number_of_triangles = self.triangles.shape[0]
            self.tetrahedron_centres = load_txt(filename=input_dir + self.name + '_tetrahedron_centres.csv')
            self.edges = load_txt(filename=input_dir + self.name + '_edges.csv')


class Fields:
    def __init__(self, name, field_type, verbose):
        self.name = name
        self.field_type = field_type
        self.dict = {}  # Placeholder.
        self.number_of_fields = len(self.dict)
        self.verbose = verbose

    def add_field(self, data, data_name, field_type):
        assert field_type == self.field_type, 'Input field type: ' + field_type + 'does not matching existing field ' \
                                                                                  'type: ' + field_type
        self.dict[data_name] = data

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        list_fields = list(self.dict.keys())
        self.number_of_fields = len(list_fields)
        for field_i in range(len(list_fields)):
            varname = list_fields[field_i]
            save_txt(filename=output_dir + self.name + '_' + self.field_type + '_' + varname + '.csv',
                     var=self.dict[varname])

    def save_to_ensight(self, output_dir, casename, geometry):
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        if not os.path.exists(output_dir + geometry.name + '.ensi.geo'):
            save_ensight_geometry(directory=output_dir, name=self.name, nodes_xyz=geometry.nodes_xyz,
                                  tetrahedrons=geometry.tetrahedrons)
        list_fields = list(self.dict.keys())
        field_dimensions = []
        list_fields_output = []
        for field_i in range(self.number_of_fields):
            varname = list_fields[field_i]
            if self.field_type == 'nodefield' or self.field_type == 'boundarynodefield' or \
                    self.field_type == 'postnodefield':
                if self.dict[varname].shape[0] == geometry.number_of_nodes:
                    save_ensight_node(directory=output_dir, name=self.name, field_name=varname, var=self.dict[varname])
                    list_fields_output.append(varname)
                    if len(self.dict[varname].shape) == 1:
                        field_dimensions.append(1)
                    elif self.dict[varname].shape[1] < 10:
                        field_dimensions.append(self.dict[varname].shape[1])
            elif self.field_type == 'elementfield' or self.field_type == 'material' or \
                    self.field_type == 'boundaryelementfield' or self.field_type == 'postelementfield':
                if self.dict[varname].shape[0] == geometry.number_of_elements:
                    save_ensight_element(directory=output_dir, name=self.name, field_name=varname,
                                         var=self.dict[varname])
                    list_fields_output.append(varname)
                    if len(self.dict[varname].shape) == 1:
                        field_dimensions.append(1)
                    elif self.dict[varname].shape[1] < 10:
                        field_dimensions.append(self.dict[varname].shape[1])
        save_ensight_case(directory=output_dir, name=casename, geometry_name=geometry.name,
                          field_names=list_fields_output,
                          field_dimensions=field_dimensions, field_type=self.field_type)

    def read_csv_to_attributes(self, input_dir, field_type):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        filenames = np.array([f for f in os.listdir(input_dir) if
                              os.path.isfile(os.path.join(input_dir, f)) and '_' + field_type + '_' in f])
        if filenames.shape[0] > 0:
            threadsNum = multiprocessing.cpu_count()
            self.dict = multiprocessing.Manager().dict()
            with pymp.Parallel(min(threadsNum, filenames.shape[0])) as p1:
                for file_i in p1.range(filenames.shape[0]):
            # if True:
            #     for file_i in range(filenames.shape[0]):
                    if self.verbose:
                        print('Reading in ' + input_dir + filenames[file_i])
                    varname = get_varname(filename=filenames[file_i], key=field_type)
                    self.dict[varname] = load_txt(filename=input_dir + filenames[file_i])
        self.number_of_fields = len(self.dict)


def save_txt(filename, var):
    if var is not None:
        np.savetxt(filename, var, delimiter=',')


def load_txt(filename):
    if os.path.exists(filename):
        return np.loadtxt(filename, delimiter=',')
    else:
        return None


def save_ensight_node(directory, name, field_name, var):
    if len(var.shape) == 1:
        with open(directory + name + '.ensi.' + field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Scalar per-node variables file\npart\n         1\ncoordinates\n')
            for node_i in range(0, var.shape[0]):
                f.write(str(var[node_i]) + '\n')
    else:
        if var.shape[1] == 1:
            with open(directory + name + '.ensi.' + field_name, 'w') as f:
                f.write('Alya Ensight Gold --- Scalar per-node variables file\npart\n         1\ncoordinates\n')
                for node_i in range(0, var.shape[0]):
                    f.write(str(var[node_i]) + '\n')
        else:
            with open(directory + name + '.ensi.' + field_name, 'w') as f:
                f.write('Alya Ensight Gold --- Vector per-node variables file\npart\n         1\ncoordinates\n')
                for c in range(var.shape[1]):
                    for i in range(var.shape[0]):
                        f.write(str(var[i, c]) + '\n')


def save_ensight_element(directory, name, field_name, var):
    if len(var.shape) == 1:
        with open(directory + name + '.ensi.' + field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Scalar per-cell variables file\npart\n         1\ntetra4\n')
            for elem_i in range(0, var.shape[0]):
                f.write(str(var[elem_i]) + '\n')
    else:
        if var.shape[1] == 1:
            with open(directory + name + '.ensi.' + field_name, 'w') as f:
                f.write('Alya Ensight Gold --- Scalar per-cell variables file\npart\n         1\ntetra4\n')
                for node_i in range(0, var.shape[0]):
                    f.write(str(var[node_i]) + '\n')
        with open(directory + name + '.ensi.' + field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Vector per-cell variables file\npart\n         1\ntetra4\n')
            for c in range(var.shape[1]):
                for elem_i in range(var.shape[0]):
                    f.write(str(var[elem_i, c]) + '\n')


def save_ensight_geometry(directory, name, nodes_xyz, tetrahedrons):
    if np.amin(tetrahedrons) == 0:
        tetrahedrons = tetrahedrons + 1  # Ensight takes node indices starting from 1.
    with open(directory + name + '.ensi.geo', 'w') as f:
        f.write(
            'Problem name:  ' + str(
                name) + '\nGeometry file\nnode id given\nelement id given\npart\n\t1\nVolume Mesh\ncoordinates\n' + str(
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
            f.write(
                str(tetrahedrons[i, 0]) + '\t' + str(tetrahedrons[i, 1]) + '\t' + str(tetrahedrons[i, 2]) + '\t' + str(
                    tetrahedrons[i, 3]) + '\n')


def save_ensight_case(directory, name, geometry_name, field_names, field_dimensions, field_type):
    with open(directory + name + '.ensi.case', 'w') as f:
        f.write(
            '#\n# Alya generated post-process files\n# Ensight Gold Format\n#\n# Problem name:\t' + str(name) + '\n#\n')
        f.write('FORMAT\ntype:\tensight gold\nGEOMETRY\nmodel:\t1\t' + geometry_name + '.ensi.geo\nVARIABLE\n')
        for field_i in range(len(field_names)):
            field_name = field_names[field_i]
            field_dimension = field_dimensions[field_i]
            if field_type == 'nodefield' or field_type == 'boundarynodefield' or field_type == 'postnodefield':
                ensight_field_type = 'node'
            elif field_type == 'elementfield' or field_type == 'boundaryelementfield' or field_type == 'postelementfield':
                ensight_field_type = 'element'
            elif field_type == 'material':
                ensight_field_type = 'element'
            else:
                raise ValueError('Field type '+ field_type + ' not found. ')
            if field_dimension == 1:
                f.write(
                    'scalar per ' + ensight_field_type + ':\t' + str(
                        field_dimension) + '\t' + field_name + '\t' + geometry_name + '.ensi.' + field_name + '\n')
            elif field_dimension > 1:
                f.write(
                    'vector per ' + ensight_field_type + ':\t' + str(
                        field_dimension) + '\t' + field_name + '\t' + geometry_name + '.ensi.' + field_name + '\n')


def get_varname(filename, key):
    return filename.split('_' + key + '_')[-1].split('.csv')[0]


if __name__ == '__main__':
    wd = os.getcwd()
    geometry = Geometry('test', verbose=False)
    geometry.save_to_csv(wd)
    node_fields = Fields('test', field_type='nodefield', verbose=False)
    node_fields.save_to_csv(wd)
    element_fields = Fields('test', field_type='elementfield', verbose=False)
    element_fields.save_to_csv(wd)
