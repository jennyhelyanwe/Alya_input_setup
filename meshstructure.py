from myformat import *
class MeshStructure:
    def __init__(self, name,  geometric_data_dir, boundary_data_dir, field_data_dir, verbose):
        if geometric_data_dir[-1] != '/':
            geometric_data_dir = geometric_data_dir + '/'
        if boundary_data_dir[-1] != '/':
            boundary_data_dir = boundary_data_dir + '/'
        if field_data_dir[-1] != '/':
            field_data_dir = field_data_dir + '/'
        if not os.path.exists(geometric_data_dir + 'ensight/'):
            os.makedirs(geometric_data_dir + 'ensight/')
        if not os.path.exists(boundary_data_dir + 'ensight/'):
            os.makedirs(boundary_data_dir+ 'ensight/')
        if not os.path.exists(field_data_dir + 'ensight/'):
            os.makedirs(field_data_dir + 'ensight/')
        self.geometric_data_dir = geometric_data_dir
        self.boundary_data_dir = boundary_data_dir
        self.field_data_dir = field_data_dir
        self.name = name
        self.geometry = Geometry(self.name, verbose=verbose)
        self.node_fields = Fields(self.name, field_type='nodefield', verbose=verbose)
        self.boundary_node_fields = Fields(self.name, field_type='nodefield', verbose=verbose)
        self.element_fields = Fields(self.name, field_type='elementfield', verbose=verbose)
        self.boundary_element_fields = Fields(self.name, field_type='elementfield', verbose=verbose)
        self.materials = Fields(self.name, field_type='material', verbose=verbose)
        # Read in existing fields
        self.geometry.read_csv_to_attributes(self.geometric_data_dir)
        self.node_fields.read_csv_to_attributes(self.geometric_data_dir, field_type='nodefield')
        self.node_fields.read_csv_to_attributes(self.field_data_dir, field_type='nodefield')
        self.element_fields.read_csv_to_attributes(self.geometric_data_dir, field_type='elementfield')
        self.element_fields.read_csv_to_attributes(self.field_data_dir, field_type='elementfield')
        self.boundary_node_fields.read_csv_to_attributes(self.boundary_data_dir, field_type='nodefield')
        self.boundary_element_fields.read_csv_to_attributes(self.boundary_data_dir, field_type='elementfield')
        self.materials.read_csv_to_attributes(self.geometric_data_dir, field_type='material')

    def save(self):
        self.geometry.save_to_csv(self.geometric_data_dir)
        self.node_fields.save_to_csv(self.geometric_data_dir)
        self.element_fields.save_to_csv(self.geometric_data_dir)
        self.materials.save_to_csv(self.geometric_data_dir)
        self.geometry.save_to_ensight(self.geometric_data_dir + 'ensight/')
        self.node_fields.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                         casename=self.name + '_nodefield', geometry=self.geometry)
        self.element_fields.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                            casename=self.name + '_elementfield', geometry=self.geometry)
        self.materials.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                       casename=self.name + '_material', geometry=self.geometry)
        self.boundary_node_fields.save_to_csv(self.boundary_data_dir)
        self.boundary_element_fields.save_to_csv(self.boundary_data_dir)
        self.geometry.save_to_ensight(self.boundary_data_dir + 'ensight/')
        self.boundary_node_fields.save_to_ensight(output_dir=self.boundary_data_dir + 'ensight/',
                                                  casename=self.name + '_boundary_nodefield',
                                                  geometry=self.geometry)
        self.boundary_element_fields.save_to_ensight(output_dir=self.boundary_data_dir + 'ensight/',
                                                     casename=self.name + '_boundary_elementfield',
                                                     geometry=self.geometry)

