from myformat import *


class MeshStructure:
    def __init__(self, name, geometric_data_dir, max_cores_used, verbose):
        if geometric_data_dir[-1] != '/':
            geometric_data_dir = geometric_data_dir + '/'
        if not os.path.exists(geometric_data_dir + 'ensight/'):
            os.makedirs(geometric_data_dir + 'ensight/')
        print('Reading in mesh structure ' + name + ' from ' + geometric_data_dir)
        self.geometric_data_dir = geometric_data_dir
        self.name = name
        self.geometry = Geometry(self.name, max_cores_used=max_cores_used, verbose=verbose)
        self.node_fields = Fields(self.name, field_type='nodefield', max_cores_used=max_cores_used, verbose=verbose)
        self.boundary_node_fields = Fields(self.name, field_type='boundarynodefield', max_cores_used=max_cores_used,verbose=verbose)
        self.element_fields = Fields(self.name, field_type='elementfield', max_cores_used=max_cores_used,verbose=verbose)
        self.boundary_element_fields = Fields(self.name, field_type='boundaryelementfield', max_cores_used=max_cores_used,verbose=verbose)
        self.materials = Fields(self.name, field_type='material', max_cores_used=max_cores_used, verbose=verbose)
        # Read in existing fields
        self.geometry.read_csv_to_attributes(self.geometric_data_dir)
        self.node_fields.read_csv_to_attributes(self.geometric_data_dir, field_type='nodefield')
        self.element_fields.read_csv_to_attributes(self.geometric_data_dir, field_type='elementfield')
        self.boundary_node_fields.read_csv_to_attributes(self.geometric_data_dir, field_type='boundarynodefield')
        self.boundary_element_fields.read_csv_to_attributes(self.geometric_data_dir, field_type='boundaryelementfield')
        self.materials.read_csv_to_attributes(self.geometric_data_dir, field_type='material')
        self.max_cores_used = max_cores_used
        print('done reading mesh structure. ')


    def save(self):
        # Save to CSV
        print('Saving geometry and fields to CSV')
        self.geometry.save_to_csv(self.geometric_data_dir)
        self.node_fields.save_to_csv(self.geometric_data_dir)
        self.element_fields.save_to_csv(self.geometric_data_dir)
        self.materials.save_to_csv(self.geometric_data_dir)
        self.boundary_node_fields.save_to_csv(self.geometric_data_dir)
        self.boundary_element_fields.save_to_csv(self.geometric_data_dir)
        # Save to ensight
        print('Saving geometry and fields to Ensight')
        self.geometry.save_to_ensight(self.geometric_data_dir + 'ensight/')
        self.node_fields.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                         casename=self.name + '_nodefield', geometry=self.geometry)
        self.element_fields.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                            casename=self.name + '_elementfield', geometry=self.geometry)
        self.materials.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                       casename=self.name + '_material', geometry=self.geometry)
        self.boundary_node_fields.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                                  casename=self.name + '_boundary_nodefield',
                                                  geometry=self.geometry)
        self.boundary_element_fields.save_to_ensight(output_dir=self.geometric_data_dir + 'ensight/',
                                                     casename=self.name + '_boundary_elementfield',
                                                     geometry=self.geometry)
