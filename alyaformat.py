from myformat import *
import json
from shutil import copy


class AlyaFormat:
    def __init__(self, name, geometry_and_fields_input_dir, simulation_json_file):
        self.version = 'alya-compbiomed2' # 'Alya_multiple_BZRZ_models'
        self.template_dir = 'alya_input_templates/'
        self.name = name
        self.output_dir = simulation_json_file.split('.')[0] + '/'
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        copy(simulation_json_file, self.output_dir)
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.geometry_and_fields_input_dir = geometry_and_fields_input_dir
        self.geometry = None
        self.node_fields = None
        self.element_fields = None
        self.materials = None
        self.import_geometry_and_fields()
        self.section_divider = '$------------------------------------------------------\n'
        self.write_dat()
        self.write_dom_dat()
        # if 'EXMEDI' in self.simulation_dict.physics:
        #     self.write_exm_dat()
        # if 'SOLIDZ' in self.simulation_dict.physics:
        #     self.write_sld_dat()
        # self.write_ker_dat()
        # self.write_post_dat()
        # self.write_cell_txt()
        # self.write_job_scripts()

    def import_geometry_and_fields(self):
        # Read geometry from geometry directory CSV files
        print('Read in geometry and fields from CSV files.')
        self.geometry = Geometry(self.name)
        self.geometry.read_csv_to_attributes(self.geometry_and_fields_input_dir)
        self.node_fields = Fields(self.name, field_type='nodefield')
        self.node_fields.read_csv_to_attributes(self.geometry_and_fields_input_dir)
        self.element_fields = Fields(self.name, field_type='elementfield')
        self.element_fields.read_csv_to_attributes(self.geometry_and_fields_input_dir)
        self.materials = Materials(self.name)
        self.materials.read_csv_to_attributes(self.geometry_and_fields_input_dir)

    def generate_fields(self):
        # Generate required fields for Alya simulations.
        print('Make use of field_evaluation_functions to generate Alya fields...')

    def write_dat(self):
        filename = self.output_dir+self.name+'.dat'
        print('Writing out '+filename)
        if self.version == 'alya-compbiomed2' or self.version == 'Alya_multiple_BZRZ_models':
            with open(self.template_dir + self.version + '.dat', 'r') as f:
                data = f.readlines()
            data = ''.join(data)
            data = data.replace('<<name>>', self.name)
            if self.simulation_dict['run_type'] == 'CONTINUE':
                run_type_str = '\tRUNTYPE: CONTINUE\n'
            elif self.simulation_dict['run_type'] == 'PRELIMINARY':
                run_type_str = '\tRUNTYPE: PRELIMINARY, FREQUENCY = ' + str(
                    self.simulation_dict['warm_start_write_out_frequency']) + '\n'
            elif self.simulation_dict['run_type'] == 'PRELIMINARY CONTINUE':
                run_type_str = '\tRUNTYPE: CONTINUE\n\tRUNTYPE: PRELIMINARY, FREQUENCY = ' + str(
                    self.simulation_dict['warm_start_write_out_frequency']) + '\n'
            else:
                raise 'write_dat: Run type ' + self.simulation_dict['run_type'] + ' not found.'
            data = data.replace('<<run_type_str>>', run_type_str)
            data = data.replace('<<time_start>>', str(self.simulation_dict['time_start']))
            data = data.replace('<<time_end>>', str(self.simulation_dict['time_end']))

            if self.simulation_dict['time_function'] is not None:
                time_function_str = ''
                time_function_str = time_function_str + '\tFUNCTION: DISCRETE\n'
                time_function_str = time_function_str + '\tDTIME_FUNCTION:\n'
                time_function_str = time_function_str + '\t\tTIME_SHAPE: DISCRETE\n'
                time_function_str = time_function_str + '\t\tSHAPE_DEFINITION:\n'
                time_function = np.array(self.simulation_dict['time_function'])
                time_function_str = time_function_str + '\t\t\t' + str(time_function.shape[0]) + '\n'
                for shape_i in range(time_function.shape[0]):
                    time_function_str = time_function_str + '\t\t\t' + str(time_function[shape_i, 0]) + '\t' + str(
                        time_function[shape_i, 1]) + '\n'
                time_function_str = time_function_str + '\t\tEND_SHAPE_DEFINITION\n'
                time_function_str = time_function_str + '\tEND_DTIME_FUNCTION\n'
                data = data.replace('<<time_function_str>>', time_function_str)
            else:
                data = data.replace('<<time_function_str>>', '')

            physics_str = ''
            for physics_i in range(len(self.simulation_dict['physics'])):
                if self.simulation_dict['physics'][physics_i] == 'SOLIDZ':
                    physics_str = physics_str + '\tSOLIDZ_PROBLEM:  ON\n'
                    physics_str = physics_str + '\tEND_SOLIDZ'
                if self.simulation_dict['physics'][physics_i] == 'EXMEDI':
                    physics_str = physics_str + '\tEXMEDI_PROBLEM:  ON\n'
                    physics_str = physics_str + '\tEND_EXMEDI'
            data = data.replace('<<physics_str>>', physics_str)
            with open(self.output_dir+self.name+'.dat', 'w') as f:
                f.write(data)

    def write_dom_dat(self):
        if self.geometry is not None:
            with open(self.template_dir + self.version + '.dom.dat', 'r') as f:
                data = f.readlines()
            data = ''.join(data)
            data = data.replace('<<number_of_nodes>>', str(self.geometry.number_of_nodes))
            data = data.replace('<<number_of_elements>>', str(self.geometry.number_of_elements))
            data = data.replace('<<spatial_dimensions>>', str(self.simulation_dict['spatial_dimensions']))
            data = data.replace('<<number_of_nodes_per_element>>', str(self.geometry.tetrahedrons[0].shape[0]))
            data = data.replace('<<number_of_boundaries>>', str(self.geometry.number_of_triangles))
            data = data.replace('<<number_of_materials>>', str(np.amax(self.materials)))
            data = data.replace('<<number_of_fields>>', str(self.node_fields.number_of_fields +
                                                            self.element_fields.number_of_fields))
            field_declaration_str = ''
            varnames = list(self.node_fields.dict.keys())
            print(varnames)
            for field_i in range(self.node_fields.number_of_fields):
                varname = varnames[field_i]
                assert self.node_fields.dict[varname].shape[0] == self.geometry.number_of_nodes, \
                    'Field: ' + varname + ' does not match in first dimension with the number of nodes ' \
                                          'in the geometry'
                field_declaration_str = field_declaration_str + '\t\tFIELD = ' + str(field_i) + ', DIMENSION = ' + str(
                    self.node_fields.dict[varname].shape[1]) + ', NODES\n'
            for field_i in range(self.element_fields.number_of_fields):
                varname = varnames[field_i]
                assert self.element_fields.dict[varname].shape[0] == self.geometry.number_of_elements, \
                    'Field: ' + varname + ' does not match in first dimension with the number of elements ' \
                                          'in the geometry'
                field_declaration_str = field_declaration_str + '\t\tFIELD = ' + str(field_i) + ', DIMENSION = ' + str(
                    self.element_fields.dict[varname].shape[1]) + ', ELEMENTS\n'
            data = data.replace('<<field_declaration_str>>', field_declaration_str)
            data = data.replace('<<x_scale>>', str(self.simulation_dict['x_scale']))
            data = data.replace('<<y_scale>>', str(self.simulation_dict['y_scale']))
            data = data.replace('<<z_scale>>', str(self.simulation_dict['z_scale']))
            data = data.replace('<<x_translation>>', str(self.simulation_dict['x_translation']))
            data = data.replace('<<y_translation>>', str(self.simulation_dict['y_translation']))
            data = data.replace('<<z_translation>>', str(self.simulation_dict['z_translation']))
            data = data.replace('<<geometry_file_name>>', self.name+'.geo')
            data = data.replace('<<surface_file_name>>', self.name + '.surface')
            data = data.replace('<<materials_file_name>>', self.name + '.materials')
            data = data.replace('<<sets_file_name>>', self.name+'.sets')
            data = data.replace('<<boundaries_file_name>>', self.name+'.boundaries')
            field_initialisation_str = ''
            for field_i in range(self.element_fields.number_of_fields):
                varname = varnames[field_i]
                field_initialisation_str = field_initialisation_str + '\tFIELD = ' + str(field_i) + '\n'
                field_initialisation_str = field_initialisation_str + '\t\tINCLUDE ' + str(self.name) + '.' + varname + '\n'
                field_initialisation_str = field_initialisation_str + '\tEND_FIELD\n'
            data = data.replace('<<field_initialisation_str>>', field_initialisation_str)
            with open(self.output_dir + self.name + '.dom.dat', 'w') as f:
                f.write(data)
        else:
            raise ValueError('Please read in geometry information first before writing .dom.dat file. ')


    def write_exm_dat(self):
        if self.version == 'alya-compbiomed2':
            print('adfaf')



    def write_sld_dat(self):
        print('fdaf')


    def write_ker_dat(self):
        pass

    def write_post_dat(self):
        pass

    def write_cell_txt(self):
        pass

    def write_job_scripts(self):
        pass

