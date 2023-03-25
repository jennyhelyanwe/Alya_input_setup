import json
from shutil import copy

from meshstructure import MeshStructure
from pre_refactoring.evaluationfunctions import *


class AlyaFormat(MeshStructure):
    def __init__(self, name, geometric_data_dir, personalisation_dir, simulation_json_file, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.version = 'alya-compbiomed2'  # 'Alya_multiple_BZRZ_models'
        self.template_dir = 'alya_input_templates/'
        self.name = name
        print('Alya version: '+self.version + ', simulation name: ', self.name)
        self.output_dir = simulation_json_file.split('.')[0] + '_' + self.name + '/'
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        copy(simulation_json_file, self.output_dir + self.name + '.json')
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.geometric_data_dir = geometric_data_dir
        self.personalisation_dir = personalisation_dir
        self.verbose = verbose
        self.write_alya_simulation_files()

    def write_alya_simulation_files(self):
        self.write_dat()
        self.write_dom_dat()
        if 'EXMEDI' in self.simulation_dict['physics']:
            self.write_exm_dat()
        # if 'SOLIDZ' in self.simulation_dict['physics']:
        #     self.write_sld_dat()
        # self.write_ker_dat()
        # self.write_post_dat()
        # self.write_cell_txt()
        # self.write_job_scripts()

    def write_dat(self):
        filename = self.output_dir + self.name + '.dat'
        print('Writing out ' + filename)
        if self.version == 'alya-compbiomed2' or self.version == 'Alya_multiple_BZRZ_models':
            with open(self.template_dir + self.version + '.dat', 'r') as f:
                data = f.readlines()
            data = ''.join(data)
            data = data.replace('<<name>>', self.name)
            if self.simulation_dict['run_type'] == 'CONTINUE':
                run_type_str = 'RUNTYPE: CONTINUE\n'
            elif self.simulation_dict['run_type'] == 'PRELIMINARY':
                run_type_str = 'RUNTYPE: PRELIMINARY, FREQUENCY = ' + str(
                    self.simulation_dict['warm_start_write_out_frequency']) + '\n'
            elif self.simulation_dict['run_type'] == 'PRELIMINARY CONTINUE':
                run_type_str = 'RUNTYPE: CONTINUE\nRUNTYPE: PRELIMINARY, FREQUENCY = ' + str(
                    self.simulation_dict['warm_start_write_out_frequency']) + '\n'
            else:
                raise 'write_dat: Run type ' + self.simulation_dict['run_type'] + ' not found.'
            data = data.replace('<<run_type_str>>', run_type_str)
            data = data.replace('<<time_start>>', str(self.simulation_dict['time_start']))
            data = data.replace('<<time_end>>', str(self.simulation_dict['time_end']))

            if self.simulation_dict['time_function'] is not None:
                time_function_str = ''
                time_function_str = time_function_str + 'FUNCTION: DISCRETE\n'
                time_function_str = time_function_str + '\tDTIME_FUNCTION:\n'
                time_function_str = time_function_str + '\tTIME_SHAPE: DISCRETE\n'
                time_function_str = time_function_str + '\tSHAPE_DEFINITION:\n'
                time_function = np.array(self.simulation_dict['time_function'])
                time_function_str = time_function_str + '\t\t' + str(time_function.shape[0]) + '\n'
                for shape_i in range(time_function.shape[0]):
                    time_function_str = time_function_str + '\t\t' + str(time_function[shape_i, 0]) + '\t' + str(
                        time_function[shape_i, 1]) + '\n'
                time_function_str = time_function_str + '\tEND_SHAPE_DEFINITION\n'
                time_function_str = time_function_str + '\tEND_DTIME_FUNCTION\n'
                data = data.replace('<<time_function_str>>', time_function_str)
            else:
                data = data.replace('<<time_function_str>>', '')

            physics_str = ''
            for physics_i in range(len(self.simulation_dict['physics'])):
                if self.simulation_dict['physics'][physics_i] == 'SOLIDZ':
                    physics_str = physics_str + 'SOLIDZ_PROBLEM:  ON\n'
                    physics_str = physics_str + 'END_SOLIDZ'
                if self.simulation_dict['physics'][physics_i] == 'EXMEDI':
                    physics_str = physics_str + 'EXMEDI_PROBLEM:  ON\n'
                    physics_str = physics_str + 'END_EXMEDI'
            data = data.replace('<<physics_str>>', physics_str)
            with open(self.output_dir + self.name + '.dat', 'w') as f:
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
            data = data.replace('<<number_of_materials>>', str(np.amax(self.materials.dict['tetra']).astype(int)))
            data = data.replace('<<number_of_fields>>', str(int(len(self.simulation_dict['field_names']))))
            field_declaration_str = ''
            field_names = self.simulation_dict['field_names']
            field_types = self.simulation_dict['field_types']
            for field_i in range(len(field_names)):
                varname = field_names[field_i]
                if field_types[field_i] == 'nodefield':
                    field = self.node_fields.dict[varname]
                    if len(field.shape) == 1:
                        field_declaration_str = field_declaration_str + '\n\t\tFIELD = ' + str(
                            field_i + 1) + ', DIMENSION = 1, NODES'
                    else:
                        field_declaration_str = field_declaration_str + '\n\t\tFIELD = ' + str(
                            field_i + 1) + ', DIMENSION = ' + str(
                            field.shape[1]) + ', NODES'
                elif field_types[field_i] == 'elementfield':
                    field = self.node_fields.dict[varname]
                    if len(field.shape) == 1:
                        field_declaration_str = field_declaration_str + '\n\t\tFIELD = ' + str(
                            field_i) + ', DIMENSION = 1, ELEMENTS'
                    else:
                        field_declaration_str = field_declaration_str + '\n\t\tFIELD = ' + str(
                            field_i) + ', DIMENSION = ' + str(
                            field.shape[1]) + ', ELEMENTS'
            data = data.replace('<<field_declaration_str>>', field_declaration_str)
            data = data.replace('<<x_scale>>', str(self.simulation_dict['x_scale']))
            data = data.replace('<<y_scale>>', str(self.simulation_dict['y_scale']))
            data = data.replace('<<z_scale>>', str(self.simulation_dict['z_scale']))
            data = data.replace('<<x_translation>>', str(self.simulation_dict['x_translation']))
            data = data.replace('<<y_translation>>', str(self.simulation_dict['y_translation']))
            data = data.replace('<<z_translation>>', str(self.simulation_dict['z_translation']))
            data = data.replace('<<geometry_tetrahedron_file_name>>', self.name + '.tetra')
            data = data.replace('<<geometry_nodes_xyz_file_name>>', self.name + '.nodes_xyz')
            data = data.replace('<<surface_file_name>>', self.name + '.surface')
            data = data.replace('<<materials_file_name>>', self.name + '.materials')
            data = data.replace('<<sets_file_name>>', self.name + '.sets')
            data = data.replace('<<boundaries_file_name>>', self.name + '.boundaries')
            field_initialisation_str = ''
            for field_i in range(len(field_names)):
                varname = field_names[field_i]
                field_initialisation_str = field_initialisation_str + '\n\tFIELD = ' + str(field_i) + '\n'
                field_initialisation_str = field_initialisation_str + '\t\tINCLUDE ' + str(
                    self.name) + '.' + varname + '\n'
                field_initialisation_str = field_initialisation_str + '\tEND_FIELD'
            data = data.replace('<<field_initialisation_str>>', field_initialisation_str)
            filename = self.output_dir + self.name + '.dom.dat'
            print('Writing out ' + filename)
            with open(self.output_dir + self.name + '.dom.dat', 'w') as f:
                f.write(data)

            # Write out required geometric fields in Alya format
            print('Writing out geometric files in Alya format')
            write_alya_field(filename=self.output_dir + self.name + '.nodes_xyz',
                             field_idx=np.arange(1, self.geometry.number_of_nodes),
                             field_data=self.geometry.nodes_xyz)
            write_alya_field(filename=self.output_dir + self.name + '.tetra',
                             field_idx=np.arange(1, self.geometry.number_of_nodes),
                             field_data=self.geometry.tetrahedrons.astype(int) + 1)
            write_alya_field(filename=self.output_dir + self.name + '.surfaces',
                             field_idx=np.arange(1, self.geometry.number_of_triangles),
                             field_data=self.geometry.triangles.astype(int) + 1)
            write_alya_field(filename=self.output_dir + self.name + '.materials',
                             field_idx=np.arange(1, self.materials.dict['tetra'].shape[0]),
                             field_data=self.materials.dict['tetra'])
            write_alya_field(filename=self.output_dir + self.name + '.boundaries',
                             field_idx=np.arange(1, self.boundary_element_fields.dict[
                                 'mechanical-element-boundary-label'].shape[0]),
                             field_data=self.boundary_element_fields.dict['mechanical-element-boundary-label'].astype(
                                 int) + 1)
            write_alya_field(filename=self.output_dir + self.name + '.sets', field_idx=np.arange(1,
                                                                                                 self.boundary_element_fields.dict[
                                                                                                     'mechanical-element-boundary-label'].shape[
                                                                                                     0]),
                             field_data=self.boundary_element_fields.dict['mechanical-element-boundary-label'].astype(
                                 int) + 1)

            # Write out required fields in Alya format
            print('Writing out field files in Alya format')
            field_names = self.simulation_dict['field_names']
            field_types = self.simulation_dict['field_types']
            for field_i in range(len(field_names)):
                varname = field_names[field_i]
                field = None
                if field_types[field_i] == 'nodefield':
                    field = self.node_fields.dict[varname]
                elif field_types[field_i] == 'elementfield':
                    field = self.element_fields.dict[varname]
                elif field_types[field_i] == 'material':
                    field = self.materials.dict[varname]
                if varname == 'endocardial-activation-times':
                    field_idx = self.boundary_node_fields.dict['endocardial-nodes'].astype(int) + 1
                else:
                    field_idx = np.arange(1, field.shape[0])
                write_alya_field(filename=self.output_dir + self.name + '.' + varname, field_idx=field_idx,
                                 field_data=field)
        else:
            raise ValueError('Please read in geometry information first before writing .dom.dat file. ')

    def write_exm_dat(self):
        if self.version == 'alya-compbiomed2':
            with open(self.template_dir + self.version + '.exm.dat', 'r') as f:
                data = f.readlines()
            data = ''.join(data)
            stimulus_field_number = self.simulation_dict['field_names'].index(
                self.simulation_dict['stimulus_field_name']) + 1
            data = data.replace('<<stimulus_field_number>>', str(stimulus_field_number))
            data = data.replace('<<monodomain_delay>>', str(self.simulation_dict['prestress_time']))
            num_electrodes = self.node_fields.dict['electrode_xyz'].shape[0]
            data = data.replace('<<number_of_ecg_electrodes>>', str(num_electrodes))
            ecg_electrode_coordinates_str = ''
            for i in range(num_electrodes):
                ecg_electrode_coordinates_str = ecg_electrode_coordinates_str + '\n\t\t' + \
                                                str(self.node_fields.dict['electrode_xyz'][i, 0]) + ' ' + \
                                                str(self.node_fields.dict['electrode_xyz'][i, 1]) + ' ' + \
                                                str(self.node_fields.dict['electrode_xyz'][i, 2])
            data = data.replace('<<ecg_electrodes_coordinates>>', ecg_electrode_coordinates_str)
            with open(self.template_dir + self.version + '.subtemplate.exmedi_property_template', 'r') as f:
                subdata = f.readlines()
            subdata = ''.join(subdata)
            exmedi_properties_str = ''
            for material_i in range(np.amax(self.materials.dict['tetra'].astype(int))):
                str_i = '\n'
                str_i = str_i + subdata.replace('<<material_idx>>', str(material_i + 1))
                str_i = str_i.replace('<<diffusivities>>',
                                      str(self.simulation_dict['sigma'][material_i][0]) + ','
                                      + str(self.simulation_dict['sigma'][material_i][1]) + ','
                                      + str(self.simulation_dict['sigma'][material_i][2]))
                str_i = str_i.replace('<<cell_model>>', self.simulation_dict['cell_model'][material_i])
                str_i = str_i.replace('<<cell_initialisation_txt_file_name>>',
                                      self.simulation_dict['cell_filename'][material_i])
                exmedi_properties_str = exmedi_properties_str + str_i
            data = data.replace('<<exmedi_properties_str>>', exmedi_properties_str)
            with open(self.template_dir + self.version + '.subtemplate.exmedi_property_template', 'r') as f:
                subdata = f.readlines()
            subdata = ''.join(subdata)
            exmedi_scaling_factor_str = ''
            sf_fields = [x for x in list(self.simulation_dict.keys()) if 'sf_' in x]
            for sf_i in range(len(sf_fields)):
                str_i = '\n'
                str_i = str_i + subdata.replace('<<scaling_name>>', sf_fields[sf_i])
                str_i = str_i.replace('<<scaling_name_field_number>>',
                                      str(list(self.simulation_dict.keys()).index(sf_fields[sf_i])+1))
                exmedi_scaling_factor_str = exmedi_scaling_factor_str + str_i
            data = data.replace('<<exmedi_scaling_factor_str>>', exmedi_scaling_factor_str)
            with open(self.template_dir + self.version + '.subtemplate.postprocess_template', 'r') as f:
                subdata = f.readlines()
            subdata = ''.join(subdata)
            post_process_exmedi_str = ''
            for var_i in range(len(self.simulation_dict['exmedi_output'])):
                str_i = '\n'
                str_i = str_i + subdata.replace('<<post_process_name>>', self.simulation_dict['exmedi_output'][var_i])
                str_i = str_i.replace('<<post_process_period>>', str(self.simulation_dict['exmedi_output_period'][var_i]))
                post_process_exmedi_str = post_process_exmedi_str + str_i
            data = data.replace('<<post_process_exmedi_str>>', post_process_exmedi_str)
            # Write out Exmedi file
            filename = self.output_dir + self.name + '.exm.dat'
            print('Writing out ' + filename)
            with open(filename, 'w') as f:
                f.write(data)

    def write_sld_dat(self):
        if self.version == 'alya-compbiomed2':
            with open(self.template_dir + self.version + '.sld.dat', 'r') as f:
                data = f.readlines()
            data = ''.join(data)
            data = data.replace()

    def write_ker_dat(self):
        pass

    def write_post_dat(self):
        pass

    def write_cell_txt(self):
        pass

    def write_job_scripts(self):
        pass


def write_alya_field(filename, field_idx, field_data):
    field_idx = np.array(field_idx).astype(int)
    field_data = np.array(field_data)
    with open(filename, 'w') as f:
        for i in range(len(field_idx)):
            f.write(str(field_idx[i]))
            if len(field_data.shape) == 1:
                f.write('\t' + str(field_data[i]))
            else:
                for j in range(field_data.shape[1]):
                    f.write('\t' + str(field_data[i, j]))
            f.write('\n')
