import json
from shutil import copy

from meshstructure import MeshStructure
from pre_refactoring.evaluationfunctions import *


class AlyaFormat(MeshStructure):
    def __init__(self, name, geometric_data_dir, personalisation_dir, simulation_json_file, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        self.version = 'alya-compbiomed2'  # 'Alya_multiple_BZRZ_models'
        self.template_dir = 'alya_input_templates/'
        self.job_template_dir = 'job_script_template/'
        self.job_version = 'jureca'
        print('Alya version: '+self.version + ', simulation name: ', self.name)
        self.output_dir = simulation_json_file.split('.')[0] + '_' + self.name + '/'
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        copy(simulation_json_file, self.output_dir + self.name + '.json')
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.name = self.simulation_dict['name']
        self.geometric_data_dir = geometric_data_dir
        self.personalisation_dir = personalisation_dir
        self.verbose = verbose
        self.write_alya_simulation_files()

    def template(self, filename, keys, data, num_duplicates):
        with open(filename, 'r') as f:
            subdata = f.readlines()
        subdata = ''.join(subdata)
        replacement_str = ''
        for duplicate_i in range(num_duplicates):
            if num_duplicates > 1:
                str_i = '\n' + subdata
            else:
                str_i = subdata
            for key_i in range(len(keys)):
                str_i = str_i.replace('<<' + keys[key_i] + '>>', str(data[duplicate_i][key_i]))
            replacement_str = replacement_str + str_i
        return replacement_str

    def write_alya_simulation_files(self):
        # self.write_dat()
        # self.write_dom_dat()
        # if 'EXMEDI' in self.simulation_dict['physics']:
        #     self.write_exm_dat()
        # if 'SOLIDZ' in self.simulation_dict['physics']:
        #     self.write_sld_dat()
        # self.write_ker_dat()
        # self.write_post_dat()
        self.write_cell_txt()
        # self.write_job_scripts()

    def write_dat(self):
        filename = self.output_dir + self.name + '.dat'
        print('Writing out ' + filename)
        if self.version == 'alya-compbiomed2' or self.version == 'Alya_multiple_BZRZ_models':
            if self.simulation_dict['run_type'] == 'CONTINUE':
                run_type_str = '\n\tRUNTYPE: CONTINUE'
            elif self.simulation_dict['run_type'] == 'PRELIMINARY':
                run_type_str = '\n\tRUNTYPE: PRELIMINARY, FREQUENCY = ' + str(
                    self.simulation_dict['warm_start_write_out_frequency'])
            elif self.simulation_dict['run_type'] == 'PRELIMINARY CONTINUE':
                run_type_str = '\n\tRUNTYPE: CONTINUE\n\tRUNTYPE: PRELIMINARY, FREQUENCY = ' + str(
                    self.simulation_dict['warm_start_write_out_frequency'])
            else:
                raise 'write_dat: Run type ' + self.simulation_dict['run_type'] + ' not found.'
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
            else:
                time_function_str = ''
            physics_str = ''
            for physics_i in range(len(self.simulation_dict['physics'])):
                if self.simulation_dict['physics'][physics_i] == 'SOLIDZ':
                    physics_str = physics_str + '\n\tSOLIDZ_PROBLEM:  ON\n'
                    physics_str = physics_str + '\tEND_SOLIDZ'
                if self.simulation_dict['physics'][physics_i] == 'EXMEDI':
                    physics_str = physics_str + '\n\tEXMEDI_PROBLEM:  ON\n'
                    physics_str = physics_str + '\tEND_EXMEDI'
            keys = ["name", "run_type_str", "time_start", "time_end", "time_function_str", "physics_str"]
            insert_data = [[self.name, run_type_str, self.simulation_dict['time_start'],
                           self.simulation_dict['time_end'],
                           time_function_str, physics_str]]
            filename = self.template_dir + self.version + '.dat'
            data = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
            with open(self.output_dir + self.name + '.dat', 'w') as f:
                f.write(data)

    def write_dom_dat(self):
        if self.geometry is not None:
            filename = self.template_dir + self.version + '.dom.dat'
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
                    field = self.element_fields.dict[varname]
                    if len(field.shape) == 1:
                        field_declaration_str = field_declaration_str + '\n\t\tFIELD = ' + str(
                            field_i + 1) + ', DIMENSION = 1, ELEMENTS'
                    else:
                        field_declaration_str = field_declaration_str + '\n\t\tFIELD = ' + str(
                            field_i + 1) + ', DIMENSION = ' + str(
                            field.shape[1]) + ', ELEMENTS'
            field_initialisation_str = ''
            for field_i in range(len(field_names)):
                varname = field_names[field_i]
                field_initialisation_str = field_initialisation_str + '\n\tFIELD = ' + str(field_i + 1) + '\n'
                field_initialisation_str = field_initialisation_str + '\t\tINCLUDE ' + str(
                    self.name) + '.' + varname + '\n'
                field_initialisation_str = field_initialisation_str + '\tEND_FIELD'
            keys = ["number_of_nodes", "number_of_elements", "spatial_dimensions", "number_of_nodes_per_element",
                    "number_of_boundaries", "number_of_materials", "number_of_fields", "field_declaration_str",
                    "x_scale", "y_scale", "z_scale", "x_translation", "y_translation", "z_translation",
                    "geometry_file_name", "surface_file_name",
                    "materials_file_name", "sets_file_name", "boundaries_file_name", "field_initialisation_str"]
            insert_data = [[self.geometry.number_of_nodes, self.geometry.number_of_elements,
                           self.simulation_dict['spatial_dimensions'], self.geometry.tetrahedrons[0].shape[0],
                           self.geometry.number_of_triangles, np.amax(self.materials.dict['tetra']).astype(int),
                           int(len(self.simulation_dict['field_names'])), field_declaration_str,
                           self.simulation_dict['x_scale'], self.simulation_dict['y_scale'],
                           self.simulation_dict['z_scale'], self.simulation_dict['x_translation'],
                           self.simulation_dict['y_translation'], self.simulation_dict['z_translation'],
                           self.name + '.geo', self.name + '.surfaces',
                           self.name + '.materials', self.name + '.sets', self.name + '.boundaries',
                           field_initialisation_str]]
            data = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
            filename = self.output_dir + self.name + '.dom.dat'
            print('Writing out ' + filename)
            with open(filename, 'w') as f:
                f.write(data)

            # Write out required geometric fields in Alya format
            print('Writing out geometric files to ' + self.output_dir + ' in Alya format')
            filename = self.output_dir + self.name + '.geo'
            with open(filename, 'w') as f:
                f.write('ELEMENTS\n')
                field_idx = np.arange(1, self.geometry.number_of_elements + 1)
                field_data = self.geometry.tetrahedrons + 1
                field_idx = np.array(field_idx).astype(int)
                field_data = np.array(field_data)
                for i in range(len(field_idx)):
                    f.write(str(field_idx[i]))
                    if len(field_data.shape) == 1:
                        f.write('\t' + str(field_data[i]))
                    else:
                        for j in range(field_data.shape[1]):
                            f.write('\t' + str(field_data[i, j]))
                    f.write('\n')
                f.write('END_ELEMENTS\n')
                f.write('COORDINATES\n')
                field_idx = np.arange(1, self.geometry.number_of_nodes + 1)
                field_data = self.geometry.nodes_xyz
                field_idx = np.array(field_idx).astype(int)
                field_data = np.array(field_data)
                for i in range(len(field_idx)):
                    f.write(str(field_idx[i]))
                    if len(field_data.shape) == 1:
                        f.write('\t' + str(field_data[i]))
                    else:
                        for j in range(field_data.shape[1]):
                            f.write('\t' + str(field_data[i, j]))
                    f.write('\n')
                f.write('END_COORDINATES\n')

            filename = self.output_dir + self.name + '.surfaces'
            with open(filename, 'w') as f:
                f.write('BOUNDARIES\n')
                field_idx = np.arange(1, self.geometry.number_of_triangles + 1)
                field_data = self.geometry.triangles.astype(int) + 1
                field_idx = np.array(field_idx).astype(int)
                field_data = np.array(field_data)
                for i in range(len(field_idx)):
                    f.write(str(field_idx[i]))
                    if len(field_data.shape) == 1:
                        f.write('\t' + str(field_data[i]))
                    else:
                        for j in range(field_data.shape[1]):
                            f.write('\t' + str(field_data[i, j]))
                    f.write('\n')
                f.write('END_BOUNDARIES\n')
            write_alya_field(filename=self.output_dir + self.name + '.materials',
                             field_idx=np.arange(1, self.materials.dict['tetra'].shape[0] + 1),
                             field_data=self.materials.dict['tetra'])
            write_alya_field(filename=self.output_dir + self.name + '.boundaries',
                             field_idx=np.arange(1, self.boundary_element_fields.dict[
                                 'mechanical-element-boundary-label'].shape[0] + 1),
                             field_data=self.boundary_element_fields.dict['mechanical-element-boundary-label'].astype(
                                 int))
            write_alya_field(filename=self.output_dir + self.name + '.sets',
                             field_idx=np.arange(1,self.boundary_element_fields.dict['mechanical-element-boundary-label'].shape[0] + 1),
                             field_data=self.boundary_element_fields.dict['mechanical-element-boundary-label'].astype(
                                 int))

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
            filename = self.template_dir + self.version + '.exm.dat'
            keys = ["stimulus_field_number", "monodomain_delay", "number_of_ecg_electrodes",
                    "ecg_electrodes_coordinates"]
            stimulus_field_number = self.simulation_dict['field_names'].index(
                self.simulation_dict['stimulus_field_name']) + 1
            num_electrodes = self.node_fields.dict['electrode_xyz'].shape[0]
            ecg_electrode_coordinates_str = ''
            for i in range(num_electrodes):
                ecg_electrode_coordinates_str = ecg_electrode_coordinates_str + '\n\t\t' + \
                                                str(self.node_fields.dict['electrode_xyz'][i, 0]) + ' ' + \
                                                str(self.node_fields.dict['electrode_xyz'][i, 1]) + ' ' + \
                                                str(self.node_fields.dict['electrode_xyz'][i, 2])
            insert_data = [[stimulus_field_number, self.simulation_dict['prestress_time'], num_electrodes,
                           ecg_electrode_coordinates_str]]
            data = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
            # Exmedi properties string
            subtemplate_filename = self.template_dir + self.version + '.subtemplate.exmedi_property_template'
            num_materials = np.amax(self.materials.dict['tetra'].astype(int))
            keys = ["material_idx", "sigma_f", "sigma_s", "sigma_n", "cell_model", "cell_initialisation_txt_file_name"]
            insert_data = []
            for material_i in range(num_materials):
                insert_data.append([material_i + 1, self.simulation_dict['sigma'][material_i][0],
                                    self.simulation_dict['sigma'][material_i][1],
                                    self.simulation_dict['sigma'][material_i][2],
                                    self.simulation_dict['cell_model'][material_i],
                                    self.simulation_dict['cell_filename'][material_i]])
            data = data.replace('<<exmedi_properties_str>>', self.template(filename=subtemplate_filename,
                                                                              keys=keys, data=insert_data,
                                                                              num_duplicates=num_materials))
            # Exmedi scaling factor string
            subtemplate_filename = self.template_dir + self.version + '.subtemplate.ionic_scaling_factors_template'
            keys = ["scaling_name", "scaling_name_field_number"]
            insert_data = []
            sf_fields = [x for x in list(self.simulation_dict['field_names']) if 'sf_' in x]
            num_sf = len(sf_fields)
            for field_i in range(num_sf):
                insert_data.append([sf_fields[field_i].upper(), list(self.simulation_dict['field_names']).index(sf_fields[field_i])+1])
            data = data.replace('<<exmedi_scaling_factor_str>>', self.template(filename=subtemplate_filename,
                                                                                  keys=keys, data=insert_data,
                                                                                  num_duplicates=num_sf))
            # Exmedi post process
            subtemplate_filename = self.template_dir + self.version + '.subtemplate.postprocess_template'
            keys = ["post_process_name", "post_process_period"]
            insert_data = []
            num_postprocess = len(self.simulation_dict['exmedi_output'])
            for var_i in range(num_postprocess):
                insert_data.append([self.simulation_dict['exmedi_output'][var_i], self.simulation_dict['exmedi_output_period'][var_i]])
            data = data.replace('<<post_process_exmedi_str>>', self.template(filename=subtemplate_filename,
                                                                                keys=keys, data=insert_data,
                                                                                num_duplicates=num_postprocess))
            # Write out Exmedi file
            filename = self.output_dir + self.name + '.exm.dat'
            print('Writing out ' + filename)
            with open(filename, 'w') as f:
                f.write(data)

    def write_sld_dat(self):
        if self.version == 'alya-compbiomed2':
            # Solidz properties
            filename = self.template_dir + self.version + '.subtemplate.solidz_property_template'
            keys = ["material_idx", "density", "Kct", "a", "b", "af", "bf", "as", "bs", "afs", "bfs", "prestress_field"]
            num_materials = np.amax(self.materials.dict['tetra']).astype(int)
            insert_data = []
            for material_i in range(num_materials):
                insert_data.append([material_i+1, self.simulation_dict['density'][material_i],
                                    self.simulation_dict['Kct'][material_i], self.simulation_dict['a'][material_i],
                                    self.simulation_dict['b'][material_i], self.simulation_dict['af'][material_i],
                                    self.simulation_dict['bf'][material_i], self.simulation_dict['as'][material_i],
                                    self.simulation_dict['bs'][material_i], self.simulation_dict['afs'][material_i],
                                    self.simulation_dict['bfs'][material_i], self.simulation_dict['field_names'].index('tv-element')])
            mechanical_properties_str = self.template(filename=filename, keys=keys, data=insert_data,
                                                      num_duplicates=num_materials)
            # Solidz cardiac cycle
            filename = self.template_dir + self.version + '.subtemplate.cavity_definition_template'
            keys = ["cavity_idx", "cavity_boundary_number", "cavity_basal_node_1", "cavity_basal_node_2", "pressure_field_idx"]
            num_cavities = len(self.simulation_dict['cavity_bcs'])
            insert_data = []
            for cavity_i in range(num_cavities):
                cavity_boundary_number = 0
                if self.simulation_dict['cavity_bcs'][cavity_i] == 'LV':
                    cavity_boundary_number = self.geometry.lv_endocardium
                    insert_data.append([cavity_i + 1, cavity_boundary_number,
                                        self.node_fields.dict['lv-cavity-nodes'][0].astype(int),
                                        self.node_fields.dict['lv-cavity-nodes'][1].astype(int), cavity_i])
                elif self.simulation_dict['cavity_bcs'][cavity_i] == 'RV':
                    cavity_boundary_number = self.geometry.rv_endocardium
                    insert_data.append([cavity_i+1, cavity_boundary_number,
                                        self.node_fields.dict['rv-cavity-nodes'][0].astype(int),
                                        self.node_fields.dict['rv-cavity-nodes'][1].astype(int), cavity_i])
            cavities_str = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=num_cavities)
            # Pressure string
            filename = self.template_dir + self.version + '.subtemplate.pressure_cycle_template'
            subkeys = ["prestress_t", "prestress_p", "diastasis_t", "diastasis_p",
                    "gain_error_prestress", "gain_derror_prestress", "end_diastole_t", "end_diastole_p",
                    "gain_error_contraction", "gain_derror_contraction", "arterial_compliance", "arterial_resistance",
                    "ejection_pressure_threshold", "gain_error_relaxation", "gain_derror_relaxation",
                    "filling_pressure_threshold"]
            keys = ["cycle_idx", "cavity_idx"] + subkeys + ["cycle_length"]
            insert_data = []
            for cavity_i in range(num_cavities):
                temp = [cavity_i + 1, cavity_i + 1]
                for key in subkeys:
                    temp.append(self.simulation_dict[key][cavity_i])
                temp.append(self.simulation_dict['cycle_length'])
                insert_data.append(temp)
            pressure_cycles_str = self.template(filename=filename, keys=keys, data=insert_data,
                                                num_duplicates=num_cavities)
            # Cardiac cycle template
            filename = self.template_dir + self.version + '.subtemplate.cardiac_cycle_template'
            keys = ["cavities_str", "pressure_cycles_str"]
            insert_data = [[cavities_str, pressure_cycles_str]]
            cardiac_cycle_str = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
            # Post processing
            filename = self.template_dir + self.version + '.subtemplate.postprocess_template'
            keys = ["post_process_name", "post_process_period"]
            insert_data = []
            num_postprocess = len(self.simulation_dict['solidz_output'])
            for var_i in range(num_postprocess):
                insert_data.append(
                    [self.simulation_dict['solidz_output'][var_i], self.simulation_dict['solidz_output_period'][var_i]])
            post_process_solidz_str = self.template(filename=filename, keys=keys, data=insert_data,
                                                    num_duplicates=num_postprocess)
            # Boundary codes
            filename = self.template_dir + self.version + '.subtemplate.solidz_bc_cavity_pressure_template'
            keys = ["cavity_boundary_number", "cavity_idx"]
            insert_data = []
            for cavity_i in range(num_cavities):
                if self.simulation_dict['cavity_bcs'][cavity_i] == 'LV':
                    cavity_boundary_number = self.geometry.lv_endocardium
                    insert_data.append([cavity_boundary_number, cavity_i + 1])
                elif self.simulation_dict['cavity_bcs'][cavity_i] == 'RV':
                    cavity_boundary_number = self.geometry.rv_endocardium
                    insert_data.append([cavity_boundary_number, cavity_i + 1])
            cavities_code = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=num_cavities)
            filename = self.template_dir + self.version + '.subtemplate.solidz_bc_pericardium_template'
            keys = ["epicardial_boundary_number","pericardial_stiffness"]
            insert_data = [[self.geometry.epicardium, self.simulation_dict['pericardial_stiffness']]]
            pericardial_code = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
            boundary_conditions_str = cavities_code + '\n' + pericardial_code
            # Replace all strings
            filename = self.template_dir + self.version + '.sld.dat'
            keys = ["mechanical_properties_str", "cardiac_cycle_str", "solidz_residual", "solidz_convergence_tolerance",
                    "post_process_solidz_str", "boundary_conditions_str"]
            insert_data = [[mechanical_properties_str, cardiac_cycle_str, self.simulation_dict['solidz_residual'],
                            self.simulation_dict['solidz_convergence_tolerance'], post_process_solidz_str,
                            boundary_conditions_str]]
            data = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
            # Write sld.dat
            filename = self.output_dir + self.name + '.sld.dat'
            print('Writing out ' + filename)
            with open(filename, 'w') as f:
                f.write(data)

    def write_ker_dat(self):
        if self.version == 'alya-compbiomed2':
            coupling_str = ''
            if 'SOLIDZ' in self.simulation_dict['physics'] and 'EXMEDI' in self.simulation_dict['physics']:
                coupling_str = '\n\tCOUPLING\n\t\tSOLIDZ EXMEDI\n\tEND_COUPLING'
            filename = self.template_dir + self.version + '.subtemplate.coupling_template'
            subkeys = ["eccoupling_model_name", "cal50", "tref_sheet_scaling", "tref_normal_scaling", "tref_scaling"]
            insert_data = []
            num_materials = np.amax(self.materials.dict['tetra']).astype(int)
            for material_i in range(num_materials):
                temp = [material_i + 1]
                for key in subkeys:
                    temp.append(self.simulation_dict[key][material_i])
                insert_data.append(temp)
            keys = ["material_idx"] + subkeys
            coupling_str_2 = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=num_materials)
            coupling_str = coupling_str + coupling_str_2
            # Replace all ker.dat
            filename = self.template_dir + self.version + '.ker.dat'
            keys = ["coupling_str", "fibre_field_number", "sheet_field_number", "normal_field_number",
                    "celltype_field_number"]
            insert_data = [[coupling_str,
                           self.simulation_dict['field_names'].index(self.simulation_dict['fibre_field_name'])+1,
                           self.simulation_dict['field_names'].index(self.simulation_dict['sheet_field_name'])+1,
                           self.simulation_dict['field_names'].index(self.simulation_dict['normal_field_name'])+1,
                           self.simulation_dict['field_names'].index(self.simulation_dict['celltype_field_name'])+1]]
            data = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
            # Write sld.dat
            filename = self.output_dir + self.name + '.ker.dat'
            print('Writing out ' + filename)
            with open(filename, 'w') as f:
                f.write(data)

    def write_post_dat(self):
        os.system('cp '+self.template_dir+self.version+'.post.alyadat '+self.output_dir+self.name+'.post.alyadat')

    def write_cell_txt(self):
        filename = self.template_dir + self.version + '.celltxt'
        subkeys = ["cell_number_of_beats", "cell_steady_state_type",
                "cell_steady_state_tolerance", "cell_model_modification_toggle"]
        insert_data = []
        num_materials = np.amax(self.materials.dict['tetra']).astype(int)
        for material_i in range(num_materials):
            temp = []
            for key in subkeys:
                temp.append(self.simulation_dict[key][material_i])
            insert_data.append(temp)
        # Get all sf_
        with open(filename, 'r') as f:
            data = f.readlines()
        sf_keys = []
        for i in range(len(data)):
            if 'sf_' in data[i]:
                sf_keys.append(data[i].split('<')[2].split('>')[0])
        keys = subkeys
        cell_types = ["ENDO", "MID", "EPI"]
        cell_type_str = list(np.array(cell_types)[np.unique(self.node_fields.dict['cell-type']).astype(int) - 1])
        cell_type_str = ' '.join(map(str, cell_type_str))
        keys = keys + sf_keys
        for material_i in range(num_materials):
            for sf_key in sf_keys:
                if sf_key in list(self.simulation_dict.keys()):
                    string = ' '.join(map(str, self.simulation_dict[sf_key][material_i]))
                    print(string)
                    insert_data[material_i][:] = insert_data[material_i][:] + [string]
                else:
                    insert_data[material_i][:] = insert_data[material_i][:] + ['1 1 1']
            insert_data[material_i][:] = insert_data[material_i][:] + \
                                         [str(int(self.simulation_dict['cycle_length'] * 1000))]
            insert_data[material_i][:] = insert_data[material_i][:] + [cell_type_str]
        keys = keys + ["cycle_length", "cell_type_str"]
        print(insert_data[1][:])
        # Write to txtfiles
        for material_i in range(num_materials):
            data = self.template(filename=filename, keys=keys, data=[insert_data[material_i][:]], num_duplicates=1)
            output_filename = self.output_dir + self.simulation_dict['cell_filename'][material_i]
            print('Writing out ' + output_filename)
            with open(output_filename, 'w') as f:
                f.write(data)

    def write_job_scripts(self):
        filename = self.job_template_dir + self.job_version + '.main.py'
        keys = ['alya_exec_path', 'casename']
        insert_data = [[self.simulation_dict['alya_exec_path'], self.simulation_dict['name']]]
        data = self.template(filename=filename, keys=keys, data=insert_data, num_duplicates=1)
        filename = self.output_dir + 'main.py'
        with open(filename, 'w') as f:
            print('Writing out ' + filename)
            f.write(data)

        filename = self.job_template_dir + self.job_version + '.run_job.cmd'
        keys = ['job_name', 'job_time', 'computational_nodes', 'tasks_per_node', 'computational_cores', 'job_type']
        tasks_per_node = 0
        job_type = ''
        if self.job_version == 'jureca':
            tasks_per_node = 128
            job_type = 'dc-cpu'
        insert_data = [[self.simulation_dict['name'], str(int(np.ceil(self.simulation_dict['time_end']*2.))),
                       np.ceil(self.simulation_dict['computational_cores']/tasks_per_node).astype(int), tasks_per_node,
                       self.simulation_dict['computational_cores'], job_type]]
        data = self.template(filename, keys=keys, data=insert_data, num_duplicates=1)
        filename = self.output_dir + 'run_job.cmd'
        with open(filename, 'w') as f:
            print('Writing out ' + filename)
            f.write(data)



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
