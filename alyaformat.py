from myformat import *
import json

class AlyaFormat:
    def __init__(self, name, geometry_and_fields_input_dir, simulation_json_file):
        self.version = 'alya-compbiomed2' # 'Alya_multiple_BZRZ_models'
        self.name = name
        self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.geometry_and_fields_input_dir = geometry_and_fields_input_dir
        self.geometry = None
        self.node_fields = None
        self.element_fields = None
        self.import_geometry_and_fields()
        # self.space_dimensions = 3 # 3D problem, 2 for 2D problem
        # self.simulation_files_dir = simulation_files_dir
        # self.simulation_start_time = 0
        # self.simulation_end_time = 0.8
        # self.time_step_function = None
        # self.run_type = 'PRELIMINARY'
        # self.warm_start_writing_steps = 50000
        # self.physics = np.array(['SOLIDZ', 'EXMEDI'])
        # self.coupling_model = None
        # self.materials= None
        # self.simulation = None
        # self.simulation_dict = None
        # self.geometric_scaling = [1., 1., 1.]
        # self.section_divider = '$------------------------------------------------------\n'


    def initialise_default_simulation_parameters(self):
        self.simulation_dict = {'start_time': 0, 'end_time': 0.8,
                'time_step_function': np.array([[0.8, 2e-5]]),
                'run_type': 'PRELIMINARY',
                'warm_start_writing_steps': 50000}

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

    def import_simulation_parameters(self):
        self.simulation = Fields(self.name, field_type='simulationfield')
        self.simulation.read_csv_to_attributes(self.simulation_files_dir)

    def generate_fields(self):
        # Generate required fields for Alya simulations.
        print('Make use of field_evaluation_functions to generate Alya fields...')


    def write_dat(self):
        with open(self.version+'.dat', 'r') as f:
            data = f.readlines()
        iterative_keys = []
        for key in self.simulation_dict.keys():
            if 'iter_' in key:
                iterative_keys.append(key)
                full_string = ''
                with open(self.simulation_dict[key][1], 'r') as f:
                    data = f.readlines()
                for key_i in range(self.simulation_dict[key][0]):
                    temp = data.replace("<<iter_idx>>", key_i)

                    to_append = temp.replace("")

                    full_string = full_string +
                with open(self.simulation_dict[key][])


    def write_dat(self):
        if self.version == 'alya-compbiomed2' or self.version == 'Alya_multiple_BZRZ_models':
            with open(self.simulation_files_dir + self.name + '.dat', 'w') as f:
                f.write(self.section_divider)
                f.write('RUN DATA\n')
                f.write('\tALYA: '+self.name+'\n')
                f.write('\tLIVE_INFO: screen\n\tPOSTPROCESS: domain\n')
                if self.run_type=='CONTINUE':
                    f.write('\tRUNTYPE: CONTINUE\n')
                elif self.run_type=='PRELIMINARY':
                    f.write('\tRUNTYPE: PRELIMINARY, FREQUENCY = ' + str(self.warm_start_writing_steps)+'\n')
                elif self.run_type == 'PRELIMINARY CONTINUE':
                    f.write('\tRUNTYPE: CONTINUE\n\tRUNTYPE: PRELIMINARY, FREQUENCY = ' + str(self.warm_start_writing_steps)+'\n')
                else:
                    raise 'write_dat: Run type '+self.run_type+' not found.'
                f.write(self.section_divider)
                f.write('PROBLEM_DATA\n')
                f.write('\tMAXIMUM_NUMBER_GLOBAL: 1\n')
                f.write('\tTIME_COUPLING: GLOBAL, PRESCRIBED\n')
                f.write('\tTIME_INTERVAL: '+str(self.simulation_start_time)+', '+str(self.simulation_end_time)+'\n')
                if not self.time_step_function is None:
                    f.write('\tFUNCTION: DISCRETE\n')
                    f.write('\tDTIME_FUNCTION:\n')
                    f.write('\t\tTIME_SHAPE: DISCRETE\n')
                    f.write('\t\tSHAPE_DEFINITION:\n')
                    f.write('\t\t\t'+str(self.time_step_function.shape[0])+'\n')
                    for shape_i in range(self.time_step_function.shape[0]):
                        f.write('\t\t\t'+str(self.time_step_function[shape_i,0])+'\t'+str(self.time_step_function[shape_i,1])+'\n')
                    f.write('\t\tEND_SHAPE_DEFINITION\n')
                    f.write('\tEND_DTIME_FUNCTION\n')
                f.write('\tNUMBER_OF_STEPS_(MAX): 1000000000\n')
                f.write(self.section_divider)
                for physics_i in range(self.physics.shape[0]):
                    if self.physics[physics_i] == 'SOLIDZ':
                        f.write('\tSOLIDZ_PROBLEM:  ON\n')
                        f.write('\tEND_SOLIDZ')
                    if self.physics[physics_i] == 'EXMEDI':
                        f.write('\tEXMEDI_PROBLEM:  ON\n')
                        f.write('\tEND_EXMEDI')
                f.write('\tPARALL_SERVICE: ON\n\t\tOUTPUT_FILE: NO\n\t\tPOSTPROCESS: MASTER\n\t\tPARTITION: FACES\n\tEND_PARALL\n')
                f.write(self.section_divider)
                f.write('END_PROBLEM_DATA\n')
                f.write(self.section_divider)

    def write_dom_dat(self):
        if not self.geometry is None:
            if self.version == 'alya-compbiomed2':
                with open(self.simulation_files_dir + self.name + '.dom.dat', 'w') as f:
                    f.write(self.section_divider)
                    f.write('DIMENSIONS\n')
                    f.write('\tNODAL_POINTS:\t'+str(self.geometry.number_of_nodes)+'\n')
                    f.write('\tELEMENTS:\t'+str(self.geometry.number_of_elements)+'\n')
                    f.write('\tSPACE_DIMENSIONS:\t' + str(self.space_dimensions) + '\n')
                    f.write('\tNODES:\t'+str(self.geometry.tetrahedrons[0].shape[0])+'\n')
                    f.write('\tBOUNDARIES:\t'+str(self.geometry.number_of_triangles)+'\n')
                    f.write('\tSKEW_SYSTEMS:\t0\n\tSLAVES:\t0\n\tNO_SETS\n')
                    f.write('\nMATERIALS = '+str(np.amax(self.materials))+'\n')
                    f.write('\nFIELDS = '+str(self.node_fields.number_of_fields +
                                              self.element_fields.number_of_fields)+'\n')
                    varnames = self.node_fields.dict.keys()
                    for field_i in range(self.node_fields.number_of_fields):
                        varname = varnames[field_i]
                        assert self.node_fields.dict[varname].shape[0] == self.geometry.number_of_nodes, \
                            'Field: ' + varname + ' does not match in first dimension with the number of nodes ' \
                                                  'in the geometry'
                        f.write('\t\tFIELD = '+str(field_i)+', DIMENSION = '+str(self.node_fields.dict[varname].shape[1])+', NODES\n')
                    varnames = self.element_fields.dict.keys()
                    for field_i in range(self.element_fields.number_of_fields):
                        varname = varnames[field_i]
                        assert self.node_fields.dict[varname].shape[0] == self.geometry.number_of_elements, \
                            'Field: ' + varname + ' does not match in first dimension with the number of elements ' \
                                                  'in the geometry'
                        f.write('\t\tFIELD = ' + str(field_i) + ', DIMENSION = ' + str(
                            self.node_fields.dict[varname].shape[1]) + ', ELEMENTS\n')
                    f.write('\tEND_FIELDS\n')
                    f.write('END_DIMENSIONS\n')
                    f.write(self.section_divider)
                    f.write('STRATEGY\n\tINTEGRATION_RULE:\tOPEN\n\tDOMAIN_INTEGRATION_POINTS: 1\n')
                    f.write('\tSCALE: XSCAL='+str(self.geometric_scaling[0])+', YSCAL='+str(self.geometric_scaling[1])+
                            ', ZSCAL='+str(self.geometric_scaling)+'\n')
                    f.write('\tTRANSLATION: XTRAN='+str(self.geometric_scaling[0])+', YTRAN='+str(self.geometric_scaling[1])+
                            ', ZTRAN='+str(self.geometric_scaling)+'\n')
                    f.write('\tEXTRAPOLATE_BOUNDARY_CONDTIONS: ON\n')
                    f.write('\tBOUNDARY_ELEMENT: OFF\n')
                    f.write('END_STRATEGY\n')
                    f.write(self.section_divider)
                    f.write('GEOMETRY\n')
                    f.write('\tINCLUDE '+str(self.name)+'.geo\n')
                    f.write('\tINCLUDE '+str(self.name)+'.surface\n')
                    f.write('\tINCLUDE ' + str(self.name) + '.materials\n')
                    f.write('END_GEOMETRY\n')
                    f.write(self.section_divider)
                    f.write('SETS\n')
                    f.write('\tINCLUDE '+str(self.name)+'.sets\n')
                    f.write('END_SETS\n')
                    f.write(self.section_divider)
                    f.write('BOUNDARY_CONDITIONS\n')
                    f.write('\tINCLUDE '+str(self.name)+'.boundaries\n')
                    f.write('END_BOUNDARY_CONDITIONS\n')
                    f.write(self.section_divider)
                    f.write('FIELDs, NUMBER = '+str(self.node_fields.number_of_fields+self.element_fields.number_of_fields)+'\n')
                    for field_i in range(self.element_fields.number_of_fields):
                        varname = varnames[field_i]
                        f.write('\tFIELD = '+str(field_i)+'\n')
                        f.write('\t\tINCLUDE '+str(self.name)+'.'+varname+'\n')
                        f.write('\tEND_FIELD\n')
                    f.write('END_FIELDS\n')
        else:
            raise 'Please read in geometry information first before writing .dom.dat file. '


    def write_exm_dat(self):




