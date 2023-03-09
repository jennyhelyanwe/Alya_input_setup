class AlyaFormat:
    def __init__(self, name, input_dir, output_dir, geometry_dir):
        self.name = name
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.geometry_dir = geometry_dir
        self.simulation_time = 0
        self.time_step_function = None
        self.run_type = 'CONTINUE'
        self.warm_start_writing_steps = 50000
        self.physics = None
        self.coupling_model = None
        self.geometry = None
        self.node_fields = None
        self.element_fields = None
        self.material_fields = None
        self.geometric_scaling = [1., 1., 1.]

    def read_geometry(self):
        # Read geometry from geometry directory CSV files
        print('Read in geometry and fields from CSV files.')

    def generate_fields(self):
        # Generate required fields for Alya simulations.
        print('Make use of field_evaluation_functions to generate Alya fields...')

    def write_dat(self):
        with open(self.output_dir+self.name+'.dat', 'w') as f:
            f.write('RUN DATA\n')
            f.write('\tALYA: '+self.name+'\n')
