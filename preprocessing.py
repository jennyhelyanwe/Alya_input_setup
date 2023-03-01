# This takes as input csv files in myformat and produces fields required for Alya simulations
# CSV inputs and CSV outputs
from myformat import *
class preprocessing:
    def __init__(self, name, input_dir, output_dir):
        self.name = name
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        self.input_dir = input_dir
        self.output_dir = output_dir
        self._read_input_from_csv()

    def _read_input_from_csv(self):
        self.geometry = geometry(self.name)
        self.geometry.read_csv_to_attributes(self.input_dir)
        self.node_fields = node_fields(self.name)
        self.node_fields.read_csv_to_attributes(self.input_dir)
        self.element_fields = element_fields(self.name)
        self.element_fields.read_csv_to_attributes(self.input_dir)
        self.surface_fields = surface_fields(self.name)
        self.surface_fields.read_csv_to_attributes(self.input_dir)





