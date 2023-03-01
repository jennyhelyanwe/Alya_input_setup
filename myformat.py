import pandas as pd
import os

class dictionary_to_object:
    def __init__(self, **entries):
        self.__dict__.update(entries)

class geometry:
    def __init__(self, name):
        # Geometry initialisations
        self.name = name
        self.number_of_nodes = 0
        self.nodes_xyz_1 = []
        self.nodes_xyz_2 = []
        self.nodes_xyz_3 = []
        self.number_of_elements = 0
        self.tetrahedrons_1 = []
        self.tetrahedrons_2 = []
        self.tetrahedrons_3 = []
        self.tetrahedrons_4 = []
        self.number_of_surfaces = 0
        self.triangles = []
    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        dataframe = pd.DataFrame([vars(self)])
        dataframe.to_csv(output_dir + self.name + '_geometry.csv')
    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        df = pd.read_csv(input_dir + self.name + '_geometry.csv')
        self = dictionary_to_object(**df.to_dict())

class node_fields:
    def __init__(self, name):
        self.name = name
        # Field initialisations
        self.fibres_1 = []
        self.fibres_2 = []
        self.fibres_3 = []
        self.sheets_1 = []
        self.sheets_2 = []
        self.sheets_3 = []
        self.normal_1 = []
        self.normal_2 = []
        self.normal_3 = []
        self.lvrv = []
        self.tm = []
        self.ab = []
        self.rt = []
        self.celltype = []
        self.ab_Gks_sf = []
        self.stimulus = []

    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        dataframe = pd.DataFrame([vars(self)])
        dataframe.to_csv(output_dir + self.name + '_node_fields.csv')

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        df = pd.read_csv(input_dir + self.name + '_node_fields.csv')
        self = dictionary_to_object(**df.to_dict())

class element_fields:
    def __init__(self, name):
        self.name = name
        self.materials = []
        self.lvrv = []
    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        dataframe = pd.DataFrame([vars(self)])
        dataframe.to_csv(output_dir + self.name + '_element_fields.csv')

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        df = pd.read_csv(input_dir + self.name + '_element_fields.csv')
        self = dictionary_to_object(**df.to_dict())


class surface_fields:
    def __init__(self, name):
        self.name = name
        self.boundaries = []
    def save_to_csv(self, output_dir):
        if output_dir[-1] != '/':
            output_dir = output_dir+'/'
        dataframe = pd.DataFrame([vars(self)])
        dataframe.to_csv(output_dir + self.name + '_surface_fields.csv')

    def read_csv_to_attributes(self, input_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        df = pd.read_csv(input_dir + self.name + '_surface_fields.csv')
        self = dictionary_to_object(**df.to_dict())


if __name__ == '__main__':
    wd = os.getcwd()
    geometry = geometry('test')
    geometry.save_to_csv(wd)
    node_fields = node_fields('test')
    node_fields.save_to_csv(wd)
    element_fields = element_fields('test')
    element_fields.save_to_csv(wd)
    surface_fields = surface_fields('test')
    surface_fields.save_to_csv(wd)