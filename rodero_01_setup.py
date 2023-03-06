from vtk_to_csv import *

vtk_dir = '/Users/Jenny/Documents/Codes/Alya/Simulations/RoderoNiedererMeshHealthy/01'
csv_dir = '/Users/Jenny/PycharmProjects/Alya_input_setup/csv/'
vtk_name = '01_bivent_only'

vtk_reader = VTKtoCSV(name=vtk_name, input_dir=vtk_dir, output_dir=csv_dir)
