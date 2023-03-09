import numpy
from myformat import Geometry, NodeFields, ElementFields
from pyDOE import lhs
from SALib.sample import saltelli
from matplotlib import pyplot as plt

class SensivityAnalysis:
    def __init__(self, vtk_name, output_name, input_dir, output_dir):
        if input_dir[-1] != '/':
            input_dir = input_dir + '/'
        if output_dir[-1] != '/':
            output_dir = output_dir + '/'
        self.vtk_name = vtk_name
        self.output_name = output_name
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.geometry = Geometry(self.output_name)
        self.node_fields = NodeFields(self.output_name)
        self.element_fields = ElementFields(self.output_name)

    def sampling(self, ranges, sample_size, ranges_names,method='lhs', verbose=False):
        if method == 'lhs':
            sample = lhs(ranges.shape[0], samples=sample_size, criterion='corr')
            for i in range(ranges.shape[0]):
                sample[:,i] = sample[:,i] * (ranges[i,1] - ranges[i,0]) + ranges[i,0]
        elif method == 'saltelli':
            problem = {
                'num_vars': ranges.shape[0],
                'names': ranges_names,
                'bounds': ranges
            }
            sample = saltelli.sample(problem, sample_size)
        else:
            raise('Sampling method: '+str(method)+' not recognised')
        if verbose:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            ax.scatter(sample[:,0], sample[:,1], sample[:,2])
            ax.set_xlabel(ranges_names[0])
            ax.set_ylabel(ranges_names[1])
            ax.set_zlabel(ranges_names[2])
            ax.set_title('Sampling method: '+str(method)+', n='+str(sample.shape[0]))
            plt.show(block=False)
        return sample

