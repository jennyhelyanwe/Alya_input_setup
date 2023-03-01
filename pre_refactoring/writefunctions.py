import os, sys
import numpy as np
import vtk
from vtk.util import numpy_support as VN

class O():
    def __init__(self, I, name, output_dir):
        self.name = name
        self.output_dir = output_dir
        self.I = I


    def _write_alya_mesh(self):
        print('# Write .geo file for Alya input')
        with open(self.output_dir + self.name + '.geo', 'w') as f:
            f.write('ELEMENTS\n')
            for i in range(0, self.I.nelems):
                f.write(str(i+1) + '\t' + str(self.I.elems[i,0]) + '\t' + str(self.I.elems[i,1])+ '\t' + str(self.I.elems[i,2]) + '\t'+ str(self.I.elems[i,3]) + '\n')
            f.write('END_ELEMENTS\n')
            f.write('COORDINATES\n')
            for i in range(0, len(self.I.nodes)):
                f.write(str(i+1) + '\t' + str(self.I.nodes[i,0]) + '\t' + str(self.I.nodes[i,1]) + '\t' + str(self.I.nodes[i,2]) + '\n')
            f.write('END_COORDINATES\n')


    def _write_alya_surface_boundary(self):
        print('# Writing to Alya surface and boundary files')
        with open(self.output_dir + self.name + '.surfaces','w') as f:
            f.write('BOUNDARIES\n')
            for i in range(0, len(self.I.faces)):
                f.write(str(i+1)+'\t'+str(self.I.faces[i,0]+1)+'\t'+str(self.I.faces[i,1]+1)+'\t'+str(self.I.faces[i,2]+1)+'\n')
            f.write('END_BOUNDARIES\n')
        with open(self.output_dir + self.name + '.boundaries','w') as f:
            f.write('ON_BOUNDARIES\n')
            for i in range(0, len(self.I.faces_label)):
                f.write(str(i+1)+'\t'+str(self.I.faces_label[i])+'\n')
            f.write('END_ON_BOUNDARIES\n')


    def _write_alya_celltype(self):
        print('# Writing to Alya cell type file')
        self.I.celltype = [1]*self.I.nnodes
        for i in range(0, self.I.nnodes):
            if self.I.transmural_coordinate[i] > 0.7:
                self.I.celltype[i] = 3 # epicardial
            elif self.I.transmural_coordinate[i] > 0.3:
                self.I.celltype[i] = 2 # mid-myocardial
            else:
                self.I.celltype[i] = 1 # endocardial
        with open(self.output_dir+self.name+'.celltype','w') as f:
            for i in range(0, self.I.nnodes):
                f.write(str(i+1)+'\t'+str(int(self.I.celltype[i]))+'\n')
