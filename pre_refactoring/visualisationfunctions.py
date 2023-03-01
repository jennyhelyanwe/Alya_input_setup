import os, sys

class V():
    def __init__(self, i, name, visualisation_dir, visualisation_format):
        self.name = name
        self.visualisation_dir = visualisation_dir
        self.visualisation_format = visualisation_format
        self.I = i

    def _ensight_export_mesh(self, nodes, elems):
        print('ENSIGHT EXPORT: '+self.visualisation_dir+self.name+'.ensi.geo')
        with open(self.visualisation_dir+self.name+'.ensi.geo', 'w') as f:
            f.write('Problem name:  '+self.name+'\nGeometry file\nnode id given\nelement id given\npart\n\t1\nVolume Mesh\ncoordinates\n'+str(len(nodes))+'\n')
            for i in range(0, len(nodes)):
                f.write(str(i+1)+'\n')
            for c in [0,1,2]:
                for i in range(0, len(nodes)):
                    f.write(str(nodes[i,c])+'\n')
            f.write('tetra4\n  '+str(len(elems))+'\n')
            for i in range(0, len(elems)):
                f.write('  '+str(i+1)+'\n')
            for i in range(0, len(elems)):
                f.write(str(elems[i,0])+'\t'+str(elems[i,1])+'\t'+str(elems[i,2])+'\t'+str(elems[i,3])+'\n')
        with open(self.visualisation_dir+self.name+'.ensi.case', 'w') as f:
            f.write('#\n# Alya generated post-process files\n# Ensight Gold Format\n#\n# Problem name:\t'+self.name+'\n#\n')
            f.write('FORMAT\ntype:\tensight gold\nGEOMETRY\nmodel:\t1\t'+self.name+'.ensi.geo\nVARIABLE\n')


    def _ensight_export_scalar_per_node(self, field_name, data):
        print('ENSIGHT EXPORT: '+self.visualisation_dir+self.name+'.ensi.'+field_name)
        with open(self.visualisation_dir+self.name+'.ensi.'+field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Scalar per-node variables file\npart\n\t1\ncoordinates\n')
            for i in range(0, len(data)):
                f.write(str(data[i])+'\n')
        if os.path.exists(self.visualisation_dir+self.name+'.ensi.case'):
            with open(self.visualisation_dir+self.name+'.ensi.case', 'a') as f:
                f.write('scalar per node:	1	'+field_name+'	'+self.name+'.ensi.'+field_name+'\n')
        else:
            with open(self.visualisation_dir+self.name+'.ensi.case', 'w') as f:
                f.write('#\n# Alya generated post-process files\n# Ensight Gold Format\n#\n# Problem name:\t'+self.name+'\n#\n')
                f.write('FORMAT\ntype:\tensight gold\nGEOMETRY\nmodel:\t1\t'+self.name+'.ensi.geo\nVARIABLE\n')
                f.write('scalar per node:	1	'+field_name+'	'+self.name+'.ensi.'+field_name+'\n')


    def _ensight_export_vector_per_node(self, field_name, data):
        print('ENSIGHT EXPORT: '+self.visualisation_dir+self.name+'.ensi.'+field_name)
        with open(self.visualisation_dir+self.name+'.ensi.'+field_name, 'w') as f:
            f.write('Alya Ensight Gold --- Vector per-node variables file\npart\n1\ncoordinates\n')
            for c in [0, 1, 2]:
                for i in range(0, len(data)):
                    f.write(str(data[i,c])+'\n')
        if os.path.exists(self.visualisation_dir+self.name+'.ensi.case'):
            with open(self.visualisation_dir+self.name+'.ensi.case', 'a') as f:
                f.write('vector per node:	1	'+field_name+'	'+self.name+'.ensi.'+field_name+'\n')
        else:
            with open(self.visualisation_dir+self.name+'.ensi.case', 'w') as f:
                f.write('#\n# Alya generated post-process files\n# Ensight Gold Format\n#\n# Problem name:\t'+self.name+'\n#\n')
                f.write('FORMAT\ntype:\tensight gold\nGEOMETRY\nmodel:\t1\t'+self.name+'.ensi.geo\nVARIABLE\n')
                f.write('vector per node:	1	'+field_name+'	'+self.name+'.ensi.'+field_name+'\n')
