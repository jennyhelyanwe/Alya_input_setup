import os
import numpy as np
import json
from scipy.optimize import minimize
import time
import gmsh
from meshstructure import MeshStructure


class CalibrateCV(MeshStructure):
    def __init__(self, name, geometric_data_dir, calibration_dir, alya_executable_path, personalisation_data_dir, verbose):
        super().__init__(name=name, geometric_data_dir=geometric_data_dir, verbose=verbose)
        edges = []
        for element_i in range(self.geometry.number_of_elements):
            edges.append(
                [self.geometry.tetrahedrons[element_i, 0], self.geometry.tetrahedrons[element_i, 1]])
            edges.append(
                [self.geometry.tetrahedrons[element_i, 1], self.geometry.tetrahedrons[element_i, 2]])
            edges.append(
                [self.geometry.tetrahedrons[element_i, 2], self.geometry.tetrahedrons[element_i, 3]])
            edges.append(
                [self.geometry.tetrahedrons[element_i, 3], self.geometry.tetrahedrons[element_i, 0]])
            edges.append(
                [self.geometry.tetrahedrons[element_i, 1], self.geometry.tetrahedrons[element_i, 3]])
            edges.append(
                [self.geometry.tetrahedrons[element_i, 0], self.geometry.tetrahedrons[element_i, 2]])
        edges = np.unique(np.sort(edges, axis=1), axis=0)
        edge_lengths = np.linalg.norm(self.geometry.nodes_xyz[edges[:,0],:] - self.geometry.nodes_xyz[edges[:,1],:],
                                      axis=1)
        mean_edge_length = np.mean(edge_lengths)
        self.lc = mean_edge_length
        best_parameters = np.loadtxt(personalisation_data_dir+name+'_inferred-best-parameter.csv', delimiter=',', skiprows=1)
        self.target_cvs = best_parameters[0:3]*1000  # convert to cm/s
        print('Targe CVs: ', self.target_cvs)
        self.name = name
        self.calibration_dir = calibration_dir + name + '/'
        # self.simulation_dict = json.load(open(simulation_json_file, 'r'))
        self.alya_executable_path = alya_executable_path

        width = 0.7
        cube_length = 1.0
        self.create_tissue_block(lc=self.lc, h=cube_length, w=cube_length, d=cube_length)
        self.local_node_list, self.local_node_list_coords, self.local_weight_list = self.get_element_probes(self.lc, width,cube_length)
        conductivities = minimize(self.objective_function, [0.001, 0.001, 0.001], method='L-BFGS-B',
                                  bounds=[[0.0005, 0.005], [0.0005, 0.005], [0.0005, 0.005]], options={'eps': 0.0005})
        print('Optimised conductivities: ')
        print(conductivities)

    def create_tissue_block(self, lc, h, w, d):
        gmsh.initialize()
        point1 = gmsh.model.geo.add_point(0, 0, 0, lc)
        point2 = gmsh.model.geo.add_point(w, 0, 0, lc)
        point3 = gmsh.model.geo.add_point(w, h, 0, lc)
        point4 = gmsh.model.geo.add_point(0, h, 0, lc)
        point5 = gmsh.model.geo.add_point(0, 0, d, lc)
        point6 = gmsh.model.geo.add_point(w, 0, d, lc)
        point7 = gmsh.model.geo.add_point(w, h, d, lc)
        point8 = gmsh.model.geo.add_point(0, h, d, lc)

        line1 = gmsh.model.geo.add_line(point1, point2)
        line2 = gmsh.model.geo.add_line(point2, point3)
        line3 = gmsh.model.geo.add_line(point3, point4)
        line4 = gmsh.model.geo.add_line(point4, point1)
        line5 = gmsh.model.geo.add_line(point5, point6)
        line6 = gmsh.model.geo.add_line(point6, point7)
        line7 = gmsh.model.geo.add_line(point7, point8)
        line8 = gmsh.model.geo.add_line(point8, point5)
        line9 = gmsh.model.geo.add_line(point1, point5)
        line10 = gmsh.model.geo.add_line(point2, point6)
        line11 = gmsh.model.geo.add_line(point3, point7)
        line12 = gmsh.model.geo.add_line(point4, point8)

        # faces of cube:
        face1 = gmsh.model.geo.add_curve_loop([line4, line1, line2, line3])
        face2 = gmsh.model.geo.add_curve_loop([line8, line7, line6, line5])
        face3 = gmsh.model.geo.add_curve_loop([line11, line7, -line12, -line3])
        face4 = gmsh.model.geo.add_curve_loop([line10, line6, -line11, -line2])
        face5 = gmsh.model.geo.add_curve_loop([line9, line5, -line10, -line1])
        face6 = gmsh.model.geo.add_curve_loop([line12, line8, -line9, -line4])

        # surfaces of cube:
        plane_surface1 = gmsh.model.geo.add_plane_surface([face1])
        plane_surface2 = gmsh.model.geo.add_plane_surface([face2])
        plane_surface3 = gmsh.model.geo.add_plane_surface([face3])
        plane_surface4 = gmsh.model.geo.add_plane_surface([face4])
        plane_surface5 = gmsh.model.geo.add_plane_surface([face5])
        plane_surface6 = gmsh.model.geo.add_plane_surface([face6])

        surface_loop = gmsh.model.geo.add_surface_loop([plane_surface1,plane_surface2,plane_surface3,plane_surface4,
                                                         plane_surface5,plane_surface6])
        volume = gmsh.model.geo.add_volume([surface_loop])
        # gmsh.model.geo.add_physical_group(1, [1], name="Line1")
        gmsh.model.geo.add_physical_group([line1], "Line1")
        gmsh.model.geo.add_physical_group([line2], "Line2")
        gmsh.model.geo.add_physical_group([line3], "Line3")
        gmsh.model.geo.add_physical_group([line4], "Line4")
        gmsh.model.geo.add_physical_group([line5], "Line5")
        gmsh.model.geo.add_physical_group([line6], "Line6")
        gmsh.model.geo.add_physical_group([line7], "Line7")
        gmsh.model.geo.add_physical_group([line8], "Line8")
        gmsh.model.geo.add_physical_group([line9], "Line9")
        gmsh.model.geo.add_physical_group([line10], "Line10")
        gmsh.model.geo.add_physical_group([line11], "Line11")
        gmsh.model.geo.add_physical_group([line12], "Line12")

        gmsh.model.geo.add_physical_group([plane_surface1], "Surface1")
        gmsh.model.geo.add_physical_group([plane_surface2], "Surface2")
        gmsh.model.geo.add_physical_group([plane_surface3], "Surface3")
        gmsh.model.geo.add_physical_group([plane_surface4], "Surface4")
        gmsh.model.geo.add_physical_group([plane_surface5], "Surface5")
        gmsh.model.geo.add_physical_group([plane_surface6], "Surface6")
        gmsh.model.geo.add_physical_group([volume], "Volume")

        gmsh.model.mesh.generate()
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
        gmsh.write("3D.msh")
        gmsh.finalize()


    def setup_alya_simulation(self):
        with open('3D_' + str(self.lc) + '.dims.dat', 'r') as f:
            lines = f.readlines()
        nnodes = int(lines[0].split()[-1])
        nelems = int(lines[1].split()[-1])
        nbound = int(lines[2].split()[-1])
        with open('3D.dom.dat', 'r') as f:
            data = f.readlines()
        data[3] = 'NODAL_POINTS    ' + str(nnodes) + '\n'
        data[4] = 'ELEMENTS        ' + str(nelems) + '\n'
        data[5] = 'BOUNDARIES      ' + str(nbound) + '\n'
        data[30] = '    INCLUDE  3D_' + str(self.lc) + '.geo.dat\n'
        data[40] = '    include 3D_' + str(self.lc) + '.fix.bou\n'
        with open('3D.dom.dat', 'w') as f:
            f.writelines(data)
        with open('3D.cell.dat', 'w') as f:
            for j in range(0, nnodes):
                f.write(str(j + 1) + '\t1\n')
        with open('3D.fibre.dat', 'w') as f:
            for j in range(0, nnodes):
                f.write(str(j + 1) + ' 0 0 1\n')

        with open('3D.sheet.dat', 'w') as f:
            for j in range(0, nnodes):
                f.write(str(j + 1) + ' 1 0 0\n')

        with open('3D.normal.dat', 'w') as f:
            for j in range(0, nnodes):
                f.write(str(j + 1) + ' 0 1 0\n')

        with open('3D.ker.dat', 'r') as f:
            data = f.readlines()
        data[8] = '$DIVISION=0\n'

        with open('3D.ker.dat', 'w') as f:
            f.writelines(data)


    def objective_function(self, conductivities):
        # Update speeds
        with open('3D.exm.dat', 'r') as f:
            data = f.readlines()
        data[27] = '                IN_DIFFUSION=  '+str(conductivities[0])+' '+str(conductivities[1])+' '+str(conductivities[2])+'\n'
        with open('3D.exm.dat', 'w') as f:
            f.writelines(data)

        # Run Alya
        print ('Running Alya simulation with conductivities: '+str(conductivities))
        st = time.time()
        os.system('rm alyalog')
        alya_executable = self.alya_executable_path
        os.system('mpirun -n 4 ' + alya_executable +' 3D >> alyalog')
        print ('Time taken: '+str((time.time()-st)/60)+' minutes.')

        # Evaluate activation time map
        print ('Evaluating activation time map')
        os.system('mpirun -n 4 alya2csvensight_mpi.py')
        if not os.path.exists('output_lc'+str(self.lc)):
            os.mkdir('output_lc'+str(self.lc))
        os.system('mv 3D-* output_lc'+str(self.lc))
        os.system('cp 3D.post.* output_lc'+str(self.lc))
        os.system('cp alya2csvensight_mpi.py output_lc'+str(self.lc))
        os.chdir('output_lc'+str(self.lc))
        #os.system('/data/Alya/Alya_multiple_BZRZ_models/Utils/user/alya2pos/alya2pos.x 3D >> ensightlog')
        os.system('rm maplog')
        os.system('mpirun -n 3 python alya2csvensight_mpi.py >> maplog')
        os.chdir('../')

        # Evaluate CVs
        print('Evaluating CVs')
        CV = np.loatxt()
        CV = self.get_AT(self.lc, self.local_node_list,self.local_node_list_coords, self.local_weight_list)

        # Target speed
        fCVt = self.target_cvs[0] # 65. # cm/s From Taggart et al. (2000)
        sCVt = self.target_cvs[1] #44. # cm/s  Inferred
        nCVt = self.target_cvs[2] # 48. # cm/s From Taggart et al. (2000)

        # Speed error:
        e = abs(CV[0] - fCVt ) + abs(CV[1] - sCVt) + abs(CV[2] - nCVt)
        print('CVs: ', CV)
        print('Error: ', e)
        return e


    def dist_to_plane(self, A, B, C, D):
        normal = np.cross((A - B), (C - B))
        dist = np.dot((D - B), normal)  # This can be negative, in which case the point is outside of the element.
        return dist

    def tet_local_coords(self, point, tet_nodes):
        local_coord = np.zeros((4))
        local_coord[3] = self.dist_to_plane(tet_nodes[2, :], tet_nodes[1, :], tet_nodes[0, :], point) / dist_to_plane(
            tet_nodes[2, :], tet_nodes[1, :], tet_nodes[0, :], tet_nodes[3, :])
        local_coord[2] = self.dist_to_plane(tet_nodes[0, :], tet_nodes[1, :], tet_nodes[3, :], point) / dist_to_plane(
            tet_nodes[0, :], tet_nodes[1, :], tet_nodes[3, :], tet_nodes[2, :])
        local_coord[1] = self.dist_to_plane(tet_nodes[0, :], tet_nodes[2, :], tet_nodes[3, :], point) / dist_to_plane(
            tet_nodes[0, :], tet_nodes[2, :], tet_nodes[3, :], tet_nodes[1, :])
        local_coord[0] = self.dist_to_plane(tet_nodes[1, :], tet_nodes[2, :], tet_nodes[3, :], point) / dist_to_plane(
            tet_nodes[1, :], tet_nodes[2, :], tet_nodes[3, :], tet_nodes[0, :])
        return local_coord  # They always add to one, but some values can be negative, if point lies outside of element, which is fine.

    def is_in_elem(self, point, tet_nodes):
        local_coord = self.tet_local_coords(point, tet_nodes)
        return np.amin(local_coord) > 0

    def get_AT(self, lc, local_node_list, local_node_list_coords, local_weight_list):
        activ = np.load(self.calibration_dir + self.name + '_lat.csv')
        CV = np.zeros((3))
        counter = 0
        for i in [0, 1, 2]:
            activ_aux = np.zeros((2))
            nodes_aux = np.zeros((2, 3))
            for j in range(0, activ_aux.shape[0]):
                activ_aux[j] = np.dot(activ[local_node_list[counter, :]],
                                      local_weight_list[counter, :])  # Weighted average
                nodes_aux[j, :] = np.dot(local_node_list_coords[counter, :, :].transpose(),
                                         local_weight_list[counter, :].transpose())
                counter = counter + 1
            CV[i] = np.linalg.norm(nodes_aux[1, :] - nodes_aux[0, :]) / (activ_aux[1] - activ_aux[0])
        return CV

    def get_element_probes(self, lc, width, cube_length):
        centre = np.array([0.5 * cube_length, 0.5 * cube_length, 0.5 * cube_length])
        nodes = np.load(self.calibration_dir + self.name + '_point_coordinates.npy')
        elems = np.load('output_lc' + str(lc) + '/ELEMS.npy')

        # Compute edge lenghts:
        elems = elems - 1
        elements_1 = np.sort(elems[:, [0, 2]], 1)
        elements_2 = np.sort(elems[:, [0, 2]], 1)
        elements_3 = np.sort(elems[:, [0, 3]], 1)
        elements_4 = np.sort(elems[:, [1, 2]], 1)
        elements_5 = np.sort(elems[:, [1, 3]], 1)
        elements_6 = np.sort(elems[:, [2, 3]], 1)
        edges = np.unique(
            np.concatenate((elements_1, elements_2, elements_3, elements_4, elements_5, elements_6), axis=0), axis=0)
        edge_length = np.linalg.norm(nodes[edges[:, 0]] - nodes[edges[:, 1]], axis=1)
        print('Mean edge length is: {:.2f} +- {:.2f} micron'.format(np.mean(edge_length) * 10000,
                                                                    np.std(edge_length) * 10000))

        probes = np.zeros((2, 3))
        local_node_list = np.zeros((6, 4), dtype='int')
        local_node_list_coords = np.zeros((6, 4, 3))
        local_weight_list = np.zeros((6, 4))
        CV = np.zeros((3))
        counter = 0
        for i in [0, 1, 2]:
            probes[0, :] = centre
            probes[1, :] = centre
            probes[0, i] = centre[i] - width / 2.
            probes[1, i] = centre[i] + width / 2.
            print(probes)
            nodes_aux = np.zeros((2, 3))
            for j in range(0, probes.shape[0]):
                in_elem = np.zeros(elems.shape[0], dtype='bool')
                for e in range(0, elems.shape[0]):
                    local_nodes = nodes[elems[e, :], :]
                    in_elem[e] = self.is_in_elem(probes[j, :], local_nodes)
                print(in_elem.nonzero())
                local_elem = in_elem.nonzero()[0][0]
                local_node_list[counter, :] = elems[local_elem, :]
                local_weight_list[counter, :] = self.tet_local_coords(probes[j, :], nodes[local_node_list[counter, :], :])
                local_node_list_coords[counter, :, :] = nodes[local_node_list[counter, :], :]
                counter = counter + 1
        return local_node_list, local_node_list_coords, local_weight_list





