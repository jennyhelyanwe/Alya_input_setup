import matplotlib
#matplotlib.use('tkagg')
from matplotlib import pyplot as plt
import numpy as np
import vtk
from vtk.util import numpy_support as VN
import pandas

def vector_angle(x, start_angle, flat_x):
    if x <= flat_x:
        return start_angle/(flat_x)**2 * (x-flat_x)**2 # This is a quadratic function that begins at start angle and goes to zero at end_x
    else:
        return 0

def plot_3D_vector(ax, coord, vector, color):
    ax.quiver(coord[0], coord[1], coord[2], vector[0], vector[1], vector[2], color=color)


def generate_cap_surface_data(edge_nodes, start_vectors, longitudinal_vector, lid_rise, flat_percentage, resolution_along_radius):
    # Create curves from each edge node to the middle
    middle_coord = np.mean(edge_nodes, axis=0)
    middle_lid_coord = middle_coord + longitudinal_vector*lid_rise
    # plot_3D_vector(ax, middle_lid_coord, longitudinal_vector, 'k')
    data = np.array([[0, 0, 0]])
    for edge_node, start_vector in zip(edge_nodes, start_vectors):
        # Create local axes
        x_vector = middle_coord - edge_node
        y_axis = np.cross(x_vector/np.linalg.norm(x_vector), longitudinal_vector)
        y_axis = y_axis/np.linalg.norm(y_axis)
        x_axis = np.cross(longitudinal_vector, y_axis)
        x_axis = x_axis/np.linalg.norm(x_axis)
        dx = resolution_along_radius * (x_vector)/np.linalg.norm(x_vector)
        # plot_3D_vector(ax, edge_node, dx, 'm')
        # plot_3D_vector(ax, edge_node, x_axis, 'g')
        # plot_3D_vector(ax, edge_node, y_axis, 'r')
        # plot_3D_vector(ax, edge_node, longitudinal_vector, 'b')
        n_points = np.floor(np.linalg.norm(x_vector)/ resolution_along_radius).astype(int)
        curve = np.zeros((n_points, 3))
        curve[0, :] = edge_node
        flat_x = flat_percentage * np.linalg.norm(x_vector)
        end_vector = x_vector / np.linalg.norm(x_vector)
        start_angle = np.arccos(np.dot(start_vector, end_vector))
        # plot_3D_vector(ax, edge_node, start_vector, 'k')
        current_coord = edge_node
        for i in range(1, n_points):
            angle = vector_angle(i * resolution_along_radius, start_angle, flat_x)
            basis_transformation = np.vstack((x_axis, y_axis, longitudinal_vector)).transpose()
            R = np.array([[np.cos(angle), 0, -np.sin(angle)], [0, 1, 0], [np.sin(angle), 0, np.cos(angle)]])
            new_vector = np.matmul(R, np.matmul(np.linalg.inv(basis_transformation),end_vector))
            current_coord = current_coord + dx
            new_vector_global = np.matmul(basis_transformation, new_vector)
            new_vector_global = new_vector_global/np.linalg.norm(new_vector_global)
            # plot_3D_vector(ax, current_coord, new_vector_global, 'b')
            # plot_3D_vector(ax, current_coord, dx, 'm')
            # scaling = resolution_along_radius / (np.dot(new_vector_global, dx/np.linalg.norm(dx)))
            curve[i, :] = curve[i-1,:] + new_vector_global * resolution_along_radius
            # ax.scatter(curve[0, 0], curve[0, 1], curve[0, 2], marker='o', c='m')
            # ax.scatter(curve[i, 0], curve[i, 1], curve[i, 2], marker='o', c='g')
            # plt.show()
            # quit()
        # lid_scale = lid_rise/np.dot(longitudinal_vector, curve[-1,:] - middle_coord)
        # x_scale = np.linalg.norm(x_vector)/np.dot(x_axis, np.linalg.norm(curve[-1,:] - edge_node))
        lid_scale = lid_rise / (np.dot(longitudinal_vector, curve[-1,:] - edge_node))
        x_scale =  np.linalg.norm(x_vector)/np.dot(x_axis, curve[-1,:] - edge_node)
        for i in range(0, n_points):
            curve_x = np.dot((curve[i,:]-edge_node), x_axis)
            curve_y = np.dot((curve[i,:]-edge_node), y_axis)
            curve_z = np.dot((curve[i,:]-edge_node), longitudinal_vector)
            curve[i, :] = (curve_x* x_axis*x_scale + curve_y*y_axis + curve_z * longitudinal_vector * lid_scale) + edge_node
        curve = curve[1:, :]
        # ax.scatter(curve[:, 0], curve[:, 1], curve[:, 2], marker='*', c='m')
        # plt.show()
        # quit()
        data = np.concatenate((data,curve))
    data = data[1:,:]
    return data

twod_profile = False
threed_profile = False
rodero_05_lvendocardium = True
rodero_05_rvendocardium = True
# 2D generation of side profile
if twod_profile:

    end_x = 1
    flat_x = 0.8
    x = np.linspace(0, 1, 100)
    print(x)
    vector_angles = np.zeros(x.shape)


    start_angle = np.pi/2 # degrees
    for i in range(0, x.shape[0]):
        vector_angles[i] = vector_angle(x[i], start_angle, flat_x)

    # print(vector_angles)
    # plt.plot(x, vector_angles, '-')
    # plt.show()

    # Integrate curve
    start_vector = np.array([0, 1])
    end_vector = np.array([1,0])
    y_vector= np.array([0, 1])
    dx = x[1] - x[0]
    curve = np.zeros((x.shape[0], 2))
    curve[:, 0] = x[:]
    for i in range(1, x.shape[0]):
        angle = vector_angle(x[i], start_angle, flat_x)
        R = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        new_vector = np.matmul(R, end_vector)
        plt.quiver(x[i], 0, new_vector[0], new_vector[1])
        curve[i, :] = curve[i-1, :] + new_vector * dx/ new_vector[0]
        prev_vector = new_vector

    plt.plot(curve[:,0], curve[:,1], '-')
    plt.show()

if threed_profile:
    # 3D generation of surface data points
    # Generate the nodes on the edge
    x = np.linspace(-2, 2, 100)
    y = np.sqrt(4 - x**2)
    y_append = - np.sqrt(4 - x**2)
    x_append = x
    final_x = np.concatenate((x, x_append))
    final_y = np.concatenate((y, y_append))
    final_z = np.zeros(final_x.shape[0])
    edge_nodes = np.vstack((final_x, final_y, final_z)).transpose()
    start_vectors = np.zeros((edge_nodes.shape[0], 3))
    start_vectors[:, 0] = 0
    start_vectors[:, 1] = 0
    start_vectors[:, 2] = 1
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(edge_nodes[:, 0], edge_nodes[:, 1], edge_nodes[:, 2], marker='*', c='b')
    ax.set_zlim(-2, 2)
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    longitudinal_vector = np.array([0, 0, 1])
    cap_data = generate_cap_surface_data(edge_nodes=edge_nodes, start_vectors=start_vectors,
                                         longitudinal_vector=longitudinal_vector, lid_rise=0.2,
                                         flat_percentage=0.8,
                                         resolution_along_radius=0.05)
    ax.scatter(cap_data[:, 0], cap_data[:, 1], cap_data[:, 2], marker='*', c='m')
    plt.show()

if rodero_05_lvendocardium:
    vtk_names = ['DTI003/lv_endocardium.vtk',
                 'DTI003/rv_endocardium.vtk',
                 'DTI003/epicardium.vtk']
    edge_names = ['DTI003/lv_endo_edge_node_start_vectors.csv',
                  'DTI003/rv_endo_edge_node_start_vectors.csv',
                  'DTI003/epicardium_edge_node_start_vectors.csv']
    save_names = ['DTI003/lvendo_cap_and_surface.txt',
                  'DTI003/rvendo_cap_and_surface.txt',
                  'DTI003/epicardium_edge_node_start_vectors.txt']
    lid_rise_percentages = [0.15, 0.15, 0.2]
    # vtk_names = ['DTI003/rv_endocardium.vtk']
    # edge_names = ['DTI003/rv_endo_edge_node_start_vectors.csv']
    # save_names = ['DTI003/rvendo_cap_and_surface.txt']
    for vtk_name, edge_name, save_name, lid_rise_percentage in zip(vtk_names, edge_names, save_names, lid_rise_percentages):
        reader = vtk.vtkUnstructuredGridReader()
        reader.SetFileName(vtk_name)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()
        data = reader.GetOutput()
        nodes_xyz = VN.vtk_to_numpy(data.GetPoints().GetData())  # Convert from mm to cm
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection='3d')
        ax.scatter(nodes_xyz[:, 0], nodes_xyz[:, 1], nodes_xyz[:, 2], marker='*', c='b')
        node_ids = VN.vtk_to_numpy(data.GetPointData().GetArray('PointIds'))
        normal_vector = np.array([-0.6767995188650355, 0.5464585613268474, 0.49328029761652614])
        existing_lv_length = np.amax(np.abs(np.sum(normal_vector * nodes_xyz, axis=1, keepdims=True)))
        lid_rise = existing_lv_length * lid_rise_percentage
        data = pandas.read_csv(edge_name)
        edge_node_idx = data['PointIds'].values
        start_vectors = np.zeros((edge_node_idx.shape[0], 3))
        start_vectors[:, 0] = data['ScalarGradient:0'].values
        start_vectors[:, 1] = data['ScalarGradient:1'].values
        start_vectors[:, 2] = data['ScalarGradient:2'].values
        # edge_node_idx = np.loadtxt('DTI003/epi_edge_node_idx.csv', delimiter=',')[:, 0].astype(int)
        edge_nodes = np.zeros((edge_node_idx.shape[0], 3))
        for i, idx in enumerate(edge_node_idx):
            edge_nodes[i, :] = nodes_xyz[np.where(node_ids==idx)[0],:]
        # ax.scatter(edge_nodes[:, 0], edge_nodes[:, 1], edge_nodes[:, 2], marker='*', c='r')
        # start_vectors = np.zeros((edge_nodes.shape[0], 3))
        # start_vectors[:, 0] = normal_vector[0]
        # start_vectors[:, 1] = normal_vector[1]
        # start_vectors[:, 2] = normal_vector[2]
        cap_data = generate_cap_surface_data(edge_nodes=edge_nodes, start_vectors=start_vectors,
                                             longitudinal_vector=normal_vector, lid_rise=lid_rise,
                                             flat_percentage=0.9,resolution_along_radius=0.05)
        cap_and_surface_data = np.concatenate((cap_data, nodes_xyz), axis=0)
        ax.scatter(cap_data[:, 0], cap_data[:, 1], cap_data[:, 2], marker='*', c='m')
        plt.show()
        np.savetxt(save_name, cap_and_surface_data)
        reader = None

# Use MeshLab to generate the surface mesh from this point on
