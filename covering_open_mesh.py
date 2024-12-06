import matplotlib
matplotlib.use('tkagg')
from matplotlib import pyplot as plt
import numpy as np

def vector_angle(x, start_angle, flat_x):
    if x <= flat_x:
        return start_angle/(flat_x)**2 * (x-flat_x)**2 # This is a quadratic function that begins at start angle and goes to zero at end_x
    else:
        return 0

twod_profile = False
threed_profile = True
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
    x = np.linspace(-1, 1, 100)
    y = np.sqrt(1 - x**2)
    y_append = - np.sqrt(1 - x**2)
    x_append = x
    final_x = np.concatenate((x, x_append))
    final_y = np.concatenate((y, y_append))
    final_z = np.zeros(final_x.shape[0])
    edge_nodes = np.vstack((final_x, final_y, final_z)).transpose()
    print(edge_nodes)
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(edge_nodes[:, 0], edge_nodes[:, 1], edge_nodes[:, 2], marker='*', c='b')
    plt.show()

    # Create curves from each edge node to the middle
    start_vectors = np.zeros((edge_nodes.shape[0], 3))
    start_vectors[:, 0] = 0
    start_vectors[:, 1] = 0
    start_vectors[:, 2] = 1
    longitudinal_vector = [0, 0, 1]
    middle_coord = np.array([0, 0, 0])
    for edge_node, start_vector in zip(edge_nodes, start_vectors):
        # Create local axes
        x_vector = middle_coord - edge_node
        y_axis = np.cross(x_vector/np.linalg.norm(x_vector), longitudinal_vector)
        y_axis = y_axis/np.linalg.norm(y_axis)
        x_axis = np.cross(longitudinal_vector, y_axis)
        x_axis = x_vector/np.linalg.norm(x_axis)
        true_x_vector = np.dot(x_vector, x_axis) * x_axis
        dx = 0.01 * true_x_vector
        n_points = np.linalg.norm(true_x_vector)/ dx
        flat_x = 0.8 * np.linalg.norm(x_vector)

        end_vector = true_x_vector / np.linalg.norm(true_x_vector)
        start_angle = np.arccos(np.dot(start_vector, end_vector))
        current_coord = edge_node
        for i in range(1, n_points):
            current_coord = current_coord + dx
            angle = vector_angle(x[i], start_angle, flat_x)
            basis_transformation = np.vstack((x_axis, y_axis, longitudinal_vector))
            R = np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
            new_vector = np.matmul(R, np.matmul(basis_transformation,end_vector))
            plt.quiver(current_coord[0], current_coord[1], current_coord[2],  new_vector[0], new_vector[1])
            curve[i, :] = curve[i - 1, :] + new_vector * dx / new_vector[0]
            prev_vector = new_vector
