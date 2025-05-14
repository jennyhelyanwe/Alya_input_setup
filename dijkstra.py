import numpy as np
import math
import numba
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
# Imported and adapted from https://github.com/juliacamps/Cardiac-Digital-Twin/blob/main/src/
# 1 May 2025


def get_nan_value():
    return -2147483648  # min value for int32 #np.array([np.nan]).astype(np.int32)[0]

class Dijkstra:
    def prepare_for_dijkstra(self, edge, node_xyz, sub_node_index=None):

        if not sub_node_index:
            sub_node_xyz = node_xyz
            sub_edge = edge.astype(int)
        else:
            sub_node_xyz = node_xyz[sub_node_index, :]
            node_xyz = None
            sub_edge = edge[np.all(np.isin(edge, sub_node_index), axis=1), :]
            edge = None
            aux_edge = sub_edge
            sub_edge[:, 0] = np.asarray([np.flatnonzero(sub_node_index == node_i)[0] for node_i in sub_edge[:, 0]]).astype(int)
            sub_edge[:, 1] = np.asarray([np.flatnonzero(sub_node_index == node_i)[0] for node_i in sub_edge[:, 1]]).astype(int)
        sub_edge_vec = sub_node_xyz[sub_edge[:, 0], :] - sub_node_xyz[sub_edge[:, 1], :]  # edge vectors
        sub_unfolded_edge = np.concatenate((sub_edge, np.flip(sub_edge, axis=1))).astype(int)
        aux = [[] for i in range(0, sub_node_xyz.shape[0], 1)]
        for i in range(0, len(sub_unfolded_edge), 1):
            aux[sub_unfolded_edge[i, 0]].append(i)
        sub_neighbour = [np.array(n, dtype=int) for n in aux]
        return sub_neighbour, sub_node_xyz, sub_unfolded_edge, sub_edge_vec


    def sort_dijkstra_by_distance(self, source_to_all_distance_mat, source_to_all_path_mat):
        for destination_index in range(source_to_all_path_mat.shape[0]):
            for source_i in range(source_to_all_path_mat.shape[1]):
                path_indexes = source_to_all_path_mat[destination_index, source_i, :]
                path_indexes = path_indexes[path_indexes != get_nan_value()]
                source_to_path_distance = source_to_all_distance_mat[path_indexes, source_i]
                sorted_by_distance_indexes = np.argsort(
                    source_to_path_distance)  # Sort nodes by distance to the reference
                source_to_all_path_mat[destination_index, source_i, :path_indexes.shape[0]] = path_indexes[
                    sorted_by_distance_indexes]
        return source_to_all_path_mat


    def sorted_dijkstra(self, source_indexes, seed_times, dijkstra_nodes_xyz, dijkstra_unfoldedEdges, dijkstra_edgeVEC,
                        dijkstra_neighbours,
                        approx_dijkstra_max_path_len):
        # Djikstra
        print('Solving Dijkstra...')
        source_to_all_distance_mat, source_to_all_path_mat = self.dijkstra(
            source_id_list=np.asarray(source_indexes, dtype=int),
            dijkstra_nodes_xyz=dijkstra_nodes_xyz,
            dijkstra_unfoldedEdges=dijkstra_unfoldedEdges,
            dijkstra_edgeVEC=dijkstra_edgeVEC,
            dijkstra_neighbours=dijkstra_neighbours,
            approx_max_path_len=approx_dijkstra_max_path_len)
        # Sort path dijkstra results by distance
        source_to_all_path_mat = self.sort_dijkstra_by_distance(source_to_all_distance_mat=source_to_all_distance_mat,
                                                           source_to_all_path_mat=source_to_all_path_mat)
        distance_mat_with_offset = source_to_all_distance_mat + seed_times
        min_distance = np.min(distance_mat_with_offset, axis=1)
        return min_distance


    def insert_sorted(self, aList, newV):
        # Function to insert element
        ini_index = 0
        end_index = len(aList)
        index = int((end_index - ini_index) / 2)
        for i in range(0, len(aList), 1):
            if newV[1] < aList[index][1]:
                if end_index - ini_index <= 1 or index + 1 == end_index:
                    index = index + 1
                    break
                else:
                    ini_index = index + 1
                    index = int(index + (end_index - ini_index) / 2 + 1)
            elif newV[1] > aList[index][1]:
                if end_index - ini_index <= 1 or index == ini_index:
                    index = index  # Place before the current position
                    break
                else:
                    end_index = index
                    index = int(index - (end_index - ini_index) / 2)
            else:
                index = ini_index
                break
        aList.insert(index, newV)

    def eikonal_part1(self, sub_node_coordinates, source_indexes, seed_times, sub_edge_indexes, sub_edgeVEC, sub_neighbour, sub_unfoldedEdges):
        root_node_indexes = source_indexes
        root_node_times = seed_times
        # TODO REIMPLEMENT THE ABOVE SECTION USING FUNCTIONS FROM conduction_system.py

        ## ISOTROPIC REGIONS - without fibre orientation
        # Compute the cost of all endocardial edges
        navigation_costs = np.empty((sub_edge_indexes.shape[0]))
        # print('2')
        for index in range(sub_edge_indexes.shape[0]):
            # Cost for the propagation in the endocardium
            # if self.geometry.is_endocardial[index]:
            navigation_costs[index] = math.sqrt(
                np.dot(sub_edgeVEC[index, :], sub_edgeVEC[index, :]))
        # Build adjacentcy costs for current activation parameter values
        adjacent_cost = numba.typed.List()
        for i in range(0, sub_node_coordinates.shape[0], 1):
            not_nan_neighbours = sub_neighbour[i][sub_neighbour[i] != get_nan_value()]
            adjacent_cost.append(np.concatenate((sub_unfoldedEdges[not_nan_neighbours][:, 1:2],
                                                 np.expand_dims(navigation_costs[
                                                                    not_nan_neighbours % sub_edge_indexes.shape[0]],
                                                                -1)), axis=1))
        return adjacent_cost, root_node_indexes, root_node_times


    def iso_eikonal(self, sub_node_coordinates, source_indexes, seed_times, sub_edge_indexes, sub_edgeVEC, sub_neighbour, sub_unfoldedEdges):
        # Initialise variables
        predicted_lat = np.zeros((sub_node_coordinates.shape[0],), np.float64)
        visited_nodes = np.zeros((sub_node_coordinates.shape[0],), dtype=np.bool_)
        # Root nodes will be activated at time ==  self.root_node_times
        temp_times = np.zeros((sub_node_coordinates.shape[0],),
                              np.float64) + 1e6  # Initialise times to largely impossible values
        adjacent_cost, eikonal_root_nodes, eikonal_root_lat = self.eikonal_part1(sub_node_coordinates, source_indexes,
                                                                              seed_times, sub_edge_indexes, sub_edgeVEC,
                                                                              sub_neighbour, sub_unfoldedEdges=sub_unfoldedEdges)
        temp_times[eikonal_root_nodes] = eikonal_root_lat
        time_sorting = np.argsort(eikonal_root_lat)
        print(time_sorting)
        eikonal_root_nodes = eikonal_root_nodes[time_sorting]
        eikonal_root_lat = eikonal_root_lat[time_sorting]
        eikonal_root_lat = eikonal_root_lat - eikonal_root_lat[0]
        cumm_cost = eikonal_root_lat[0]
        print('cumm_cost ', cumm_cost)
        initial_root_nodes_indexes = eikonal_root_lat <= cumm_cost
        initial_rootNodes = eikonal_root_nodes[initial_root_nodes_indexes]
        initial_rootActivationTimes = eikonal_root_lat[initial_root_nodes_indexes]
        later_rootNodes = eikonal_root_nodes[np.logical_not(initial_root_nodes_indexes)]
        later_rootActivationTimes = eikonal_root_lat[np.logical_not(initial_root_nodes_indexes)]
        ## Run the code for the root nodes
        visited_nodes[initial_rootNodes] = True  # Not simultaneous activation anymore
        predicted_lat[
            initial_rootNodes] = initial_rootActivationTimes  # Not simultaneous activation anymore
        print('initial_rootActivationTimes ', initial_rootActivationTimes)
        next_nodes = (np.vstack([adjacent_cost[initial_rootNodes[rootNode_i]]
                                 + np.array([0, initial_rootActivationTimes[rootNode_i]]) for rootNode_i in
                                 range(initial_rootNodes.shape[
                                           0])])).tolist()  # Not simultaneous activation anymore
        for rootNode_i in range(later_rootNodes.shape[0]):
            next_nodes.append(np.array([later_rootNodes[rootNode_i], later_rootActivationTimes[rootNode_i]]))

        activeNode_i = initial_rootNodes[0]  # eikonal_root_nodes[0]
        sortSecond = lambda x: x[1]
        next_nodes.sort(key=sortSecond, reverse=True)

        print('next_nodes ', next_nodes)
        while visited_nodes[activeNode_i]:
            nextEdge = next_nodes.pop()
            activeNode_i = int(nextEdge[0])
        cumm_cost = nextEdge[1]
        print('cumm_cost ', cumm_cost)
        if next_nodes:  # Check if the list is empty, which can happen while everything being Ok
            temp_times[np.array(next_nodes, dtype=int)[:, 0]] = np.array(next_nodes)[:,
                                                                1]  # 04/10/2022 Why is this happening?

        ## Run the whole algorithm
        # for i in range(0, node_coordinates.shape[0] - np.sum(visited_nodes)+1, 1):
        for i in range(0, sub_node_coordinates.shape[0] - np.sum(visited_nodes), 1):
            # if activeNode_i == 40:
            #   print('SUPER HAPPY ', cumm_cost)
            # if not visited_nodes[activeNode_i]: # This is not necessary if the loop is done well.
            visited_nodes[activeNode_i] = True
            predicted_lat[
                activeNode_i] = cumm_cost  # Instead of using cumCost, I could use the actual time cost for each node
            adjacents = (adjacent_cost[activeNode_i] + np.array(
                [0, cumm_cost])).tolist()  # Instead of using cumCost, I could use the actual time cost for each node
            # If I use the actual costs, I only have to do it the first time and then it will just propagate, I will have to use decimals though, so no more uint type arrays.
            for adjacent_i in range(0, len(adjacents), 1):
                if (not visited_nodes[int(adjacents[adjacent_i][0])]
                        and (temp_times[int(adjacents[adjacent_i][0])] >
                             adjacents[adjacent_i][1])):
                    self.insert_sorted(next_nodes, adjacents[adjacent_i])
                    temp_times[int(adjacents[adjacent_i][0])] = adjacents[adjacent_i][1]
            while visited_nodes[activeNode_i] and len(next_nodes) > 0:
                nextEdge = next_nodes.pop()
                activeNode_i = int(nextEdge[0])
            cumm_cost = nextEdge[1]

        # Clean Memory
        adjacent_cost = None  # Clear Mem
        visited_nodes = None  # Clear Mem
        tempTimes = None  # Clear Mem
        next_nodes = None  # Clear Mem
        tempVisited = None  # Clear Mem
        navigationCosts = None  # Clear Mem
        return predicted_lat

    def dijkstra(self, source_id_list, dijkstra_nodes_xyz, dijkstra_unfoldedEdges, dijkstra_edgeVEC, dijkstra_neighbours,
                 approx_max_path_len):
        """approx_max_path_len: this parameter is not a hard constraint, but an approximate value to speed up the method.
        In practice, the method will increase its own approx_max_path_len value until it has explored the whole range of
        possibilities."""
        distances = np.zeros((dijkstra_nodes_xyz.shape[0], source_id_list.shape[0])).astype(float)
        paths = np.full((dijkstra_nodes_xyz.shape[0], source_id_list.shape[0], approx_max_path_len), get_nan_value(),
                        np.int32)
        for source_id_index in range(source_id_list.shape[0]):
            source_id = source_id_list[source_id_index]
            distances_per_source = np.zeros((dijkstra_nodes_xyz.shape[0]))
            previous_node_indexes_temp = np.zeros((dijkstra_nodes_xyz.shape[0])).astype(
                int)  # 2022/01/10 This object starts with a wrong solution and it makes it better and better until in the end the
            # solution is the correct one. Idea by me and Jenny :-)
            visitedNodes = np.zeros((dijkstra_nodes_xyz.shape[0])).astype(bool)

            # Compute the cost of all endocardial edges
            navigationCosts = np.zeros((int(dijkstra_unfoldedEdges.shape[0] / 2)))
            for index in range(0, navigationCosts.shape[0]):
                # Cost for the propagation in the endocardium
                navigationCosts[index] = math.sqrt(np.dot(dijkstra_edgeVEC[index, :], dijkstra_edgeVEC[index, :]))

            # Build adjacentcy costs
            # TODO remove this:
            # if True:
            # for i in range(0, dijkstra_nodes_xyz.shape[0], 1):
            #     print('dijkstra_neighbours ', dijkstra_neighbours)
            adjacentCost = [np.concatenate((dijkstra_unfoldedEdges[dijkstra_neighbours[i]][:, 1][:, np.newaxis],
                                            navigationCosts[dijkstra_neighbours[i] % navigationCosts.shape[0]][:,
                                            np.newaxis]), axis=1) for i in range(0, dijkstra_nodes_xyz.shape[0], 1)]

            cummCost = 0.  # Distance from a node to itself is zero
            tempDists = np.zeros((dijkstra_nodes_xyz.shape[0],), float) + 1000

            ## Run the code for the root nodes
            visitedNodes[source_id] = True
            distances_per_source[source_id] = cummCost
            nextNodes = (adjacentCost[source_id] + np.array([0, cummCost])).tolist()
            activeNode_i = source_id
            sortSecond = lambda x: x[1]
            nextNodes.sort(key=sortSecond, reverse=True)
            previous_node_indexes_temp[activeNode_i] = activeNode_i  # 2022/01/10
            for nextEdge_aux in nextNodes:  # 2022/01/10
                previous_node_indexes_temp[int(nextEdge_aux[0])] = activeNode_i  # 2022/01/10
            while visitedNodes[activeNode_i] and len(nextNodes) > 0:
                nextEdge = nextNodes.pop()
                activeNode_i = int(nextEdge[0])
            cummCost = nextEdge[1]
            if nextNodes:  # Check if the list is empty, which can happen while everything being Ok
                tempDists[(np.array(nextNodes)[:, 0]).astype(int)] = np.array(nextNodes)[:, 1]

            ## Run the whole algorithm
            for i in range(distances_per_source.shape[0]):
                visitedNodes[activeNode_i] = True
                distances_per_source[activeNode_i] = cummCost
                adjacents = (adjacentCost[activeNode_i] + np.array([0, cummCost])).tolist()
                for adjacent_i in range(0, len(adjacents), 1):
                    if (not visitedNodes[int(adjacents[adjacent_i][0])] and (
                            tempDists[int(adjacents[adjacent_i][0])] > adjacents[adjacent_i][1])):
                        self.insert_sorted(nextNodes, adjacents[adjacent_i])
                        tempDists[int(adjacents[adjacent_i][0])] = adjacents[adjacent_i][1]
                        previous_node_indexes_temp[int(adjacents[adjacent_i][0])] = activeNode_i  # 2022/01/10
                while visitedNodes[activeNode_i] and len(nextNodes) > 0:
                    nextEdge = nextNodes.pop()
                    activeNode_i = int(nextEdge[0])
                cummCost = nextEdge[1]

            distances[:, source_id_index] = distances_per_source
            for dijkstra_node_id in range(0, dijkstra_nodes_xyz.shape[0], 1):  # 2022/01/10
                path_per_source = np.full((dijkstra_nodes_xyz.shape[0]), get_nan_value(), np.int32)  # 2022/01/10
                path_node_id = dijkstra_node_id  # 2022/01/10
                path_node_id_iter = 0  # 2022/01/10
                path_per_source[path_node_id_iter] = path_node_id  # 2022/01/14
                path_node_id_iter = path_node_id_iter + 1  # 2022/01/14
                while path_node_id != source_id:  # 2022/01/10
                    path_node_id = previous_node_indexes_temp[path_node_id]  # 2022/01/10
                    path_per_source[path_node_id_iter] = path_node_id  # 2022/01/10
                    path_node_id_iter = path_node_id_iter + 1  # 2022/01/10
                # If the path is longer than the current size of the matrix, make the matrix a little bigger and continue
                if path_node_id_iter + 1 > approx_max_path_len:  # 2022/01/11
                    paths_aux = np.full((dijkstra_nodes_xyz.shape[0], source_id_list.shape[0], path_node_id_iter + 10),
                                        get_nan_value(), np.int32)  # 2022/01/11
                    paths_aux[:, :, :approx_max_path_len] = paths  # 2022/01/11
                    paths = paths_aux  # 2022/01/11
                    approx_max_path_len = path_node_id_iter + 10  # 2022/01/11
                paths[dijkstra_node_id, source_id_index, :] = path_per_source[:approx_max_path_len]  # 2022/01/10
        return distances, paths  # 2022/01/

    def plot_geometry(self, xyz, tetrahedrons, lat_simulation):

        def list_faces(t):
            t.sort(axis=1)
            n_t, m_t = t.shape
            f = np.empty((4 * n_t, 3), dtype=int)
            i = 0
            for j in range(4):
                f[i:i + n_t, 0:j] = t[:, 0:j]
                f[i:i + n_t, j:3] = t[:, j + 1:4]
                i = i + n_t
            return f

        def extract_unique_triangles(t):
            _, indxs, count = np.unique(t, axis=0, return_index=True, return_counts=True)
            return t[indxs[count == 1]]

        def extract_surface(t):
            f = list_faces(t)
            f = extract_unique_triangles(f)
            return f

        surf = extract_surface(tetrahedrons)
        fig = go.Figure(data=[
            go.Mesh3d(
                x=xyz[:, 0],
                y=xyz[:, 1],
                z=xyz[:, 2],
                colorbar=dict(title=dict(text='z')),
                colorscale=px.colors.sequential.Viridis,
                # Intensity of each vertex, which will be interpolated and color-coded
                intensity=lat_simulation,
                # i, j and k give the vertices of triangles
                # here we represent the 4 triangles of the tetrahedron surface
                i=surf[:, 0],
                j=surf[:, 1],
                k=surf[:, 2],
                name='y',
                showscale=True
            )
        ])
        fig.update_layout(width=1000,
                          height=1000)

        fig.show()