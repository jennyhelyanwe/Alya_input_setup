import os, sys
import numpy as np


def evaluate_mesh_characteristics(geometry):
    tetrahedrons = geometry.tetrahedrons
    number_of_nodes = geometry.number_of_nodes
    # edges = []
    # for i in range(0, geometry.number_of_elements):
    #     edges.append([tetrahedrons[i, 0], tetrahedrons[i, 1]])
    #     edges.append([tetrahedrons[i, 1], tetrahedrons[i, 2]])
    #     edges.append([tetrahedrons[i, 2], tetrahedrons[i, 3]])
    #     edges.append([tetrahedrons[i, 3], tetrahedrons[i, 0]])
    #     edges.append([tetrahedrons[i, 1], tetrahedrons[i, 3]])
    #     edges.append([tetrahedrons[i, 0], tetrahedrons[i, 2]])
    # edges = np.unique(np.sort(edges, axis=1), axis=0)
    edges = geometry.edges
    unfolded_edges = np.concatenate((edges, np.flip(edges, axis=1))).astype(int)
    aux = [[] for i in range(0, number_of_nodes, 1)]
    for i in range(0, len(unfolded_edges)):
        aux[unfolded_edges[i, 0]].append(unfolded_edges[i, 1])
    neighbours = [np.array(n) for n in aux]  # Node numbers starting 0
    return neighbours, edges, unfolded_edges



def evaluate_dijkstra_endocardial_activation(number_of_nodes, number_of_faces, face_fields, root_node_locations):
    lv_endocardial_nodes = face_fields


def evaluate_endocardial_activation_map(activation_times, geometry):
    lv_activation_times = activation_times[geometry.dict['lvnodes']]
    rv_activation_times = activation_times[geometry.dict['rvnodes']]
    stimulus = np.vstack((np.concatenate((geometry.dict['lvnodes'], geometry.dict['rvnodes'])), np.concatenate((lv_activation_times, rv_activation_times))))
    print(stimulus.shape)
    quit()

def evaluate_celltype(number_of_nodes,uvc_transmural, endo_mid_divide, mid_epi_divide):
    print('Evaluating nodal cell type...')
    celltype = [1]*number_of_nodes
    for node_i in range(number_of_nodes):
        if uvc_transmural[node_i] <= endo_mid_divide:
            celltype[node_i] = 1 # Endocardial cell type in Alya.
        if uvc_transmural[node_i] > endo_mid_divide and uvc_transmural[node_i] < mid_epi_divide: # It is possible to have no mid-myocardial cells.
            celltype[node_i] = 2 # Midmyocardial cell type in Alya.
        if uvc_transmural[node_i] >= mid_epi_divide:
            celltype[node_i] = 3 # Epicardial cell type in Alya.
    return np.array(celltype)

def evaluate_ab_Gks_scaling(number_of_nodes, uvc_longitudinal, max_sf, min_sf):
    gks_scaling = [0] * number_of_nodes
    for node_i in range(number_of_nodes):
        z_scale = uvc_longitudinal[node_i] * 2.0 - 1.0
        gks_scaling[node_i] = pow(0.2, z_scale)
        if (uvc_longitudinal[node_i] < 0):
            gks_scaling[node_i] = 0
    return np.array(gks_scaling)

def evaluate_materials(number_of_elements, lvrv_element_field):
    materials = [0] * number_of_elements
    for elem_i in range(number_of_elements):
        if lvrv_element_field[elem_i] == 3:
            materials[elem_i] = 2
        else:
            materials[elem_i] = 1
    return np.array(materials)


# def evaluate_hybrid_rodero_fibres(number_of_nodes, number_of_elements, uvc_transmural, uvc_longitudinal, fibres):
#

def normalise_vector(vector):
    m = np.linalg.norm(vector)
    if m > 0:
        return vector / m
    else:
        return vector
def evaluate_hybrid_rodero_fibres(geometry, node_fields, element_fields, neighbours):
    transmural_vector = node_fields.dict['transmural-vector']
    fibres_nodes = node_fields.dict['fibres']
    number_of_nodes = geometry.number_of_nodes
    ortho = np.zeros([number_of_nodes, 9])
    for i in range(0, number_of_nodes):
        f = fibres_nodes[i, :]
        f = f / np.linalg.norm(f)
        s = transmural_vector[i, :]
        n = np.cross(f, s)
        n = normalise_vector(n)
        ortho[i, :] = list(f) + list(s) + list(n)
    plug_fibres(geometry=geometry, element_fields=element_fields, neighbours=neighbours, ortho=ortho)
    # Sanity check:
    zero_fvector_nodes = np.where(~ortho[:, 0:3].any(axis=1))[0]  # Vectors all zero
    zero_svector_nodes = np.where(~ortho[:, 3:6].any(axis=1))[0]
    zero_nvector_nodes = np.where(~ortho[:, 6:9].any(axis=1))[0]
    nan_vector_nodes = np.where(np.isnan(ortho).any(axis=1))[0]
    print('Ortho calculation sanity check: ')
    print('Number of zero fibre vectors: ', len(zero_fvector_nodes))
    print('Number of zero sheet vectors: ', len(zero_fvector_nodes))
    print('Number of zero normal vectors: ', len(zero_fvector_nodes))
    print('Number of vectors with NaNs: ', len(nan_vector_nodes))
    return ortho

def plug_fibres(geometry, element_fields, neighbours, ortho ):
    print('Evaluating valvular plug fibre, sheet, and normal vectors...')
    # Identify all plug elements
    lvrvbase_elems = element_fields.dict['tv-element']
    elems = geometry.tetrahedrons
    nodes = geometry.nodes_xyz
    nnodes = geometry.number_of_nodes
    for k in range(0, 4):
        # Find all elements and nodes in current valve plug
        plug_elems = np.nonzero(lvrvbase_elems == k + 7)[0]  # Plug element is labelled as 7 to 10.
        plug_nodes = np.unique(elems[plug_elems])
        plug_nodes_coords = nodes[plug_nodes - 1, :]

        non_plug_elems = elems[np.nonzero(lvrvbase_elems <= 2)[0]]
        non_plug_nodes = np.unique(non_plug_elems)
        border_elems = np.nonzero(np.any(np.isin(non_plug_elems, plug_nodes), axis=1))
        border_nodes = np.unique(non_plug_elems[border_elems, :])
        ring_nodes = np.array(list(set(border_nodes).difference(set(plug_nodes))), dtype=int)
        ring_nodes_coords = nodes[ring_nodes - 1, :]
        non_plug_nodes = np.array(list(set(non_plug_nodes).difference(set(plug_nodes))), dtype=int)

        # Ordering all plug nodes according to how far they are from the nearest ventricular node
        dists = np.zeros(len(plug_nodes))
        for i in range(0, len(plug_nodes)):
            dists[i] = np.amin(np.linalg.norm(ring_nodes_coords - plug_nodes_coords[i, :], axis=1))
        sorted_plug_nodes = plug_nodes[np.argsort(dists)]
        mask = np.zeros(nnodes)
        mask[non_plug_nodes - 1] = 1
        for i in range(0, len(sorted_plug_nodes)):
            all_closest_nodes = neighbours[sorted_plug_nodes[i] - 1]
            closest_nodes = all_closest_nodes[np.where(mask[all_closest_nodes] > 0)]
            if (len(closest_nodes) == 0):
                print('No viable neighbours found for node number: ', sorted_plug_nodes[i],
                      ' plug number: ' + str(k + 1))
                print('finding next closest nodes...')
                temp = []
                for i in range(0, len(all_closest_nodes)):
                    temp = temp + list(neighbours[all_closest_nodes[i]])
                all_closest_nodes = np.unique(np.array(temp))
                closest_nodes = all_closest_nodes[np.where(mask[all_closest_nodes] > 0)]
            weight = 0.4
            fibres = np.concatenate(
                (ortho[closest_nodes, :], weight * np.array([[1, 0, 0, 0, 1, 0, 0, 0, 1]], dtype=float)))
            new_fibres = np.mean(fibres, axis=0)
            for j in [0, 3, 6]:
                new_fibres[j:j + 3] = normalise_vector(new_fibres[j:j + 3])
            ortho[sorted_plug_nodes[i] - 1, :] = new_fibres
            mask[sorted_plug_nodes[i] - 1] = 1