import os, sys
import numpy as np


def evaluate_dijkstra_endocardial_activation(number_of_nodes, number_of_faces, face_fields, root_node_locations):

    lv_endocardial_nodes = face_fields


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
    return np.arary(celltype)

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
