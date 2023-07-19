import numpy as np
from meshstructure import MeshStructure
import os
import pymp, multiprocessing
from matplotlib import pyplot as plt

output_dir = '/users/jenang/ForAdrienJune2023/'
mesh_number = '05'
workdir = os.getcwd()
simulation_name = 'rodero_' + mesh_number + '_fine'
if 'icei' in workdir:
    system = 'jureca'
else:
    system = 'heart'
if system == 'jureca':
    meta_data_dir = '/p/project/icei-prace-2022-0003/wang1/Alya_pipeline/meta_data/'
elif system == 'heart':
    meta_data_dir = '/data/Personalisation_projects/meta_data/'
geometric_data_dir = meta_data_dir + 'geometric_data/rodero_'+mesh_number+'/rodero_'+mesh_number+'_fine/'
mesh = MeshStructure(name=simulation_name, geometric_data_dir=geometric_data_dir, verbose=False)

f_invariant_projection = pymp.shared.array(mesh.node_fields.dict['fibre'].shape)
s_invariant_projection = pymp.shared.array(mesh.node_fields.dict['fibre'].shape)
n_invariant_projection = pymp.shared.array(mesh.node_fields.dict['fibre'].shape)
transmural_vector = pymp.shared.array(mesh.node_fields.dict['transmural-vector'].shape)
transmural_vector[:,:] = mesh.node_fields.dict['transmural-vector']
longitudinal_vector = pymp.shared.array(mesh.node_fields.dict['longitudinal-vector'].shape)
longitudinal_vector[:,:] = mesh.node_fields.dict['longitudinal-vector']
t_gs = pymp.shared.array(mesh.node_fields.dict['transmural-vector'].shape)
t_gs[:,:] = 0
c_gs = pymp.shared.array(mesh.node_fields.dict['transmural-vector'].shape)
c_gs[:,:] = 0
l_gs = pymp.shared.array(mesh.node_fields.dict['longitudinal-vector'].shape)
l_gs[:,:] = 0
fibre = pymp.shared.array(mesh.node_fields.dict['fibre'].shape)
sheet = pymp.shared.array(mesh.node_fields.dict['sheet'].shape)
normal = pymp.shared.array(mesh.node_fields.dict['normal'].shape)
fibre[:,:] = mesh.node_fields.dict['fibre']
sheet[:,:] = mesh.node_fields.dict['sheet']
normal[:,:] = mesh.node_fields.dict['normal']
threadsNum = multiprocessing.cpu_count()
print('Generating new fibre projections...')
zero_tcl_counter = 0
tcl_identical_counter = 0
fig = plt.figure()
ax = fig.add_subplot(131, projection='3d')
p = ax.scatter(mesh.geometry.nodes_xyz[::10,0], mesh.geometry.nodes_xyz[::10,1], mesh.geometry.nodes_xyz[::10,2],
           c=mesh.node_fields.dict['cell-type'][::10], marker='o', s=1)
fig.colorbar(p, ax=ax, label='cell-type')
problem_tcl = pymp.shared.array(mesh.geometry.nodes_xyz.shape[0])
problem_tcl[:] = 0
# if True:
#     for node_i in range(mesh.geometry.number_of_nodes):
with pymp.Parallel(min(threadsNum, mesh.geometry.number_of_nodes)) as p1:
    for node_i in p1.range(mesh.geometry.number_of_nodes):
        # Gram-Schmidt orthogonalisation: needed otherwise the reconstruction will never recover the original vector
        t = transmural_vector[node_i,:]
        l = longitudinal_vector[node_i,:]
        t_gs[node_i,:] = t
        l_gs[node_i, :] = l - np.dot(l, t_gs[node_i, :]) * l
        c_gs[node_i, :] = np.cross(l_gs[node_i,:], t_gs[node_i,:])
        if (np.linalg.norm(t_gs[node_i, :]) < 1e-3) | (np.linalg.norm(c_gs[node_i, :]) < 1e-3) | (np.linalg.norm(l_gs[node_i, :]) < 1e-3):
            f_invariant_projection[node_i, :] = [0., 0., 0.]
            s_invariant_projection[node_i, :] = [0., 0., 0.]
            n_invariant_projection[node_i, :] = [0., 0., 0.]
            problem_tcl[node_i] = 1
            continue
        else:
            t_gs[node_i,:] = t_gs[node_i,:] / np.linalg.norm(t_gs[node_i,:])
            c_gs[node_i, :] = c_gs[node_i, :] / np.linalg.norm(c_gs[node_i, :])
            l_gs[node_i, :] = l_gs[node_i, :] / np.linalg.norm(l_gs[node_i, :])
            M = np.vstack([t_gs[node_i,:], c_gs[node_i,:], l_gs[node_i,:]]).transpose()
            if np.linalg.det(M) <=0:
                f_invariant_projection[node_i, :] = [0., 0., 0.]
                s_invariant_projection[node_i, :] = [0., 0., 0.]
                n_invariant_projection[node_i, :] = [0., 0., 0.]
                problem_tcl[node_i] = 2
            else:
                f_invariant_projection[node_i, :] = np.matmul(M, fibre[node_i, :])
                s_invariant_projection[node_i, :] = np.matmul(M, sheet[node_i, :])
                n_invariant_projection[node_i, :] = np.matmul(M, normal[node_i, :])
        # else:
        #     t_gs[node_i,:] = t_gs[node_i,:]/np.linalg.norm(t_gs[node_i,:])
        #     l_gs[node_i, :] = l_gs[node_i, :] / np.linalg.norm(l_gs[node_i, :])
        #     c_gs[node_i, :] = c_gs[node_i, :] / np.linalg.norm(c_gs[node_i, :])
        # M = np.vstack([t_gs[node_i,:], c_gs[node_i,:], l_gs[node_i,:]]).transpose()

        # else:
        #     f_invariant_projection[node_i, :] = np.matmul(np.linalg.inv(M), fibre[node_i,:])
        #     s_invariant_projection[node_i, :] = np.matmul(np.linalg.inv(M), sheet[node_i, :])
        #     n_invariant_projection[node_i, :] = np.matmul(np.linalg.inv(M), normal[node_i, :])
ax2 = fig.add_subplot(132, projection='3d')
p2 = ax2.scatter(mesh.geometry.nodes_xyz[::10,0], mesh.geometry.nodes_xyz[::10,1], mesh.geometry.nodes_xyz[::10,2],
           c=problem_tcl[::10], marker='o', s=1)
fig.colorbar(p2, ax=ax2, label='Problematic basis vectors: 1 - zero t, c, l; 2 - negative det(M), identical t, c, l')

np.savetxt(output_dir+'f_invariant_projection.csv', f_invariant_projection, delimiter=',')
np.savetxt(output_dir+'s_invariant_projection.csv', s_invariant_projection, delimiter=',')
np.savetxt(output_dir+'n_invariant_projection.csv', n_invariant_projection, delimiter=',')
np.savetxt(output_dir+'node_xyz.csv', mesh.geometry.nodes_xyz, delimiter=',')
np.savetxt(output_dir+'edges.csv', mesh.geometry.edges, delimiter=',')
np.savetxt(output_dir+'ab.csv', mesh.node_fields.dict['ab'], delimiter=',')
np.savetxt(output_dir+'tm.csv', mesh.node_fields.dict['tm'], delimiter=',')
np.savetxt(output_dir+'rt.csv', mesh.node_fields.dict['rt'], delimiter=',')
np.savetxt(output_dir+'tv.csv', mesh.node_fields.dict['tv'], delimiter=',')
lvendo_mask = (mesh.boundary_node_fields.dict['mechanical-node-label-global'] == 3).astype(int)
rvendo_mask = (mesh.boundary_node_fields.dict['mechanical-node-label-global'] == 2).astype(int)
epi_mask = (mesh.boundary_node_fields.dict['mechanical-node-label-global'] == 1).astype(int)
np.savetxt(output_dir+'node_mask_lvendo.csv', lvendo_mask, fmt='%i', delimiter=',')
np.savetxt(output_dir+'node_mask_rvendo.csv', rvendo_mask, fmt='%i', delimiter=',')
np.savetxt(output_dir+'node_mask_epi.csv', epi_mask, fmt='%i', delimiter=',')

########################################################################################################################
# Check that fibre information can be recovered.
f_global = pymp.shared.array(mesh.node_fields.dict['fibre'].shape)
s_global = pymp.shared.array(mesh.node_fields.dict['fibre'].shape)
n_global = pymp.shared.array(mesh.node_fields.dict['fibre'].shape)
threadsNum = multiprocessing.cpu_count()
print('Reconstructing fibres for comparison...')
with pymp.Parallel(min(threadsNum, mesh.geometry.number_of_nodes)) as p1:
    for node_i in p1.range(mesh.geometry.number_of_nodes):
        if problem_tcl[node_i] == 0:
            M = np.vstack([t_gs[node_i,:], c_gs[node_i,:], l_gs[node_i,:]]).transpose()
            f_global[node_i, :] = np.matmul(np.linalg.inv(M), f_invariant_projection[node_i, :])
            s_global[node_i, :] = np.matmul(np.linalg.inv(M), s_invariant_projection[node_i, :])
            n_global[node_i, :] = np.matmul(np.linalg.inv(M), n_invariant_projection[node_i, :])
# for node_i in range(mesh.geometry.number_of_nodes):
#             f_global[node_i, :] = t_gs[node_i,:] * f_invariant_projection[node_i, 0] + \
#                                   + c_gs[node_i, :] * f_invariant_projection[node_i, 1] \
#                                   + l_gs[node_i,:] * f_invariant_projection[node_i, 2]
#             s_global[node_i, :] = t_gs[node_i, :] * s_invariant_projection[node_i, 0] + \
#                           + c_gs[node_i, :] * s_invariant_projection[node_i, 1] \
#                           + l_gs[node_i, :] * s_invariant_projection[node_i, 2]
#             n_global[node_i, :] = t_gs[node_i, :] * n_invariant_projection[node_i, 0] + \
#                           + c_gs[node_i, :] * n_invariant_projection[node_i, 1] \
#                           + l_gs[node_i, :] * n_invariant_projection[node_i, 2]
        else:
            f_global[node_i,:] = fibre[node_i, :]
            s_global[node_i, :] = sheet[node_i, :]
            n_global[node_i, :] = normal[node_i, :]
print('fibre: ', fibre[-1,:])
print('t; ', transmural_vector[-1,:])
print('l: ', longitudinal_vector[-1,:])
print('fhat: ', f_invariant_projection[-1,:])
print('f_recon: ', f_global[-1,:])
print('Fibre reconstruction RMSE: ', np.linalg.norm(np.linalg.norm((f_global - mesh.node_fields.dict['fibre']), axis=1), axis=0))
# print('Max magnitude error in fibre reconstruction: ', np.amax(np.linalg.norm(f_global - mesh.node_fields.dict['fibre']), axis=1))
print(f_global - mesh.node_fields.dict['fibre'])
f_error = f_global - mesh.node_fields.dict['fibre']
print(f_error.shape)
print(np.linalg.norm(f_global - mesh.node_fields.dict['fibre'], axis=1))
f_error_magnitude = np.linalg.norm(f_error, axis=1)
print('Max magnitude error in fibre reconstruction: ', np.amax(f_error_magnitude))
ax3 = fig.add_subplot(133, projection='3d')
p3 = ax3.scatter(mesh.geometry.nodes_xyz[::10,0], mesh.geometry.nodes_xyz[::10,1], mesh.geometry.nodes_xyz[::10,2],
           c=f_error_magnitude[::10], marker='o', s=1)
fig.colorbar(p3, ax=ax3, label='Fibre vector magnitude error')
plt.show()