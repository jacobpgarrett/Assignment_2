import numpy as np
from scipy.linalg import eig
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(current_dir)
from mat_struct import *

def local_geometric_stiffness_matrix_3D_beam(L, A, I_rho, Fx2, Mx2, My1, Mz1, My2, Mz2):
    """
    local element geometric stiffness matrix
    source: p. 258 of McGuire's Matrix Structural Analysis 2nd Edition
    Given:
        material and geometric parameters:
            L, A, I_rho (polar moment of inertia)
        element forces and moments:
            Fx2, Mx2, My1, Mz1, My2, Mz2
    Context:
        load vector:
            [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]
        DOF vector:
            [u1, v1, w1, th_x1, th_y1, th_z1, u2, v2, w2, th_x2, th_y2, th_z2]
        Equation:
            [load vector] = [stiffness matrix] @ [DOF vector]
    Returns:
        12 x 12 geometric stiffness matrix k_g
    """
    k_g = np.zeros((12, 12))
    # upper triangle off diagonal terms
    k_g[0, 6] = -Fx2 / L
    k_g[1, 3] = My1 / L
    k_g[1, 4] = Mx2 / L
    k_g[1, 5] = Fx2 / 10.0
    k_g[1, 7] = -6.0 * Fx2 / (5.0 * L)
    k_g[1, 9] = My2 / L
    k_g[1, 10] = -Mx2 / L
    k_g[1, 11] = Fx2 / 10.0
    k_g[2, 3] = Mz1 / L
    k_g[2, 4] = -Fx2 / 10.0
    k_g[2, 5] = Mx2 / L
    k_g[2, 8] = -6.0 * Fx2 / (5.0 * L)
    k_g[2, 9] = Mz2 / L
    k_g[2, 10] = -Fx2 / 10.0
    k_g[2, 11] = -Mx2 / L
    k_g[3, 4] = -1.0 * (2.0 * Mz1 - Mz2) / 6.0
    k_g[3, 5] = (2.0 * My1 - My2) / 6.0
    k_g[3, 7] = -My1 / L
    k_g[3, 8] = -Mz1 / L
    k_g[3, 9] = -Fx2 * I_rho / (A * L)
    k_g[3, 10] = -1.0 * (Mz1 + Mz2) / 6.0
    k_g[3, 11] = (My1 + My2) / 6.0
    k_g[4, 7] = -Mx2 / L
    k_g[4, 8] = Fx2 / 10.0
    k_g[4, 9] = -1.0 * (Mz1 + Mz2) / 6.0
    k_g[4, 10] = -Fx2 * L / 30.0
    k_g[4, 11] = Mx2 / 2.0
    k_g[5, 7] = -Fx2 / 10.0
    k_g[5, 8] = -Mx2 / L
    k_g[5, 9] = (My1 + My2) / 6.0
    k_g[5, 10] = -Mx2 / 2.0
    k_g[5, 11] = -Fx2 * L / 30.0
    k_g[7, 9] = -My2 / L
    k_g[7, 10] = Mx2 / L
    k_g[7, 11] = -Fx2 / 10.0
    k_g[8, 9] = -Mz2 / L
    k_g[8, 10] = Fx2 / 10.0
    k_g[8, 11] = Mx2 / L
    k_g[9, 10] = (Mz1 - 2.0 * Mz2) / 6.0
    k_g[9, 11] = -1.0 * (My1 - 2.0 * My2) / 6.0
    # add in the symmetric lower triangle
    k_g = k_g + k_g.transpose()
    # add diagonal terms
    k_g[0, 0] = Fx2 / L
    k_g[1, 1] = 6.0 * Fx2 / (5.0 * L)
    k_g[2, 2] = 6.0 * Fx2 / (5.0 * L)
    k_g[3, 3] = Fx2 * I_rho / (A * L)
    k_g[4, 4] = 2.0 * Fx2 * L / 15.0
    k_g[5, 5] = 2.0 * Fx2 * L / 15.0
    k_g[6, 6] = Fx2 / L
    k_g[7, 7] = 6.0 * Fx2 / (5.0 * L)
    k_g[8, 8] = 6.0 * Fx2 / (5.0 * L)
    k_g[9, 9] = Fx2 * I_rho / (A * L)
    k_g[10, 10] = 2.0 * Fx2 * L / 15.0
    k_g[11, 11] = 2.0 * Fx2 * L / 15.0
    return k_g

def local_geometric_stiffness_matrix_3D_beam_without_interaction_terms(L, A, I_rho, Fx2):
    """
    local element geometric stiffness matrix
    source: p. 257 of McGuire's Matrix Structural Analysis 2nd Edition
    Given:
        material and geometric parameters:
            L, A, I_rho (polar moment of inertia)
        element forces and moments:
            Fx2
    Context:
        load vector:
            [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]
        DOF vector:
            [u1, v1, w1, th_x1, th_y1, th_z1, u2, v2, w2, th_x2, th_y2, th_z2]
        Equation:
            [load vector] = [stiffness matrix] @ [DOF vector]
    Returns:
        12 x 12 geometric stiffness matrix k_g
    """
    k_g = np.zeros((12, 12))
    # upper triangle off diagonal terms
    k_g[0, 6] = -Fx2 / L
    k_g[1, 5] = Fx2 / 10.0
    k_g[1, 7] = -6.0 * Fx2 / (5.0 * L)
    k_g[1, 11] = Fx2 / 10.0
    k_g[2, 4] = -Fx2 / 10.0
    k_g[2, 8] = -6.0 * Fx2 / (5.0 * L)
    k_g[2, 10] = -Fx2 / 10.0
    k_g[3, 9] = -Fx2 * I_rho / (A * L)
    k_g[4, 8] = Fx2 / 10.0
    k_g[4, 10] = -Fx2 * L / 30.0
    k_g[5, 7] = -Fx2 / 10
    k_g[5, 11] = -Fx2 * L / 30.0
    k_g[7, 11] = -Fx2 / 10.0
    k_g[8, 10] = Fx2 / 10.0
    # add in the symmetric lower triangle
    k_g = k_g + k_g.transpose()
    # add diagonal terms
    k_g[0, 0] = Fx2 / L
    k_g[1, 1] = 6.0 * Fx2 / (5.0 * L)
    k_g[2, 2] = 6.0 * Fx2 / (5.0 * L)
    k_g[3, 3] = Fx2 * I_rho / (A * L)
    k_g[4, 4] = 2.0 * Fx2 * L / 15.0
    k_g[5, 5] = 2.0 * Fx2 * L / 15.0
    k_g[6, 6] = Fx2 / L
    k_g[7, 7] = 6.0 * Fx2 / (5.0 * L)
    k_g[8, 8] = 6.0 * Fx2 / (5.0 * L)
    k_g[9, 9] = Fx2 * I_rho / (A * L)
    k_g[10, 10] = 2.0 * Fx2 * L / 15.0
    k_g[11, 11] = 2.0 * Fx2 * L / 15.0
    return k_g

def e_crit(k_e, k_g):
    """
    critical load factor
    Given:
        element stiffness matrix k_e
        element geometric stiffness matrix k_g
    Returns:
        critical load factor e
    """
    # calculate the eigenvalues of the matrix
    eigenvalues, eigenvectors = eig(k_e, k_g)
    
    positive_indicies = np.where(eigenvalues > 0)[0]
    if len(positive_indicies) > 0:
        min_index = positive_indicies[np.argmin(eigenvalues[positive_indicies])]

        lambda_crit = eigenvalues[min_index]

        deformation_vector = eigenvectors[:, min_index]
        return lambda_crit, deformation_vector
    else:
        raise ValueError("No positive eigenvalues found")    
    
def plot_structure_3D(nodes, elements, deformations, scale_factor=10):
    """
    Plots the original and deformed structure in 3D.

    :param nodes: N x 3 array of node coordinates
    :param elements: List of element connectivity (node indices)
    :param deformations: N x 3 array of nodal displacements
    :param scale_factor: Scaling factor for visibility
    """
    # Compute deformed node positions
    deformed_nodes = nodes + scale_factor * deformations

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot original structure
    for i, elem in enumerate(elements):
        x_orig = [nodes[elem[0], 0], nodes[elem[1], 0]]
        y_orig = [nodes[elem[0], 1], nodes[elem[1], 1]]
        z_orig = [nodes[elem[0], 2], nodes[elem[1], 2]]
        ax.plot(x_orig, y_orig, z_orig, 'b-', label="Original" if i == 0 else "")

    for i, elem in enumerate(elements):
        x_def = [deformed_nodes[elem[0], 0], deformed_nodes[elem[1], 0]]
        y_def = [deformed_nodes[elem[0], 1], deformed_nodes[elem[1], 1]]
        z_def = [deformed_nodes[elem[0], 2], deformed_nodes[elem[1], 2]]
        ax.plot(x_def, y_def, z_def, 'r--', label="Deformed" if i == 0 else "")


    # Labels and legend
    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")
    ax.set_zlabel("Z Coordinate")
    ax.set_title("Original vs. Deformed 3D Structure")
    ax.legend()
    plt.show()

def nonlinear_analysis(nodes, element_connect, f_appl, supports):
    del_vec, F_vec, nodevals, elementvals = mat_struct(nodes, element_connect, f_appl, supports)

    # Define the forces and moments at each node
    for i in range(len(nodevals)):
        nodevals[i].Fx = F_vec[6 * i]
        nodevals[i].Fy = F_vec[6 * i + 1]
        nodevals[i].Fz = F_vec[6 * i + 2]
        nodevals[i].Mx = F_vec[6 * i + 3]
        nodevals[i].My = F_vec[6 * i + 4]
        nodevals[i].Mz = F_vec[6 * i + 5]

    # Define the forces and moments at the end nodes of each element
    for i in range(len(elementvals)):
        elementvals[i].Fx1 = nodevals[element_connect[i][0]].Fx
        elementvals[i].Fy1 = nodevals[element_connect[i][0]].Fy
        elementvals[i].Fz1 = nodevals[element_connect[i][0]].Fz
        elementvals[i].Mx1 = nodevals[element_connect[i][0]].Mx
        elementvals[i].My1 = nodevals[element_connect[i][0]].My
        elementvals[i].Mz1 = nodevals[element_connect[i][0]].Mz
        elementvals[i].Fx2 = nodevals[element_connect[i][1]].Fx
        elementvals[i].Fy2 = nodevals[element_connect[i][1]].Fy
        elementvals[i].Fz2 = nodevals[element_connect[i][1]].Fz
        elementvals[i].Mx2 = nodevals[element_connect[i][1]].Mx
        elementvals[i].My2 = nodevals[element_connect[i][1]].My
        elementvals[i].Mz2 = nodevals[element_connect[i][1]].Mz
    
    # Define the local geometric stiffness matrix for each element
    for i in range(len(elementvals)):
        elementvals[i].k_g = local_geometric_stiffness_matrix_3D_beam(elementvals[i].L, elementvals[i].A, elementvals[i].Ip, elementvals[i].Fx2, elementvals[i].Mx2, elementvals[i].My1, elementvals[i].Mz1, elementvals[i].My2, elementvals[i].Mz2)
        elementvals[i].lambda_crit, elementvals[i].deformation_vector = e_crit(elementvals[i].k_e, elementvals[i].k_g)

    # Extract displacement components (u, v, w) for visualization
    num_nodes = len(nodes)
    deformations = np.array([[del_vec[6 * i], del_vec[6 * i + 1], del_vec[6 * i + 2]] for i in range(num_nodes)])

    # Call visualization function
    plot_structure_3D(np.array(nodes), element_connect, deformations)
    
    return del_vec, F_vec