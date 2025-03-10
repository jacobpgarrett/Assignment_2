import numpy as np

# Establish Node Class
class node:
    def __init__(self, x, y, z, F=None, disp=None, Fx=0, Fy=0, Fz=0, Mx=0, My=0, Mz=0):
        self.x = x
        self.y = y
        self.z = z
        self.F = np.zeros(6) if F is None else np.array(F)
        self.disp = np.zeros(6) if disp is None else np.array(disp)
        self.supported_dofs = []

        # Define forces and moments at each node
        self.Fx = Fx
        self.Fy = Fy
        self.Fz = Fz
        self.Mx = Mx
        self.My = My
        self.Mz = Mz

    def set_support(self, support):
        self.supported_dofs = [i for i, val in enumerate(support) if val]

# Establish Element Class
class element:
    def __init__(self, node1, node2, E, nu, A, Iz, Iy, Ip, J, z_axis, Fx1=0, Fy1=0, Fz1=0, Mx1=0, My1=0, Mz1=0, Fx2=0, Fy2=0, Fz2=0, Mx2=0, My2=0, Mz2=0, k_g=None, lambda_crit=None, deformation_vector=None):
        self.node1 = node1 # First node and subsequent coordinates
        self.x1 = node1.x
        self.y1 = node1.y
        self.z1 = node1.z
        self.node2 = node2 # Second node and subsequent coordinates
        self.x2 = node2.x
        self.y2 = node2.y
        self.z2 = node2.z    
        self.E = E # Modulus of elasticity
        self.nu = nu # Poisson's Ratio
        self.A = A # Cross-sectional area
        self.Iz = Iz # Moment of inertia about the z-axis
        self.Iy = Iy # Moment of inertia about the y-axis
        self.Ip = Ip # Polar moment of inertia
        self.J = J 
        self.z_axis = z_axis # Local z axis
        self.L = np.sqrt((node1.x - node2.x)**2 + (node1.y - node2.y)**2 + (node1.z - node2.z)**2) # Computes length

         # Check for duplicate nodes
        if self.L == 0:
            raise ValueError("Duplicate nodes detected: Nodes cannot be the same.")

        # Calculates element properties
        self.k_e = self.local_elastic_stiffness_matrix_3D_beam() # Obtains local elastic stiffness matrix    
        self.gamma = self.rotation_matrix_3D() # Rotation Matrix
        self.Gamma = self.transformation_matrix_3D() # Transformation Matrix
        self.k_global = self.Gamma.T @ self.k_e @ self.Gamma # Global Stiffness Matrix

        # Define forces and moments at each node in global coordinates
        self.Fx1 = Fx1
        self.Fy1 = Fy1
        self.Fz1 = Fz1
        self.Mx1 = Mx1
        self.My1 = My1
        self.Mz1 = Mz1
        self.Fx2 = Fx2
        self.Fy2 = Fy2
        self.Fz2 = Fz2
        self.Mx2 = Mx2
        self.My2 = My2
        self.Mz2 = Mz2

        self.k_g = k_g
        self.lambda_crit = lambda_crit
        self.deformation_vector = deformation_vector
    
    def local_elastic_stiffness_matrix_3D_beam(self) -> np.ndarray:
        """
        local element elastic stiffness matrix
        source: p. 73 of McGuire's Matrix Structural Analysis 2nd Edition
        Given:
            material and geometric parameters:
                A, L, Iy, Iz, J, nu, E
        Context:
            load vector:
                [Fx1, Fy1, Fz1, Mx1, My1, Mz1, Fx2, Fy2, Fz2, Mx2, My2, Mz2]
            DOF vector:
                [u1, v1, w1, th_x1, th_y1, th_z1, u2, v2, w2, th_x2, th_y2, th_z2]
            Equation:
                [load vector] = [stiffness matrix] @ [DOF vector]
        Returns:
            12 x 12 elastic stiffness matrix k_e
        """
        k_e = np.zeros((12, 12))
        # Axial terms - extension of local x axis
        axial_stiffness = self.E * self.A / self.L
        k_e[0, 0] = axial_stiffness
        k_e[0, 6] = -axial_stiffness
        k_e[6, 0] = -axial_stiffness
        k_e[6, 6] = axial_stiffness
        # Torsion terms - rotation about local x axis
        torsional_stiffness = self.E * self.J / (2.0 * (1 + self.nu) * self.L)
        k_e[3, 3] = torsional_stiffness
        k_e[3, 9] = -torsional_stiffness
        k_e[9, 3] = -torsional_stiffness
        k_e[9, 9] = torsional_stiffness
        # Bending terms - bending about local z axis
        k_e[1, 1] = self.E * 12.0 * self.Iz / self.L ** 3.0
        k_e[1, 7] = self.E * -12.0 * self.Iz / self.L ** 3.0
        k_e[7, 1] = self.E * -12.0 * self.Iz / self.L ** 3.0
        k_e[7, 7] = self.E * 12.0 * self.Iz / self.L ** 3.0
        k_e[1, 5] = self.E * 6.0 * self.Iz / self.L ** 2.0
        k_e[5, 1] = self.E * 6.0 * self.Iz / self.L ** 2.0
        k_e[1, 11] = self.E * 6.0 * self.Iz / self.L ** 2.0
        k_e[11, 1] = self.E * 6.0 * self.Iz / self.L ** 2.0
        k_e[5, 7] = self.E * -6.0 * self.Iz / self.L ** 2.0
        k_e[7, 5] = self.E * -6.0 * self.Iz / self.L ** 2.0
        k_e[7, 11] = self.E * -6.0 * self.Iz / self.L ** 2.0
        k_e[11, 7] = self.E * -6.0 * self.Iz / self.L ** 2.0
        k_e[5, 5] = self.E * 4.0 * self.Iz / self.L
        k_e[11, 11] = self.E * 4.0 * self.Iz / self.L
        k_e[5, 11] = self.E * 2.0 * self.Iz / self.L
        k_e[11, 5] = self.E * 2.0 * self.Iz / self.L
        # Bending terms - bending about local y axis
        k_e[2, 2] = self.E * 12.0 * self.Iy / self.L ** 3.0
        k_e[2, 8] = self.E * -12.0 * self.Iy / self.L ** 3.0
        k_e[8, 2] = self.E * -12.0 * self.Iy / self.L ** 3.0
        k_e[8, 8] = self.E * 12.0 * self.Iy / self.L ** 3.0
        k_e[2, 4] = self.E * -6.0 * self.Iy / self.L ** 2.0
        k_e[4, 2] = self.E * -6.0 * self.Iy / self.L ** 2.0
        k_e[2, 10] = self.E * -6.0 * self.Iy / self.L ** 2.0
        k_e[10, 2] = self.E * -6.0 * self.Iy / self.L ** 2.0
        k_e[4, 8] = self.E * 6.0 * self.Iy / self.L ** 2.0
        k_e[8, 4] = self.E * 6.0 * self.Iy / self.L ** 2.0
        k_e[8, 10] = self.E * 6.0 * self.Iy / self.L ** 2.0
        k_e[10, 8] = self.E * 6.0 * self.Iy / self.L ** 2.0
        k_e[4, 4] = self.E * 4.0 * self.Iy / self.L
        k_e[10, 10] = self.E * 4.0 * self.Iy / self.L
        k_e[4, 10] = self.E * 2.0 * self.Iy / self.L
        k_e[10, 4] = self.E * 2.0 * self.Iy / self.L
        return k_e
    
    def rotation_matrix_3D(self):
        """
        3D rotation matrix
        source: Chapter 5.1 of McGuire's Matrix Structural Analysis 2nd Edition
        Given:
            x, y, z coordinates of the ends of two beams: x1, y1, z1, x2, y2, z2
            optional: reference z vector direction v_temp to orthonormalize the local y and z axis
                if v_temp is not given, VVVV
        Compute:
            where l, m, n are defined as direction cosines:
            gamma = [[lx'=cos alpha_x', mx'=cos beta_x', nx'=cos gamma_x'],
                    [ly'=cos alpha_y', my'=cos beta_y', ny'=cos gamma_y'],
                    [lz'=cos alpha_z', mz'=cos beta_z', nz'=cos gamma_z']]
        """
        L = np.sqrt((self.x2 - self.x1) ** 2.0 + (self.y2 - self.y1) ** 2.0 + (self.z2 - self.z1) ** 2.0)
        lxp = (self.x2 - self.x1) / L
        mxp = (self.y2 - self.y1) / L
        nxp = (self.z2 - self.z1) / L
        local_x = np.asarray([lxp, mxp, nxp])

        v_temp = self.z_axis
        # choose a vector to orthonormalize the y axis if one is not given
        if v_temp is None:
            # if the beam is oriented vertically, switch to the global y axis
            if np.isclose(lxp, 0.0) and np.isclose(mxp, 0.0):
                v_temp = np.array([0, 1.0, 0.0])
            else:
                # otherwise use the global z axis
                v_temp = np.array([0, 0, 1.0])
        else:
            # check to make sure that given v_temp is a unit vector
            check_unit_vector(v_temp)
            # check to make sure that given v_temp is not parallel to the local x axis
            check_parallel(local_x, v_temp)
        
        # compute the local y axis
        local_y = np.cross(v_temp, local_x)
        local_y = local_y / np.linalg.norm(local_y)

        # compute the local z axis
        local_z = np.cross(local_x, local_y)
        local_z = local_z / np.linalg.norm(local_z)

        # assemble R
        gamma = np.vstack((local_x, local_y, local_z))
        
        return gamma

    def transformation_matrix_3D(self) -> np.ndarray:
        """
        3D transformation matrix
        source: Chapter 5.1 of McGuire's Matrix Structural Analysis 2nd Edition
        Given:
            gamma -- the 3x3 rotation matrix
        Compute:
            Gamma -- the 12x12 transformation matrix
        """
        Gamma = np.zeros((12, 12))
        Gamma[0:3, 0:3] = self.gamma
        Gamma[3:6, 3:6] = self.gamma
        Gamma[6:9, 6:9] = self.gamma
        Gamma[9:12, 9:12] = self.gamma
        return Gamma

def check_unit_vector(vec: np.ndarray):
    """
        Checks to ensure that a vector is a unit vector
    """
    if np.isclose(np.linalg.norm(vec), 1.0):
        return
    else:
        raise ValueError("Expected a unit vector for reference vector.")

def check_parallel(vec_1: np.ndarray, vec_2: np.ndarray):
    """
        Checks to ensure that two vectors are not parallel
    """
    if np.isclose(np.linalg.norm(np.cross(vec_1, vec_2)), 0.0):
        raise ValueError("Reference vector is parallel to beam axis.")
    else:
        return

def check_singular(matrix: np.ndarray):
    """
    Checks to see if a matrix is singular
    """
    if np.linalg.cond(matrix) > 1e10:
        raise ValueError("Singular stiffness matrix detected. Check supports.")

def mat_struct(nodes, element_connect, f_appl, supports):
    '''
    This function performs matrix structural analysis on an input geometry

    inputs:
    nodevals: i x 4 matrix where i is the number of nodes in the geometry
        (i, 0) = x coordinate of node
        (i, 1) = y coordinate of node
        (i, 2) = z coordinate of node
        (i, 3) = force applied to node, as a 6x1 vector [Fx, Fy, Fz, Mx, My, Mz]
    element_connect: i x 10 matrix where i is the number of elements in the geometry
        (i, 1) = first node of element
        (i, 2) = second node of element
        (i, 3) = E modulus of elasticity
        (i, 4) = nu Poisson's ratio
        (i, 5) = A cross-sectional area
        (i, 6) = Iz moment of inertia about the z-axis
        (i, 7) = Iy moment of inertia about the y-axis
        (i, 8) = J polar moment of inertia
        (i, 9) = set z-axis of the element
    f_appl: i x 6 matrix where i is the number of nodes in the geometry
        (i, 0) = force applied to node in the x direction
        (i, 1) = force applied to node in the y direction
        (i, 2) = force applied to node in the z direction
        (i, 3) = moment applied to node about the x axis
        (i, 4) = moment applied to node about the y axis
        (i, 5) = moment applied to node about the z axis
    support: i x 7 matrix where i is the number of nodes in the geometry
        (i, 0) = node number
        (i, 1) = Fx DOF (1 if restricted 0 if not)
        (i, 2) = Fy DOF (1 if restricted 0 if not)
        (i, 3) = Fz DOF (1 if restricted 0 if not)
        (i, 4) = Mx DOF (1 if restricted 0 if not)
        (i, 5) = My DOF (1 if restricted 0 if not)
        (i, 6) = Mz DOF (1 if restricted 0 if not) 
    '''

    # create nodes
    nodevals = [node(nodes[i][0], nodes[i][1], nodes[i][2], f_appl[i]) for i in range(len(nodes))]
    for support in supports:
        node_index = support[0]
        nodevals[node_index].set_support(support[1:])
    
    # create elements
    elementvals = [element(nodevals[element_connect[i][0]], nodevals[element_connect[i][1]], element_connect[i][2], 
                          element_connect[i][3], element_connect[i][4], element_connect[i][5], element_connect[i][6], 
                          element_connect[i][7], element_connect[i][8], element_connect[i][9]) for i in range(len(element_connect))]
    
    # create local stiffness matricies for each element
    lmat_local = [elementvals[i].k_e for i in range(len(elementvals))]

    # convert local stiffness matricies to global coordinates
    gamma = [elementvals[i].gamma for i in range(len(elementvals))]
    Gamma = [elementvals[i].Gamma for i in range(len(elementvals))]
    lmat_global = [elementvals[i].k_global for i in range(len(elementvals))]

    # Create global stiffness matrix
    n_squares = int(6*len(nodevals))
    k_global = np.zeros((n_squares, n_squares))
    for i, elem in enumerate(elementvals):
        node1_index, node2_index = element_connect[i][0], element_connect[i][1]
        dof_indices = np.array([node1_index * 6 + j for j in range(6)] + [node2_index * 6 + j for j in range(6)])
        for row_local, row_global in enumerate(dof_indices):
            for col_local, col_global in enumerate(dof_indices):
                k_global[row_global, col_global] += lmat_global[i][row_local, col_local]

     # Apply boundary conditions
    supported_dofs = []
    unsupported_dofs = [i for i in range(len(nodevals) * 6)]

    for support in supports:
        node_index = support[0]
        dof_indices = np.array([node_index * 6 + j for j in range(6)])
        for i, is_restricted in enumerate(support[1:]):
            if is_restricted:
                supported_dofs.append(dof_indices[i])

    unsupported_dofs = [dof for dof in unsupported_dofs if dof not in supported_dofs]

    # Partition stiffness matrix
    k_uu = k_global[np.ix_(unsupported_dofs, unsupported_dofs)]
    f_u = np.concatenate([np.array(node.F).flatten() for node in nodevals])[unsupported_dofs]
    check_singular(k_uu)
    
    # Solve for displacements
    del_u = np.linalg.solve(k_uu, f_u)
    del_f = np.zeros(len(nodevals) * 6)
    del_f[unsupported_dofs] = del_u

    # Calculate forces at all nodes
    f_all = k_global @ del_f

    return del_f, f_all, nodevals, elementvals