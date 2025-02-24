import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from src.mat_struct import mat_struct

# Define node coordinates
node0 = [0, 0, 0]  # Pinned
node1 = [2, 0, 0]  # Midpoint (force applied)
node2 = [6, 0, 0]  # Pin
nodes = np.array([node0, node1, node2])

# Define elements (connectivity: [node1, node2, E, v, A, Iz, Iy, Ip, J, z_axis])
el_1 = [0, 1, 200e9, 0.3, 0.01, 8.33e-6, 8.33e-6, 0, 0, [0, 0, 1]]
el_2 = [1, 2, 200e9, 0.3, 0.01, 8.33e-6, 8.33e-6, 0, 0, [0, 0, 1]]
element_connect = np.array([el_1, el_2], dtype=object)

# Apply force at node 1 (in the -y direction)
f_appl = np.array([[0,0,0,0,0,0],   # Node 0: No force
                   [1,1,0,0,0,0],  # Node 1: Force of -1 in y-direction
                   [0,0,0,0,0,0]])  # Node 2: No force

# Define support conditions (Node index, Fixed, Pinned)
support_0 = [0, 1, 0]  # Pinned: Restrains all translation (x, y, z)
support_2 = [2, 0, 1]  # Roller: Restrains only y-translation
supports = np.array([support_0, support_2])

# Run analysis
del_vec, F_vec = mat_struct(nodes, element_connect, f_appl, supports)

print("Displacements:\n", del_vec)
print("Forces:\n", F_vec)