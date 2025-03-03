import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.nonlinear_analysis import *
import numpy as np
import pytest

def test_e_crit_error():
    # Tests to see if input matricies have negative eigenvalues
    k_e = [[0,0],[0,0]]
    k_g = [[0,0],[0,0]]
    with pytest.raises(ValueError):
        e_crit(k_e, k_g)

# def test_e_crit():
#     # Tests to see if e_crit function works
#     k_e = [[1,0],[0,1]]
#     k_g = [[1,0],[0,1]]
#     assert e_crit(k_e, k_g) == 1

def test_example():
    # Tests to see if example works
    node0 = [0, 0, 10]  # Fixed
    node1 = [15, 0, 10]  # Midpoint (force applied)
    node2 = [15, 0, 0]  # Pin
    nodes = np.array([node0, node1, node2])

    el_1 = [0, 1, 1000, 0.3, 0.5, 0.041667, 0.010416, 0.16667, 0.02861, None]
    el_2 = [1, 2, 1000, 0.3, 0.5, 0.041667, 0.010416, 0.16667, 0.02861, [1, 0, 0]]
    element_connect = np.array([el_1, el_2], dtype=object)

    f_appl = np.array([[0,0,0,0,0,0],   # Node 0: No force
                   [0.1,0.05,-0.07,0.05,-0.1,0.25],   # Node 1: Applied force in x and y direction
                   [0,0,0,0,0,0]])  # Node 2: No force
    
    support_0 = [0, 1, 1, 1, 1, 1, 1]  # Fixed: Restraints all DOF
    support_1 = [1, 0, 0, 0, 0, 0, 0]  # Free: No restraints
    support_2 = [2, 1, 1, 1, 0, 0, 0]  # Pinned: Restraints all translation (x, y, z)
    supports = np.array([support_0, support_1, support_2])

    del_vec, F_vec, nodevals, elementvals = mat_struct(nodes, element_connect, f_appl, supports)
    
    known = np.array([0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
                      0.00000000e+00,  0.00000000e+00,  2.84049957e-03,  1.59842256e+00,
                      -1.30609177e-03, -1.47203405e-01, -1.67304182e-02,  1.82342076e-01, 
                      0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -1.66161682e-01,
                      8.79128406e-03,  1.82342076e-01])

    assert np.allclose(del_vec, known, atol=1e-6)