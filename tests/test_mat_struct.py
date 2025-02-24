import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from src.mat_struct import *
import numpy as np
import pytest

def test_duplicate():
    # Define node coordinates
    node0 = [0, 0, 0]  # Pinned
    test_element = element(node0, node0, 1, 1, 1, 1, 1, 1, 1, [0, 0, 1])
    with pytest.raises(ValueError):
        test_element.check_duplicate()

def test_vector():
    # Tests to see if a vector is a unit vector
    test_vector = [0, 3, 1]
    with pytest.raises(ValueError):
        check_unit_vector(test_vector)

def test_parallel():
    # Tests to see if two vectors are parallel
    vector1 = np.array([1, 0, 0])
    vector2 = np.array([0, 1, 0])
    with pytest.raises(ValueError):
        check_parallel(vector1, vector2)