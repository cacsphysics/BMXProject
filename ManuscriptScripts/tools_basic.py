import numpy as np


def compute_mag(x1, x2, x3):
    """ Return the magnitude of the vector
    Inputs:
        X1-component (numpy array)
        X2-component 
        X3-component
    Outputs:
        mag (numpy array)
    """
    return np.sqrt(x1**2 + x2**2 + x3**2)
