'''Functions for comparing PDB structures.'''
import numpy as np


def rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    '''Calculate Root Mean Squard Deviation between two sets of coordinates.'''
    diff = coords2 - coords1
    distance = np.sqrt(np.sum(diff * diff, axis=1))

    return np.sqrt(np.sum(distance * distance) / coords1.shape[0])
