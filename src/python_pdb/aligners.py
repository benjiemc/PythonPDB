'''Functions for aligning protein strcutures.'''
from typing import Type

import numpy as np
import pandas as pd

from python_pdb.entities import Entity


def align_sequences(seq1: str, seq2: str,
                    match_score: float = 1.0, mismatch_score: float = -1.0, indel_score: float = -1.0) -> tuple:
    '''Align two sequences using the Needleman-Wunch algorithm.

    See here https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm for more details.

    Args:
        seq1: first string to compare
        seq2: second string to compare
        match_score: score awarded for a match (Default: 1.0)
        mismatch_score: score awarded for a mismtach (Default: -1.0)
        indel_score: score awarded for an insertion or a deletion (Default: -1.0)

    Returns:
        tuple containing the alignment (represented as a list with tuple for each pair) followed by the score given to
        the alignment.

    '''
    num_cols = len(seq1) + 1
    num_rows = len(seq2) + 1

    # Initialize the matrix
    alignment_matrix = np.empty((num_rows, num_cols))
    alignment_matrix[:] = np.nan

    alignment_matrix[0, 0] = 0

    alignment_matrix[1:, 0] = np.cumsum(np.repeat(indel_score, num_rows - 1))
    alignment_matrix[0, 1:] = np.cumsum(np.repeat(indel_score, num_cols - 1))

    direction_matrix = np.empty((num_rows, num_cols), dtype=str)

    direction_matrix[1:, 0] = 't'
    direction_matrix[0, 1:] = 'l'

    # Fill the matrix
    for i in range(1, num_rows):
        for j in range(1, num_cols):
            top_score = alignment_matrix[i - 1, j] + indel_score
            left_score = alignment_matrix[i, j - 1] + indel_score
            diagonal_score = alignment_matrix[i - 1, j - 1] + (match_score if seq1[j - 1] == seq2[i - 1] else
                                                               mismatch_score)
            scores = np.array((top_score, left_score, diagonal_score))
            score_taken = np.argmax(scores)

            alignment_matrix[i, j] = np.max(scores)
            direction_matrix[i, j] = ['t', 'l', 'd'][score_taken]

    # Trace optimal path
    reverse_alignment = []

    row_pos = num_rows - 1
    col_pos = num_cols - 1

    while (row_pos, col_pos) != (0, 0):
        if direction_matrix[row_pos, col_pos] == 't':
            reverse_alignment.append(('-', seq2[row_pos - 1]))

            row_pos -= 1

        elif direction_matrix[row_pos, col_pos] == 'l':
            reverse_alignment.append((seq1[col_pos - 1], '-'))

            col_pos -= 1

        else:
            reverse_alignment.append((seq1[col_pos - 1], seq2[row_pos - 1]))

            row_pos -= 1
            col_pos -= 1

    return (list(reversed(reverse_alignment)), alignment_matrix[-1, -1])


def align_entities(mobile_coords: np.ndarray, target_coords: np.ndarray, entity_to_move: Type[Entity]) -> Type[Entity]:
    '''Align entities structure in 3D space.

    Args:
        mobile: coordinates of atoms in movable region that will align to target.
        target: coordinates of atoms to align mobile region on to.
        entity_to_move: the entity to move based on the alignment of mobile and target region.

    Returns:
        A new view of the entity (Structure, Model, etc) aligned with the target region.

        The return type will be an entity with the same level as the `entity_to_move`, or if that is not provided than
        the same as the mobile region.

    '''
    def _apply_transform(entity, rotation, translation):
        if hasattr(entity, 'position'):
            new_coords = np.dot(np.array(entity.position), rotation) + translation

            entity.pos_x = new_coords[0]
            entity.pos_y = new_coords[1]
            entity.pos_z = new_coords[2]

        else:
            for child in entity:
                _apply_transform(child, rotation, translation)

    rotation, translation = _compute_rot_trans(mobile_coords, target_coords)

    region = entity_to_move.copy()

    _apply_transform(region, rotation, translation)

    return region


def align_pandas_structure(mobile_coords: np.ndarray,
                           target_coords: np.ndarray,
                           df_to_move: pd.DataFrame) -> pd.DataFrame:
    '''Align entities structure in 3D space.

    Args:
        mobile: coordinates of atoms in movable region that will align to target.
        target: coordinates of atoms to align mobile region on to.
        df_to_move: the dataframe to update coordinates based on the alignment of mobile and target region.

    Returns:
        A new dataframe with updated coordinates.

    '''
    def _transform(rotation, translation, row):
        return pd.Series(np.dot(row.to_numpy(), rotation) + translation)

    rotation, translation = _compute_rot_trans(mobile_coords, target_coords)

    new_df = df_to_move.copy()

    new_coords = new_df[['pos_x', 'pos_y', 'pos_z']].apply(lambda row: _transform(rotation, translation, row), axis=1)
    new_df[['pos_x', 'pos_y', 'pos_z']] = new_coords

    return new_df


def _compute_rot_trans(mobile_coords, target_coords):
    '''Compute the rotation matrix and translation vector for an alignment.'''
    # center on centroid
    mobile_average = sum(mobile_coords) / len(mobile_coords)
    target_average = sum(target_coords) / len(target_coords)

    mobile_coords = mobile_coords - mobile_average
    target_coords = target_coords - target_average

    # correlation matrix
    a = np.dot(np.transpose(mobile_coords), target_coords)
    u, d, vt = np.linalg.svd(a)
    rotation = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))

    # check if we have found a reflection
    if np.linalg.det(rotation) < 0:
        vt[2] = -vt[2]
        rotation = np.transpose(np.dot(np.transpose(vt), np.transpose(u)))

    translation = target_average - np.dot(mobile_average, rotation)

    return rotation, translation
