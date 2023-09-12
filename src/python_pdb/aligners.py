'''Functions for aligning protein strcutures.'''
from typing import Type

import numpy as np

from python_pdb.entities import Entity


def align(mobile: Type[Entity], target: Type[Entity], entity_to_move: Type[Entity] = None) -> Type[Entity]:
    '''Align entities structure in 3D space.

    Args:
        mobile: movable region that will align to target.
        target: target entity to align mobile region on to.
        entity_to_move: the entity to move based on the alignment of mobile and target region (optional).

                        It may be desirable to move an entire protein complex into a new reference frame based on the
                        alignment of domains. If this is not provided, the new entity created will just be the mobile
                        region moved onto the target region.

    Returns:
        A new view of the entity (Structure, Model, etc) aligned with the target region.

        The return type will be an entity with the same level as the `entity_to_move`, or if that is not provided than
        the same as the mobile region.

    '''
    mobile_coords = mobile.get_coordinates()
    target_coords = target.get_coordinates()
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

    region = entity_to_move.copy() if entity_to_move is not None else mobile.copy()

    _apply_transform(region, rotation, translation)

    return region


def _apply_transform(entity, rotation, translation):
    if hasattr(entity, 'position'):
        new_coords = np.dot(np.array(entity.position), rotation) + translation

        entity.pos_x = new_coords[0]
        entity.pos_y = new_coords[1]
        entity.pos_z = new_coords[2]

    else:
        for child in entity:
            _apply_transform(child, rotation, translation)
