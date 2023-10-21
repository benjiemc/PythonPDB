'''Functions for parsing PDB files into python objects and vice versa.'''
import warnings

import pandas as pd

from python_pdb.entities import Structure, StructureConstructionWarning
from python_pdb.formats.pdb import generate_records_from_pdb


def parse_pdb(contents: str, silent: bool = False) -> Structure:
    '''Create structure object from PDB file.

    This function is a wrapper around the Structure.from_pdb(...) classmethod that provides additional functionality of
    suppressing warnings.

    Args:
        contents: contents of a PDB file in string format.
        silent: suppress any warnings while building the structure (Default: False)

    Return:
        Structure object parsed from the PDB file.

    '''
    with warnings.catch_warnings():
        if silent:
            warnings.filterwarnings('ignore', category=StructureConstructionWarning)

        return Structure.from_pdb(contents)


def parse_pdb_to_pandas(contents: str) -> pd.DataFrame:
    '''Create pandas dataframe from the contents of a pdb file.

    Args:
        contents: the contents of a pdb file.

    Returns:
        dataframe with the following columns representing a structure:

            - record_type
            - atom_number
            - atom_name
            - alt_loc
            - residue_name
            - chain_id
            - residue_seq_id
            - residue_insert_code
            - pos_x
            - pos_y
            - pos_z
            - occupancy
            - b_factor
            - element
            - charge

        and optionally `model_index` if there are multiple models.

    '''
    entries = []
    model_id = None
    multiple_models = False

    for record in generate_records_from_pdb(contents):
        if record.record_type == 'MODEL':
            model_id = record.serial_number
            multiple_models = True

        elif record.record_type == 'ENDMDL':
            model_id = None
            multiple_models = True

        if record.record_type == 'ATOM' or record.record_type == 'HETATM':
            entries.append({
                'record_type': record.record_type,
                'atom_number': record.atom_num,
                'atom_name': record.atom_name,
                'alt_loc': record.alt_loc,
                'residue_name': record.res_name,
                'chain_id': record.chain_id,
                'residue_seq_id': record.seq_id,
                'residue_insert_code': record.insert_code,
                'pos_x': record.x_pos,
                'pos_y': record.y_pos,
                'pos_z': record.z_pos,
                'occupancy': record.occupancy,
                'b_factor': record.b_factor,
                'element': record.element,
                'charge': record.charge,
            })

            if model_id:
                entries[-1]['model_index'] = model_id

    column_names = ['record_type', 'atom_number', 'atom_name', 'alt_loc', 'residue_name', 'chain_id',
                    'residue_seq_id', 'residue_insert_code', 'pos_x', 'pos_y', 'pos_z', 'occupancy',
                    'b_factor', 'element', 'charge']

    if multiple_models:
        column_names.append('model_index')

    return pd.DataFrame(entries, columns=column_names)


def stringify_structure(structure: Structure) -> str:
    '''Format structure as pdb atom records.'''
    return str(structure)
