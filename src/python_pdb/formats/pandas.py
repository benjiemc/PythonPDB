'''Formatting for pandas formated structures.

The pandas dataframe should have the following columns to represent a PDB structure:

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

and optionaly `model_index` if there are multiple models.

'''
from typing import Generator, Type

import pandas as pd

from python_pdb.records import AtomRecord, EndModelRecord, HetAtomRecord, ModelRecord, Record


def generate_records_from_pandas(df: pd.DataFrame) -> Generator[Type[Record], None, None]:
    '''Generate pdb record from a pandas dataframe representation.'''
    handle_models = 'model_index' in df.columns
    unended_model = False
    prev_model_index = None

    for row in df.itertuples():
        if handle_models:
            if row.model_index != prev_model_index:
                if unended_model:
                    yield EndModelRecord()
                    unended_model = False

                yield ModelRecord(row.model_index)
                unended_model = True
                prev_model_index = row.model_index

        if row.record_type == 'ATOM':
            yield AtomRecord(row.atom_number,
                             row.atom_name,
                             row.alt_loc,
                             row.residue_name,
                             row.chain_id,
                             row.residue_seq_id,
                             row.residue_insert_code,
                             row.pos_x,
                             row.pos_y,
                             row.pos_z,
                             row.occupancy,
                             row.b_factor,
                             row.element,
                             row.charge)

        if row.record_type == 'HETATM':
            yield HetAtomRecord(row.atom_number,
                                row.atom_name,
                                row.alt_loc,
                                row.residue_name,
                                row.chain_id,
                                row.residue_seq_id,
                                row.residue_insert_code,
                                row.pos_x,
                                row.pos_y,
                                row.pos_z,
                                row.occupancy,
                                row.b_factor,
                                row.element,
                                row.charge)

    if handle_models and unended_model:
        yield EndModelRecord()
