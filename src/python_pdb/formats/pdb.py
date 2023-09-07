'''Information and tools for working with PDB file formats.

Range contants for parsing PDB files.

https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html

.. list-table:: RECORD formats
   :widths: 15 15 15 40
   :header-rows: 1

   * - COLUMNS
     - DATA TYPE
     - FIELD
     - DEFINITION
   * - 1 to 6
     - Record name
     - "ATOM  " or "HETATM"
     -
   * - 7 to 11
     - Integer
     - serial
     - Atom  serial number.
   * - 13 to 16
     - Atom
     - name
     - Atom name.
   * - 17
     - Character
     - altLoc
     - Alternate location indicator.
   * - 18 to 20
     - Residue name
     - resName
     - Residue name.
   * - 22
     - Character
     - chainID
     - Chain identifier.
   * - 23 to 26
     - Integer
     - resSeq
     - Residue sequence number.
   * - 27
     - AChar
     - iCode
     - Code for insertion of residues.
   * - 31 to 38
     - Real(8.3)
     - x
     - Orthogonal coordinates for X in Angstroms.
   * - 39 to 46
     - Real(8.3)
     - y
     - Orthogonal coordinates for Y in Angstroms.
   * - 47 to 54
     - Real(8.3)
     - z
     - Orthogonal coordinates for Z in Angstroms.
   * - 55 to 60
     - Real(6.2)
     - occupancy
     - Occupancy.
   * - 61 to 66
     - Real(6.2)
     - tempFactor
     - Temperature  factor.
   * - 77 to 78
     - LString(2)
     - element
     - Element symbol, right-justified.
   * - 79 to 80
     - LString(2)
     - charge
     - Charge  on the atom.

.. code-block:: bash

    Eg:
    ATOM     13  CA  GLU A   3      42.492 -46.706  53.609  1.00  0.00           C
    ATOM     10 2HE2 GLN A   2      42.647 -51.666  50.359  1.00  0.00           H
    ...

'''
from typing import Generator, Type

from python_pdb.records import AtomRecord, EndModelRecord, HetAtomRecord, ModelRecord, Record

RECORD_TYPE_RANGE = (0, 6)
'''Range indicating the record type (1 - 6)'''
ATOM_NUMBER_RANGE = (6, 11)
'''Range indicating the atom number (7 - 11)'''
ATOM_NAME_RANGE = (12, 16)
'''Range indicating the atom name (13 - 16)'''
ALT_LOC_RANGE = (16, 17)
'''Range indicating the alternate location details (17)'''
RESIDUE_NAME_RANGE = (17, 20)
'''Range indicating the residue name (18 - 20)'''
CHAIN_ID_RANGE = (21, 22)
'''Range indicating the chain ID (22)'''
SEQ_ID_RANGE = (22, 26)
'''Range indicating the sequence ID (23 - 26)'''
INSERT_CODE_RANGE = (26, 27)
'''Range indicating the insert code (27)'''
X_POS_RANGE = (30, 38)
'''Range indicating the X position (31 - 38)'''
Y_POS_RANGE = (38, 46)
'''Range indicating the Y position (39 - 46)'''
Z_POS_RANGE = (46, 54)
'''Range indicating the Y position (47 - 54)'''
OCCUPANCY_RANGE = (54, 60)
'''Range indicating the occupancy (55 - 60)'''
B_FACTOR_RANGE = (60, 66)
'''Range indicating the B factor (61 - 66)'''
ELEMENT_RANGE = (76, 78)
'''Range indicating the element (77 - 78)'''
CHARGE_RANGE = (78, 80)
'''Range indicating the charge (79 - 80)'''
MODEL_SERIAL_RANGE = (10, 14)
'''Range for MODEL records serial number.

.. list-table:: RECORD formats
   :widths: 15 15 15 40
   :header-rows: 1

   * - COLUMNS
     - DATA TYPE
     - FIELD
     - DEFINITION
   * - 1 to 6
     - Record name
     - "MODEL "
     -
   * - 11 - 14
     - Integer
     - serial
     - Model serial number.

'''


def format_atom_record(chain, residue, atom):
    '''Format information as a record entry for pdb files. Function will handle both ATOM and HETATM entries.'''
    # Format atom name - Taken from BioPython PDBIO module
    # Pad if:
    #     - smaller than 4 characters
    # AND - is not C, N, O, S, H, F, P, ..., one letter elements
    # AND - first character is NOT numeric (funky hydrogen naming rules)
    atom_name = atom.name.strip()
    if len(atom_name) < 4 and atom_name[:1].isalpha() and len(atom.element.strip()) < 2:
        atom_name = ' ' + atom_name

    record_type = 'ATOM' if not atom.het_atom else 'HETATM'

    return (f"{record_type:<6}{atom.number:>5} {atom_name:<4}{atom.alt_loc if atom.alt_loc else '':>1}"
            f"{residue.name:>3} {chain.name:>1}"
            f"{residue.seq_id:>4}{residue.insert_code if residue.insert_code else '':>1}   "
            f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}"
            f"{atom.occupancy:>6.2f}{atom.b_factor:>6.2f}          "
            f"{atom.element:>2}{atom.charge if atom.charge else '':>2}")


def format_model_record(model_number):
    '''Format model_number into a MODEL record.'''
    return f"MODEL     {model_number:>4}"


def generate_records_from_pdb(pdb_contents: str) -> Generator[Type[Record], None, None]:
    '''Generator to create records from a PDB file.

    Args:
        pdb_contents: contents of a PDB file.

    Yields:
        individual records from the PDB file.

    '''
    def steralize(input):
        input = input.strip()
        return None if input == '' else input

    for line in pdb_contents.split('\n'):
        record_type = steralize(line[RECORD_TYPE_RANGE[0]:RECORD_TYPE_RANGE[1]])

        if record_type == 'MODEL':
            serial_number = int(steralize(line[MODEL_SERIAL_RANGE[0]:MODEL_SERIAL_RANGE[1]]))
            yield ModelRecord(serial_number)

        if record_type == 'ENDMDL':
            yield EndModelRecord()

        if record_type == 'ATOM' or record_type == 'HETATM':
            atom_num = int(steralize(line[ATOM_NUMBER_RANGE[0]:ATOM_NUMBER_RANGE[1]]))
            atom_name = steralize(line[ATOM_NAME_RANGE[0]:ATOM_NAME_RANGE[1]])
            alt_loc = steralize(line[ALT_LOC_RANGE[0]:ALT_LOC_RANGE[1]])
            res_name = steralize(line[RESIDUE_NAME_RANGE[0]:RESIDUE_NAME_RANGE[1]])
            chain_id = steralize(line[CHAIN_ID_RANGE[0]:CHAIN_ID_RANGE[1]])
            seq_id = int(steralize(line[SEQ_ID_RANGE[0]:SEQ_ID_RANGE[1]]))
            insert_code = steralize(line[INSERT_CODE_RANGE[0]:INSERT_CODE_RANGE[1]])
            x_pos = float(steralize(line[X_POS_RANGE[0]:X_POS_RANGE[1]]))
            y_pos = float(steralize(line[Y_POS_RANGE[0]:Y_POS_RANGE[1]]))
            z_pos = float(steralize(line[Z_POS_RANGE[0]:Z_POS_RANGE[1]]))
            occupancy = float(steralize(line[OCCUPANCY_RANGE[0]:OCCUPANCY_RANGE[1]]))
            b_factor = float(steralize(line[B_FACTOR_RANGE[0]:B_FACTOR_RANGE[1]]))
            element = steralize(line[ELEMENT_RANGE[0]:ELEMENT_RANGE[1]])
            charge = steralize(line[CHARGE_RANGE[0]:CHARGE_RANGE[1]])

            if record_type == 'ATOM':
                yield AtomRecord(atom_num,
                                 atom_name,
                                 alt_loc,
                                 res_name,
                                 chain_id,
                                 seq_id,
                                 insert_code,
                                 x_pos,
                                 y_pos,
                                 z_pos,
                                 occupancy,
                                 b_factor,
                                 element,
                                 charge)

            if record_type == 'HETATM':
                yield HetAtomRecord(atom_num,
                                    atom_name,
                                    alt_loc,
                                    res_name,
                                    chain_id,
                                    seq_id,
                                    insert_code,
                                    x_pos,
                                    y_pos,
                                    z_pos,
                                    occupancy,
                                    b_factor,
                                    element,
                                    charge)
