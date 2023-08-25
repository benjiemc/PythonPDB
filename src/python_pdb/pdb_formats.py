'''Information and tools for working with PDB file formats.'''

RECORD_TYPE_RANGE = (0, 6)
ATOM_NUMBER_RANGE = (6, 11)
ATOM_NAME_RANGE = (12, 16)
ALT_LOC_RANGE = (16, 17)
RESIDUE_NAME_RANGE = (17, 20)
CHAIN_ID_RANGE = (21, 22)
SEQ_ID_RANGE = (22, 26)
INSERT_CODE_RANGE = (26, 27)
X_POS_RANGE = (30, 38)
Y_POS_RANGE = (38, 46)
Z_POS_RANGE = (46, 54)
OCCUPANCY_RANGE = (54, 60)
B_FACTOR_RANGE = (60, 66)
ELEMENT_RANGE = (76, 78)
CHARGE_RANGE = (78, 80)
'''Range contants for parsing PDB files.

https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

Eg:
ATOM     13  CA  GLU A   3      42.492 -46.706  53.609  1.00  0.00           C
ATOM     10 2HE2 GLN A   2      42.647 -51.666  50.359  1.00  0.00           H
...

'''

MODEL_SERIAL_RANGE = (10, 14)
'''Range for MODEL records serial number.

Record Format

COLUMNS        DATA  TYPE    FIELD          DEFINITION
---------------------------------------------------------------------------------------
 1 -  6        Record name   "MODEL "
11 - 14        Integer       serial         Model serial number.
'''


def format_atom_record(chain, residue, atom):
    '''Format information as a record entry for pdb files.'''
    # Format atom name - Taken from BioPython PDBIO module
    # Pad if:
    #     - smaller than 4 characters
    # AND - is not C, N, O, S, H, F, P, ..., one letter elements
    # AND - first character is NOT numeric (funky hydrogen naming rules)
    atom_name = atom.name.strip()
    if len(atom_name) < 4 and atom_name[:1].isalpha() and len(atom.element.strip()) < 2:
        atom_name = ' ' + atom_name

    return (f"{'ATOM':<6}{atom.number:>5} {atom_name:<4}{atom.alt_loc if atom.alt_loc else '':>1}"
            f"{residue.name:>3} {chain.name:>1}"
            f"{residue.seq_id:>4}{residue.insert_code if residue.insert_code else '':>1}   "
            f"{atom.position[0]:>8.3f}{atom.position[1]:>8.3f}{atom.position[2]:>8.3f}"
            f"{atom.occupancy:>6.2f}{atom.b_factor:>6.2f}          "
            f"{atom.element:>2}{atom.charge if atom.charge else '':>2}")
