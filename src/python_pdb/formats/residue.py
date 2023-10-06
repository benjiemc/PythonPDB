'''Mappings for converting between residue formats.'''


_amino_acids = [
    ('GLY', 'G'),
    ('PRO', 'P'),
    ('ALA', 'A'),
    ('VAL', 'V'),
    ('LEU', 'L'),
    ('ILE', 'I'),
    ('MET', 'M'),
    ('CYS', 'C'),
    ('PHE', 'F'),
    ('TYR', 'Y'),
    ('TRP', 'W'),
    ('HIS', 'H'),
    ('LYS', 'K'),
    ('ARG', 'R'),
    ('GLN', 'Q'),
    ('ASN', 'N'),
    ('GLU', 'E'),
    ('ASP', 'D'),
    ('SER', 'S'),
    ('THR', 'T'),
]

AMINO_ACIDS = [tlc for tlc, _ in _amino_acids]
'''List of amino acid residue one letter codes.'''

THREE_TO_ONE_CODE = {tlc: olc for tlc, olc in _amino_acids}
'''Mapping from three letter code (TLC) to one letter code (OLC) for each amino acid residue.'''


ONE_TO_THREE_CODE = {olc: tlc for tlc, olc in _amino_acids}
'''Mapping from one letter code (OLC) to three letter code (TLC) for each amino acid residue.'''


class UnknownResidueError(Exception):
    '''Raised when an unknown residue is input.'''

    def __init__(self, res: str):
        self.res = res
