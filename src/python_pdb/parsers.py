'''Functions for parsing PDB files into python objects and vice versa.'''

from python_pdb.entities import Structure


def parse_pdb(contents: str) -> Structure:
    '''Create structure object from PDB file.'''
    return Structure.from_pdb(contents)


def stringify_structure(structure: Structure) -> str:
    '''Format structure as pdb atom records.'''
    return str(structure)
