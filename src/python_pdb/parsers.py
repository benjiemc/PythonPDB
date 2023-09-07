'''Functions for parsing PDB files into python objects and vice versa.'''
import warnings

from python_pdb.entities import Structure, StructureConstructionWarning


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


def stringify_structure(structure: Structure) -> str:
    '''Format structure as pdb atom records.'''
    return str(structure)
