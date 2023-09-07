'''Dataclasses to store PDB record information.'''
from dataclasses import dataclass


@dataclass(frozen=True)
class Record:
    '''Base record from which all others are inherited.'''


@dataclass(frozen=True)
class ModelRecord(Record):
    '''MODEL record.

    Attributes:
        serial_number: can optionally have a serial number associated with the model.

    '''
    serial_number: int | None

    @property
    def record_type(self) -> str:
        '''Type of record.'''
        return 'MODEL'


@dataclass(frozen=True)
class EndModelRecord(Record):
    @property
    def record_type(self) -> str:
        '''Type of record.'''
        return 'ENDMDL'


@dataclass(frozen=True)
class ParticleRecord(Record):
    atom_num: int
    atom_name: str
    alt_loc: str | None
    res_name: str
    chain_id: str
    seq_id: int
    insert_code: str | None
    x_pos: float
    y_pos: float
    z_pos: float
    occupancy: float
    b_factor: float
    element: str
    charge: str | None


@dataclass(frozen=True)
class AtomRecord(ParticleRecord):
    @property
    def record_type(self) -> str:
        '''Type of record.'''
        return 'ATOM'


@dataclass(frozen=True)
class HetAtomRecord(ParticleRecord):
    @property
    def record_type(self) -> str:
        '''Type of record.'''
        return 'HETATM'
