'''Classes for representing objects contained within PDB files.'''
import itertools
import warnings
from typing import Iterable, Type

import numpy as np
import pandas as pd

from python_pdb.formats.pandas import generate_records_from_pandas
from python_pdb.formats.pdb import format_atom_record, format_model_record, generate_records_from_pdb
from python_pdb.formats.residue import AMINO_ACIDS, THREE_TO_ONE_CODE
from python_pdb.records import Record


class StructureConstructionWarning(Warning):
    '''Raised if something is flagged for attention while building a structure.'''


class Entity:
    '''Base class that Atom, Residue, Chain, Model, and Structure inherit from.

    Attributes:
        children (list[Entity]): Entities that form this entity.
        parent (Entity): Entity that this entity belongs too.

    '''

    def __init__(self):
        self.children = []
        self.parent = None

    def add_child(self, child: Type['Entity']):
        '''Add child to entitity.'''
        self.children.append(child)
        child.parent = self

    def remove_child(self, child: Type['Entity']):
        '''Remove child from entity.'''
        removed_child = False

        new_children = []

        for parent_child in self:
            if child == parent_child:
                removed_child = True
                continue

            new_children.append(parent_child)

        self.children = new_children

        if not removed_child:
            warnings.warn(f'{child} not found in {self}.')

    def get_children(self) -> list[Type['Entity']]:
        '''Get children of this entity.'''
        return self.children

    def children_equal(self, other: list[Type['Entity']]) -> bool:
        return self.children == other.children

    def get_coordinates(self) -> np.ndarray:
        '''Get n x 3 matrix of atomic coordinates based on all the atoms in the Entitity. Ordering is preserved.'''
        if hasattr(self, 'position'):
            return np.array([self.position])

        coords = []
        for child_entity in self:
            coords.append(child_entity.get_coordinates())

        return np.concatenate(coords)

    def __len__(self):
        '''Number of children.'''
        return len(self.get_children())

    def __contains__(self, child: Type['Entity']):
        for parent_child in self:
            if child == parent_child:
                return True

        return False

    def __repr__(self):
        return str(self)

    def __iter__(self):
        yield from self.get_children()


class Atom(Entity):
    '''Atom representation.

    Attributes:
        name (str): atom name
        number (int): atom number in PDB file
        alt_loc (str | None): alternate location of the atom
        pos_x (float): atom's x-coordinate
        pos_y (float): atom's y-coordinate
        pos_z (float): atom's z-coordinate
        occupancy (float): atom's occupancy
        b_factor (float): b_factor from experimentally derived data
        element (str): str corresponding to the element id of this atom (periodic table style)
        charge (str | None): +/- charge on the atom
        het_atom (bool): True if the atom comes from a hetero atom record

    '''
    def __init__(self,
                 name: str,
                 number: int,
                 alt_loc: str | None,
                 pos_x: float,
                 pos_y: float,
                 pos_z: float,
                 occupancy: float,
                 b_factor: float,
                 element: str,
                 charge: str | None,
                 is_het_atom: bool = False):
        super().__init__()

        self.name = name
        self.number = number
        self.alt_loc = alt_loc
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.pos_z = pos_z
        self.occupancy = occupancy
        self.b_factor = b_factor
        self.element = element
        self.charge = charge
        self.parent = None
        self.het_atom = is_het_atom

    @property
    def position(self) -> tuple[float, float, float]:
        '''the x,y,z coordinates of the atom'''
        return (self.pos_x, self.pos_y, self.pos_z)

    def copy(self) -> 'Atom':
        '''Create a copy of the atom with the same propeties. Note the parent will be reset of the copy.'''
        return Atom(self.name,
                    self.number,
                    self.alt_loc,
                    self.pos_x,
                    self.pos_y,
                    self.pos_z,
                    self.occupancy,
                    self.b_factor,
                    self.element,
                    self.charge)

    def get_children(self):
        raise NotImplementedError("Atom's have no children.")

    def add_child(self):
        raise NotImplementedError("Atom's have no children.")

    def remove_child(self):
        raise NotImplementedError("Atom's have no children.")

    def children_equal(self, other):
        raise NotImplementedError("Atom's have no children.")

    def __contains__(self):
        raise NotImplementedError("Atom's have no children.")

    def __len__(self):
        raise NotImplementedError('Atom has no length property.')

    def __iter__(self):
        raise NotImplementedError('Atom has no length property.')

    def __str__(self):
        return f'Atom({self.name})'

    def __eq__(self, other):
        return (self.name, self.number, self.position) == (other.name, other.number, other.position)

    def __hash__(self):
        return hash((self.name, self.number, self.position))


class Residue(Entity):
    '''Residue representation.

    Attributes:
        name (str): Name of the residue (three letter code)
        seq_id (int): Sequence number of the residue.
        insert_code (str | None): Insert code of the residue if available.

    '''
    def __init__(self,
                 name: str,
                 seq_id: int,
                 insert_code: str | None):
        super().__init__()

        self.name = name
        self.seq_id = seq_id
        self.insert_code = insert_code

    def add_atom(self, atom: Atom):
        '''Add an atom to the residue.'''
        self.add_child(atom)

    def remove_atom(self, atom: Atom):
        '''Remove atom from residue.'''
        self.remove_child(atom)

    def get_atoms(self) -> list[Atom]:
        '''Return a list of all atoms in the residue, including atoms with alternate locations.'''
        return self.get_children()

    @property
    def tlc(self) -> str:
        '''Three letter code of the residue.'''
        return self.name

    @property
    def olc(self) -> str:
        '''One letter code of residue.'''
        return THREE_TO_ONE_CODE[self.name]

    def copy(self) -> 'Residue':
        '''Create a copy of the Residue. Note: the parent of the copied residue will reset.'''
        new_residue = Residue(self.name, self.seq_id, self.insert_code)

        for atom in self:
            new_residue.add_atom(atom.copy())

        return new_residue

    def __str__(self):
        return f"Residue({self.name}, {self.seq_id}{self.insert_code if self.insert_code else ''})"

    def __eq__(self, other):
        return ((self.name,
                self.seq_id,
                self.insert_code) == (other.name, other.seq_id, other.insert_code) and self.children_equal(other))

    def __hash__(self):
        return hash((self.name, self.seq_id, self.insert_code))

    def __getitem__(self, atom_name: str) -> Atom:
        '''Get atom from residue.

        Raises:
            KerError: if atom can not be located.

        '''
        for atom in self.get_atoms():
            if atom_name == atom.name:
                return atom

        raise KeyError(f"No atom {atom_name}")


class Chain(Entity):
    '''Chain representation.

    Attributes:
        name (str): name of the chain

    '''
    def __init__(self, name: str):
        super().__init__()

        self.name = name

    def add_residue(self, residue: Residue):
        '''Add residue to chain.'''
        self.add_child(residue)

    def remove_residue(self, residue: Residue):
        '''Remove residue from chain.'''
        self.remove_child(residue)

    def get_residues(self) -> list['Residue']:
        '''Get all residues on a chain.'''
        return self.get_children()

    def copy(self) -> 'Chain':
        '''Create a deep copy of the chain of the chain. Note the parent will be reset.'''
        new_chain = Chain(self.name)

        for residue in self:
            new_chain.add_residue(residue.copy())

        return new_chain

    def __str__(self):
        return f'Chain({self.name})'

    def __eq__(self, other):
        return (self.name == other.name) and self.children_equal(other)

    def __hash__(self):
        return hash(self.name)

    def __getitem__(self, seq_id: str | int) -> Residue:
        '''Get residue from chain by sequence identifier (id and insert code if applicable).

        Raises:
            KerError: if residue can not be located.

        '''
        if type(seq_id) is int:
            seq_id = str(seq_id)

        for res in self.get_residues():
            if seq_id == f"{res.seq_id}{res.insert_code if res.insert_code else ''}":
                return res

        raise KeyError(f"No residue {seq_id}")


class Model(Entity):
    '''Representation of a PDB model.

    Attributes:
        serial_number (int | None): optional serial number to identify the model.

    '''
    def __init__(self, serial_number: int | None = None):
        super().__init__()

        self.serial_number = serial_number

    def add_chain(self, chain: Chain):
        '''Add chain to model.'''
        self.add_child(chain)

    def remove_chain(self, chain: Chain):
        '''Remove chain from model.'''
        self.remove_child(chain)

    def get_chains(self) -> list[Chain]:
        '''Get all the chains in the model.'''
        return self.get_children()

    def copy(self) -> 'Model':
        '''Create a deep copy of the Model. Note parent will be reset in copy.'''
        new_model = Model(self.serial_number)

        for chain in self:
            new_model.add_chain(chain.copy())

        return new_model

    def __getitem__(self, chain_id):
        for chain in self.get_chains():
            if chain_id == chain.name:
                return chain

        raise KeyError(f'No chain named {chain_id}')

    def __eq__(self, other: 'Model'):
        return (self.serial_number == other.serial_number) and self.children_equal(other)


class Structure(Entity):
    '''PDB Structure representation.'''

    def add_model(self, model: Model):
        '''Add model to the structure.'''
        self.add_child(model)

    def remove_model(self, model: Model):
        '''Remove model from the structure.'''
        self.remove_child(model)

    def get_models(self) -> list[Model]:
        '''Get models from the structure.'''
        return self.get_children()

    def to_pandas(self) -> pd.DataFrame:
        '''Convert Structure into pandas dataframe.

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
        records = []

        multiple_models = len(self) > 1

        for model in self:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        records.append({
                            'record_type': 'ATOM' if not atom.het_atom else 'HETATM',
                            'atom_number': atom.number,
                            'atom_name': atom.name,
                            'alt_loc': atom.alt_loc,
                            'residue_name': residue.name,
                            'chain_id': chain.name,
                            'residue_seq_id': residue.seq_id,
                            'residue_insert_code': residue.insert_code,
                            'pos_x': atom.position[0],
                            'pos_y': atom.position[1],
                            'pos_z': atom.position[2],
                            'occupancy': atom.occupancy,
                            'b_factor': atom.b_factor,
                            'element': atom.element,
                            'charge': atom.charge,
                        })

                        if multiple_models:
                            records[-1]['model_index'] = model.serial_number

            column_names = ['record_type', 'atom_number', 'atom_name', 'alt_loc', 'residue_name', 'chain_id',
                            'residue_seq_id', 'residue_insert_code', 'pos_x', 'pos_y', 'pos_z', 'occupancy',
                            'b_factor', 'element', 'charge']
        if multiple_models:
            column_names.append('model_index')

        return pd.DataFrame(records, columns=column_names)

    @classmethod
    def from_pandas(cls, df: pd.DataFrame) -> 'Structure':
        '''Create Structure object from a pandas dataframe.

        Args:
            df: pandas dataframe that has the following columns representing a PDB structure:

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

                and optionaly `model_index` if there are multiple models

        Returns:
            Structure with chains, residues, and atoms from the input dataset.

        '''
        return cls._construct_from_records(generate_records_from_pandas(df))

    @classmethod
    def from_pdb(cls, contents: str) -> 'Structure':
        '''Create structure from a PDB formatted file.'''
        return cls._construct_from_records(generate_records_from_pdb(contents))

    @classmethod
    def _construct_from_records(cls, records: Iterable[Type[Record]]):
        '''Construction logic for building structure from records.'''
        structure = cls()
        model = None
        prev_res_id = None
        prev_chain_id = None

        for record in records:
            if record.record_type == 'MODEL':
                model = Model(record.serial_number)
                structure.add_model(model)

            elif record.record_type == 'ENDMDL':
                model = None
                prev_chain_id = None
                prev_res_id = None

            if record.record_type == 'ATOM' or record.record_type == 'HETATM':
                if model is None:
                    model = Model()
                    structure.add_model(model)

                if prev_chain_id != record.chain_id:
                    if record.chain_id not in [chain.name for chain in model]:
                        chain = Chain(record.chain_id)
                        model.add_chain(chain)
                    else:
                        chain = model[record.chain_id]

                    prev_chain_id = record.chain_id

                if prev_res_id != (record.chain_id, record.seq_id, record.insert_code):
                    residue = Residue(record.res_name, record.seq_id, record.insert_code)
                    chain.add_residue(residue)
                    prev_res_id = (record.chain_id, record.seq_id, record.insert_code)

                if record.alt_loc:
                    warnings.warn((f'{record}: ATOM contains alternate location. '
                                   'The `split_states` method can be used to separate these possible conformations '
                                   'into separate models.'),
                                  StructureConstructionWarning)

                residue.add_atom(Atom(record.atom_name,
                                      record.atom_num,
                                      record.alt_loc,
                                      record.x_pos,
                                      record.y_pos,
                                      record.z_pos,
                                      record.occupancy,
                                      record.b_factor,
                                      record.element,
                                      record.charge,
                                      is_het_atom=record.record_type == 'HETATM'))

        return structure

    def copy(self) -> 'Structure':
        '''Create a deep copy of the structure.'''
        new_structure = Structure()

        for model in self:
            new_structure.add_model(model.copy())

        return new_structure

    def split_states(self, all_combinations=False):
        '''Split multiple alternate locations into separate models.

        If there are multiple locations for an atom (as designated by alt_loc's) split the view of the protein into
        multiple models.

        .. warning::
           This method modifies the structure object by overriding the existing models property.

        Args:
            all_combinations: states will be grouped by alt code if False (ie all 'A's are together) or else every
                              different combination of alt codes will be used (Default: False).

        '''
        # Identify all residues with multiple conformations
        residues_with_alt_locs = []

        for model in self:
            for chain in model:
                for residue in chain:
                    alt_locs = set()

                    for atom in residue:
                        if atom.alt_loc:
                            alt_locs.add(atom.alt_loc)

                    if len(alt_locs) >= 2:
                        residues_with_alt_locs.append((model.serial_number, chain.name, residue))

        if len(residues_with_alt_locs) == 0:
            return

        # Create new residues for each conformation of residues with multiple confromations
        conformations = []

        for model_serial_number, chain_name, residue in residues_with_alt_locs:
            residue_states = {}

            none_states = []
            for atom in residue:
                if atom.alt_loc is None:
                    none_states.append(atom.copy())
                    continue
                if atom.alt_loc not in residue_states:
                    residue_states[atom.alt_loc] = Residue(residue.name, residue.seq_id, residue.insert_code)

                residue_states[atom.alt_loc].add_atom(atom.copy())

            # add atoms in residues with alt_locs that don't have alt loc themselves
            for k, v in residue_states.items():
                for atom in none_states:
                    residue_states[k].add_atom(atom)

            conformations.append([(model_serial_number, chain_name, res) for res in residue_states.values()])

        if all_combinations:
            # Create states by finding all combinations of residues with multiple conformations
            states = list(itertools.product(*conformations))
        else:
            # Pad conformations to be the same length with first conformation
            longest = len(max(conformations, key=len))
            for conformation in conformations:
                while len(conformation) < longest:
                    conformation.append(conformation[0])
            # transposes list
            states = [[row[i] for row in conformations] for i in range(len(conformations[0]))]

        # Create new models for each state
        new_models = []

        for state in states:
            res_counter = 0

            for model in self:
                new_model = Model()

                for chain in model:
                    new_chain = Chain(chain.name)

                    for residue in chain:
                        if ((model.serial_number,
                             chain.name,
                             residue.name,
                             residue.seq_id) == (state[res_counter][0],
                                                 state[res_counter][1],
                                                 state[res_counter][2].name,
                                                 state[res_counter][2].seq_id)):
                            new_chain.add_residue(state[res_counter][2])

                            if res_counter < len(state) - 1:
                                res_counter += 1

                        else:
                            new_chain.add_residue(residue.copy())

                    new_model.add_chain(new_chain)

            new_models.append(new_model)

        # Populate THIS structure with models representing each state
        self.children = new_models

    def dehydrate(self):
        '''Remove all water molecules from the Structure.'''
        for model in self:
            for chain in model:
                for res in chain:
                    for atom in res:
                        if atom.het_atom and res.name == 'HOH':
                            chain.remove_residue(res)

    def remove_non_amino_acids(self):
        '''Remove all non amino acids from the Structure.'''
        for model in self:
            for chain in model:
                for res in chain:
                    for atom in res:
                        if atom.het_atom and res.name not in AMINO_ACIDS:
                            chain.remove_residue(res)

    def __getitem__(self, index):
        return self.children[index]

    def __str__(self):
        records = []
        for model_number, model in enumerate(self, 1):
            if len(self) > 1:
                output_number = model.serial_number if model.serial_number else model_number
                records.append(format_model_record(output_number))

            for chain in model:
                for residue in chain:
                    for atom in residue:
                        records.append(format_atom_record(chain, residue, atom))

            if len(self) > 1:
                records.append('ENDMDL')

        return '\n'.join(records)

    def __eq__(self, other: 'Structure'):
        return self.children_equal(other)
