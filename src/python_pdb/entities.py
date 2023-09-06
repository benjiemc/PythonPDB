'''Classes for representing objects contained within PDB files.'''
import itertools
import warnings

import pandas as pd

from python_pdb.formats.pdb import (ALT_LOC_RANGE, ATOM_NAME_RANGE, ATOM_NUMBER_RANGE, B_FACTOR_RANGE, CHAIN_ID_RANGE,
                                    CHARGE_RANGE, ELEMENT_RANGE, INSERT_CODE_RANGE, MODEL_SERIAL_RANGE, OCCUPANCY_RANGE,
                                    RECORD_TYPE_RANGE, RESIDUE_NAME_RANGE, SEQ_ID_RANGE, X_POS_RANGE, Y_POS_RANGE,
                                    Z_POS_RANGE, format_atom_record, format_model_record)
from python_pdb.formats.residue import THREE_TO_ONE_CODE


class StructureConstructionWarning(Warning):
    '''Raised if something is flagged for attention while building a structure.'''


class Atom:
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
        parent (Residue): Residue to which this atom belongs.
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

    def __str__(self):
        return f'Atom({self.name})'

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (self.name, self.number, self.position) == (other.name, other.number, other.position)

    def __hash__(self):
        return hash((self.name, self.number, self.position))


class Residue:
    '''Residue representation.

    Attributes:
        name (str): Name of the residue (three letter code)
        seq_id (int): Sequence number of the residue.
        insert_code (str | None): Insert code of the residue if available.
        atoms (list[Atom]): Atoms that form the residue.
        parent (Chain): Chain that this residue belongs to.

    '''
    def __init__(self,
                 name: str,
                 seq_id: int,
                 insert_code: str | None):
        self.name = name
        self.seq_id = seq_id
        self.insert_code = insert_code

        self.atoms = []
        self.parent = None

    def add_atom(self, atom: Atom):
        '''Add an atom to the residue.'''
        self.atoms.append(atom)
        atom.parent = self

    def get_atoms(self) -> list[Atom]:
        '''Return a list of all atoms in the residue, including atoms with alternate locations.'''
        return self.atoms

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

    def __repr__(self):
        return str(self)

    def __iter__(self):
        yield from self.get_atoms()

    def __eq__(self, other):
        return (self.name,
                self.seq_id,
                self.insert_code,
                self.atoms) == (other.name, other.seq_id, other.insert_code, other.atoms)

    def __hash__(self):
        return hash((self.name, self.seq_id, self.insert_code))

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, atom_name: str) -> Atom:
        '''Get atom from residue.

        Raises:
            KerError: if atom can not be located.

        '''
        for atom in self.get_atoms():
            if atom_name == atom.name:
                return atom

        raise KeyError(f"No atom {atom_name}")


class Chain:
    '''Chain representation.

    Attributes:
        name (str): name of the chain
        residues (list[Residue]): residues contained in the chain
        parent (Model | None): Model that the chain belongs to.

    '''
    def __init__(self, name: str):
        self.name = name

        self.parent = None
        self.residues = []

    def add_residue(self, residue):
        '''Add residue to chain.'''
        self.residues.append(residue)
        residue.parent = self

    def get_residues(self):
        return self.residues

    def copy(self) -> 'Chain':
        '''Create a deep copy of the chain of the chain. Note the parent will be reset.'''
        new_chain = Chain(self.name)

        for residue in self:
            new_chain.add_residue(residue.copy())

        return new_chain

    def __iter__(self):
        yield from self.get_residues()

    def __str__(self):
        return f'Chain({self.name})'

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (self.name == other.name) and (self.residues == other.residues)

    def __hash__(self):
        return hash(self.name)

    def __len__(self):
        return len(self.residues)

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


class Model:
    '''Representation of a PDB model.

    Attributes:
        chains (list[Chain]): Chains contained by the model.
        parent (Structure | None): PDB structure that the Model comes from.
        serial_number (int | None): optional serial number to identify the model.

    '''
    def __init__(self, serial_number: int | None = None):
        self.parent = None
        self.chains = []
        self.serial_number = serial_number

    def add_chain(self, chain: Chain):
        '''Add chain to model.'''
        self.chains.append(chain)
        chain.parent = self

    def get_chains(self):
        return self.chains

    def copy(self):
        '''Create a deep copy of the Model. Note parent will be reset in copy.'''
        new_model = Model(self.serial_number)

        for chain in self:
            new_model.add_chain(chain.copy())

        return new_model

    def __getitem__(self, chain_id):
        for chain in self.chains:
            if chain_id == chain.name:
                return chain

        raise KeyError(f'No chain named {chain_id}')

    def __iter__(self):
        yield from self.get_chains()

    def __len__(self):
        return len(self.chains)

    def __eq__(self, other: 'Model'):
        return (self.serial_number == other.serial_number) and (self.chains == other.chains)


class Structure:
    '''PDB Structure representation.

    Attributes:
        models (list[Model]): models contained within the PDB structure.

    '''
    def __init__(self):
        self.models = []

    def add_model(self, model: Model):
        '''Add model to structure.'''
        self.models.append(model)

    def to_pandas(self) -> pd.DataFrame:
        '''Convert Structure into pandas dataframe.

        Returns:
            pd.Dataframe: dataframe with the following columns representing a structure:

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
    def from_pandas(cls, df: pd.DataFrame):
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
        new_structure = cls()

        if 'model_index' in df.columns:
            model = None

        else:
            model = Model()
            new_structure.add_model(model)

        prev_res_id = None
        prev_chain_id = None
        prev_model_index = None

        for record in df.itertuples():

            if hasattr(record, 'model_index'):
                if record.model_index != prev_model_index:
                    model = Model(record.model_index)
                    new_structure.add_model(model)
                    prev_model_index = record.model_index

                    prev_chain_id = None
                    prev_res_id = None

            if prev_chain_id != record.chain_id:
                chain = Chain(record.chain_id)
                model.add_chain(chain)
                prev_chain_id = record.chain_id

            if prev_res_id != (record.chain_id, record.residue_seq_id, record.residue_insert_code):
                residue = Residue(record.residue_name,
                                  record.residue_seq_id,
                                  record.residue_insert_code if record.residue_insert_code != '' else None)
                chain.add_residue(residue)

                prev_res_id = (record.chain_id, record.residue_seq_id, record.residue_insert_code)

            if record.alt_loc:
                warnings.warn((f'{record}: ATOM contains alternate location. '
                               'The `split_states` method can be used to separate these possible conformations '
                               'into separate models.'),
                              StructureConstructionWarning)

            residue.add_atom(Atom(record.atom_name,
                                  record.atom_number,
                                  record.alt_loc if record.alt_loc != '' else None,
                                  record.pos_x,
                                  record.pos_y,
                                  record.pos_z,
                                  record.occupancy,
                                  record.b_factor,
                                  record.element,
                                  record.charge if record.charge != '' else None))

        return new_structure

    @classmethod
    def from_pdb(cls, contents: str) -> 'Structure':
        '''Create structure from a PDB formatted file.'''
        structure = cls()
        model = None
        prev_res_id = None
        prev_chain_id = None

        for record in contents.split('\n'):
            record_type = record[RECORD_TYPE_RANGE[0]:RECORD_TYPE_RANGE[1]].strip()

            if record_type == 'MODEL':
                serial_number = int(record[MODEL_SERIAL_RANGE[0]:MODEL_SERIAL_RANGE[1]].strip())

                model = Model(serial_number)
                structure.add_model(model)

            elif record_type == 'ENDMDL':
                model = None
                prev_chain_id = None
                prev_res_id = None

            if record_type == 'ATOM' or record_type == 'HETATM':
                if model is None:
                    model = Model()
                    structure.add_model(model)

                atom_num = int(record[ATOM_NUMBER_RANGE[0]:ATOM_NUMBER_RANGE[1]].strip())
                atom_name = record[ATOM_NAME_RANGE[0]:ATOM_NAME_RANGE[1]].strip()
                alt_loc = record[ALT_LOC_RANGE[0]:ALT_LOC_RANGE[1]].strip()
                res_name = record[RESIDUE_NAME_RANGE[0]:RESIDUE_NAME_RANGE[1]].strip()
                chain_id = record[CHAIN_ID_RANGE[0]:CHAIN_ID_RANGE[1]].strip()
                seq_id = int(record[SEQ_ID_RANGE[0]:SEQ_ID_RANGE[1]].strip())
                insert_code = record[INSERT_CODE_RANGE[0]:INSERT_CODE_RANGE[1]].strip()
                x_pos = float(record[X_POS_RANGE[0]:X_POS_RANGE[1]].strip())
                y_pos = float(record[Y_POS_RANGE[0]:Y_POS_RANGE[1]].strip())
                z_pos = float(record[Z_POS_RANGE[0]:Z_POS_RANGE[1]].strip())
                occupancy = float(record[OCCUPANCY_RANGE[0]:OCCUPANCY_RANGE[1]].strip())
                b_factor = float(record[B_FACTOR_RANGE[0]:B_FACTOR_RANGE[1]].strip())
                element = record[ELEMENT_RANGE[0]:ELEMENT_RANGE[1]].strip()
                charge = record[CHARGE_RANGE[0]:CHARGE_RANGE[1]].strip()

                if prev_chain_id != chain_id:
                    chain = Chain(chain_id)
                    model.add_chain(chain)
                    prev_chain_id = chain_id

                if prev_res_id != (chain_id, seq_id, insert_code):
                    residue = Residue(res_name, seq_id, insert_code if insert_code != '' else None)
                    chain.add_residue(residue)
                    prev_res_id = (chain_id, seq_id, insert_code)

                if alt_loc:
                    warnings.warn((f'{record}: ATOM contains alternate location. '
                                   'The `split_states` method can be used to separate these possible conformations '
                                   'into separate models.'),
                                  StructureConstructionWarning)

                residue.add_atom(Atom(atom_name,
                                      atom_num,
                                      alt_loc if alt_loc != '' else None,
                                      x_pos,
                                      y_pos,
                                      z_pos,
                                      occupancy,
                                      b_factor,
                                      element,
                                      charge if charge != '' else None,
                                      is_het_atom=record_type == 'HETATM'))

        return structure

    def copy(self) -> 'Structure':
        '''Create a deep copy of the structure.'''
        new_structure = Structure()

        for model in self:
            new_structure.add_model(model.copy())

        return new_structure

    def split_states(self):
        '''Split multiple alternate locations into separate models.

        If there are multiple locations for an atom (as designated by alt_loc's) split the view of the protein into
        multiple models.

        ..warning::
           This method modifies the structure object by overriding the existing models property.

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

            for atom in residue:
                if atom.alt_loc not in residue_states:
                    residue_states[atom.alt_loc] = Residue(residue.name, residue.seq_id, residue.insert_code)

                residue_states[atom.alt_loc].add_atom(atom.copy())

            conformations.append([(model_serial_number, chain_name, res) for res in residue_states.values()])

        # Create states by finding all combinations of residues with multiple conformations
        states = list(itertools.product(*conformations))

        # Create new models for each state
        new_models = []

        for state in states:
            res_counter = 0

            for model in self:
                new_model = Model()

                for chain in model:
                    new_chain = Chain(chain.name)

                    for residue in chain:
                        if (model.serial_number, chain.name, residue) == state[res_counter]:
                            new_chain.add_residue(state[res_counter][-1])

                            if res_counter < len(state) - 1:
                                res_counter += 1

                        else:
                            new_chain.add_residue(residue.copy())

                    new_model.add_chain(new_chain)

            new_models.append(new_model)

        # Populate THIS structure with models representing each state
        self.models = new_models

    def __getitem__(self, index):
        return self.models[index]

    def __iter__(self):
        yield from self.models

    def __len__(self):
        return len(self.models)

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
        return self.models == other.models
