import warnings
from unittest import TestCase

import numpy as np
import pandas as pd

from python_pdb.entities import Atom, Chain, Model, Residue, Structure, StructureConstructionWarning
from python_pdb.records import AtomRecord, EndModelRecord, ModelRecord


class TestAtom(TestCase):
    def test_copy(self):
        mock_atom = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)

        copy = mock_atom.copy()

        self.assertEqual(mock_atom, copy)
        mock_atom.pos_x += 1.00
        self.assertNotEqual(mock_atom, copy)

    def test_get_coordinates(self):
        mock_atom = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)
        np.testing.assert_array_equal(mock_atom.get_coordinates(), np.array([[0.00, 0.00, 0.00]]))


class TestResidue(TestCase):
    def test_copy(self):
        mock_residue = Residue('ALA', 1, None)

        copy = mock_residue.copy()

        self.assertEqual(mock_residue, copy)
        mock_residue.name = 'GLY'
        self.assertNotEqual(mock_residue, copy)

    def test_equals_children(self):
        mock_atom_1 = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)

        mock_res_1 = Residue('ALA', 1, None)
        mock_res_2 = Residue('ALA', 1, None)

        mock_res_1.add_atom(mock_atom_1)
        mock_res_2.add_atom(mock_atom_1)

        self.assertEqual(mock_res_1, mock_res_2)

    def test_not_equals_children(self):
        mock_atom_1 = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)
        mock_atom_2 = Atom('O', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'O', None)

        mock_res_1 = Residue('ALA', 1, None)
        mock_res_2 = Residue('ALA', 1, None)

        mock_res_1.add_atom(mock_atom_1)
        mock_res_2.add_atom(mock_atom_2)

        self.assertNotEqual(mock_res_1, mock_res_2)

    def test_remove_atom(self):
        mock_atom_1 = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)

        mock_res_1 = Residue('ALA', 1, None)
        mock_res_1.add_atom(mock_atom_1)

        self.assertEqual(len(mock_res_1), 1)

        mock_res_1.remove_atom(mock_atom_1)

        self.assertEqual(len(mock_res_1), 0)

    def test_remove_atom_warning(self):
        mock_atom_1 = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)

        mock_res_1 = Residue('ALA', 1, None)

        with self.assertWarns(Warning):
            mock_res_1.remove_atom(mock_atom_1)

    def test_in(self):
        mock_atom = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)
        mock_res = Residue('ALA', 1, None)

        mock_res.add_atom(mock_atom)

        self.assertTrue(mock_atom in mock_res)

    def test_get_coordinates(self):
        mock_atom_1 = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'N', None)
        mock_atom_2 = Atom('O', 1, None, 0.00, 0.00, 0.00, 1.00, 7.9, 'O', None)

        mock_res = Residue('ALA', 1, None)
        mock_res.add_atom(mock_atom_1)
        mock_res.add_atom(mock_atom_2)

        np.testing.assert_array_equal(mock_res.get_coordinates(), np.array([[0.00, 0.00, 0.00],
                                                                            [0.00, 0.00, 0.00]]))


class TestChain(TestCase):
    def test_copy(self):
        mock_chain = Chain('A')

        copy = mock_chain.copy()

        self.assertEqual(mock_chain, copy)

        mock_chain.name = 'B'

        self.assertNotEqual(mock_chain, copy)

    def test_remove_residue(self):
        mock_chain = Chain('A')
        mock_res = Residue('ALA', 1, None)

        mock_chain.add_residue(mock_res)

        self.assertEqual(len(mock_chain), 1)

        mock_chain.remove_residue(mock_res)

        self.assertEqual(len(mock_chain), 0)

    def test_remove_residue_warning(self):
        mock_chain = Chain('A')
        mock_res = Residue('ALA', 1, None)

        with self.assertWarns(Warning):
            mock_chain.remove_residue(mock_res)

    def test_in(self):
        mock_chain = Chain('A')
        mock_res = Residue('ALA', 1, None)

        mock_chain.add_residue(mock_res)

        self.assertTrue(mock_res in mock_chain)

    def test_get_coordinates(self):
        mock_atom_1 = Atom('N', 1, None, 1.00, 0.00, 0.00, 1.00, 7.9, 'N', None)
        mock_atom_2 = Atom('O', 1, None, 1.00, 0.00, 0.00, 1.00, 7.9, 'O', None)

        mock_res = Residue('ALA', 1, None)
        mock_res.add_atom(mock_atom_1.copy())
        mock_res.add_atom(mock_atom_2.copy())

        mock_chain = Chain('A')
        mock_chain.add_residue(mock_res.copy())
        mock_chain.add_residue(mock_res.copy())

        np.testing.assert_array_equal(mock_chain.get_coordinates(), np.array([[1.00, 0.00, 0.00],
                                                                              [1.00, 0.00, 0.00],
                                                                              [1.00, 0.00, 0.00],
                                                                              [1.00, 0.00, 0.00]]))


class TestModel(TestCase):
    def test_in(self):
        mock_model = Model()
        mock_chain = Chain('A')

        mock_model.add_chain(mock_chain)

        self.assertTrue(mock_chain in mock_model)

    def test_get_coordinates(self):
        mock_atom_1 = Atom('N', 1, None, 1.00, 0.00, 0.00, 1.00, 7.9, 'N', None)
        mock_atom_2 = Atom('O', 1, None, 0.00, 1.00, 0.00, 1.00, 7.9, 'O', None)

        mock_res = Residue('ALA', 1, None)
        mock_res.add_atom(mock_atom_1.copy())
        mock_res.add_atom(mock_atom_2.copy())

        mock_chain_1 = Chain('A')
        mock_chain_1.add_residue(mock_res.copy())
        mock_chain_1.add_residue(mock_res.copy())

        mock_model = Model()

        mock_model.add_chain(mock_chain_1.copy())
        mock_model.add_chain(mock_chain_1.copy())

        np.testing.assert_array_equal(mock_model.get_coordinates(), np.array([[1.00, 0.00, 0.00],
                                                                              [0.00, 1.00, 0.00],
                                                                              [1.00, 0.00, 0.00],
                                                                              [0.00, 1.00, 0.00],
                                                                              [1.00, 0.00, 0.00],
                                                                              [0.00, 1.00, 0.00],
                                                                              [1.00, 0.00, 0.00],
                                                                              [0.00, 1.00, 0.00]]))


class TestStructure(TestCase):
    def test_to_pandas_models(self):
        structure = Structure()

        model1 = Model(1)
        model2 = Model(2)

        chain1 = Chain('A')
        chain2 = Chain('A')
        residue1 = Residue('ALA', 1, None)
        residue2 = Residue('ALA', 1, None)

        atom = Atom('N', 1, None, 0.00, 0.00, 0.00, 1.00, 0.00, 'N', None)

        residue1.add_atom(atom)
        residue2.add_atom(atom.copy())

        chain1.add_residue(residue1)
        chain2.add_residue(residue2)

        model1.add_chain(chain1)
        model2.add_chain(chain2)

        structure.add_model(model1)
        structure.add_model(model2)

        test_df = pd.DataFrame([
            {'record_type': 'ATOM',
             'atom_number': 1,
             'atom_name': 'N',
             'alt_loc': None,
             'residue_name': 'ALA',
             'chain_id': 'A',
             'residue_seq_id': 1,
             'residue_insert_code': None,
             'pos_x': 0.00,
             'pos_y': 0.00,
             'pos_z': 0.00,
             'occupancy': 1.00,
             'b_factor': 0.00,
             'element': 'N',
             'charge': None,
             'model_index': 1},
            {'record_type': 'ATOM',
             'atom_number': 1,
             'atom_name': 'N',
             'alt_loc': None,
             'residue_name': 'ALA',
             'chain_id': 'A',
             'residue_seq_id': 1,
             'residue_insert_code': None,
             'pos_x': 0.00,
             'pos_y': 0.00,
             'pos_z': 0.00,
             'occupancy': 1.00,
             'b_factor': 0.00,
             'element': 'N',
             'charge': None,
             'model_index': 2}])

        pd.testing.assert_frame_equal(structure.to_pandas(), test_df)

    def test_from_pdb_model(self):
        mock_pdb = ('MODEL        1 \n'
                    'ATOM      1  N   ALA A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      2  CA  ALA A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      3  C   ALA A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      4  O   ALA A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      5  CB  ALA A   3     -11.858   2.289   4.034  1.00  0.00           C \n'
                    'ATOM      6  H1  ALA A   3     -14.406   2.827   4.625  1.00  0.00           H \n'
                    'ATOM      7  H2  ALA A   3     -14.140   1.259   4.031  1.00  0.00           H \n'
                    'ATOM      8  H3  ALA A   3     -15.041   2.374   3.119  1.00  0.00           H \n'
                    'ENDMDL ')

        structure = Structure.from_pdb(mock_pdb)

        self.assertEqual(structure[0].serial_number, 1)
        self.assertEqual(structure[0]['A'].name, 'A')
        self.assertEqual(structure[0]['A'][3]['N'].element, 'N')

    def test_from_pandas(self):
        mock_df = pd.DataFrame(
            [['ATOM', 10570, 'N',  '', 'VAL', 'C', 4, '', -15.115, 2.602, -43.993, 1.00, 47.21, 'N', ''],
             ['ATOM', 10571, 'CA', '', 'VAL', 'C', 4, '', -14.470, 1.464, -43.338, 1.00, 41.45, 'C', ''],
             ['ATOM', 10572, 'C',  '', 'VAL', 'C', 4, '', -13.540, 1.977, -42.246, 1.00, 36.39, 'C', ''],
             ['ATOM', 10586, 'N',  '', 'GLN', 'd', 5, '', -13.811, 1.612, -40.993, 1.00, 36.76, 'N', ''],
             ['ATOM', 10587, 'CA', '', 'GLN', 'd', 5, '', -13.021, 2.105, -39.861, 1.00, 35.56, 'C', ''],
             ['ATOM', 10588, 'C',  '', 'GLN', 'd', 5, '', -12.114, 1.030, -39.260, 1.00, 30.89, 'C', ''],
             ['ATOM', 10589, 'O',  '', 'GLN', 'd', 5, '', -12.567, 0.013, -38.744, 1.00, 23.49, 'O', '']],
            columns=['record_type', 'atom_number', 'atom_name', 'alt_loc', 'residue_name', 'chain_id',
                     'residue_seq_id', 'residue_insert_code', 'pos_x', 'pos_y', 'pos_z', 'occupancy', 'b_factor',
                     'element', 'charge'])

        structure = Structure.from_pandas(mock_df)

        self.assertEqual(len(structure[0]), 2)
        self.assertEqual(structure[0]['C'][4].name, 'VAL')
        self.assertEqual(structure[0]['d'][5].name, 'GLN')

    def test_from_pandas_model(self):
        mock_df = pd.DataFrame(
            [['ATOM', 10570, 'N',  '', 'VAL', 'C', 4, '', -15.115, 2.602, -43.993, 1.00, 47.21, 'N', '', 0],
             ['ATOM', 10571, 'CA', '', 'VAL', 'C', 4, '', -14.470, 1.464, -43.338, 1.00, 41.45, 'C', '', 0],
             ['ATOM', 10572, 'C',  '', 'VAL', 'C', 4, '', -13.540, 1.977, -42.246, 1.00, 36.39, 'C', '', 0],
             ['ATOM', 10586, 'N',  '', 'GLN', 'd', 5, '', -13.811, 1.612, -40.993, 1.00, 36.76, 'N', '', 0],
             ['ATOM', 10587, 'CA', '', 'GLN', 'd', 5, '', -13.021, 2.105, -39.861, 1.00, 35.56, 'C', '', 0],
             ['ATOM', 10588, 'C',  '', 'GLN', 'd', 5, '', -12.114, 1.030, -39.260, 1.00, 30.89, 'C', '', 0],
             ['ATOM', 10589, 'O',  '', 'GLN', 'd', 5, '', -12.567, 0.013, -38.744, 1.00, 23.49, 'O', '', 0],
             ['ATOM', 10570, 'N',  '', 'VAL', 'C', 4, '', -15.115, 2.602, -43.993, 1.00, 47.21, 'N', '', 1],
             ['ATOM', 10571, 'CA', '', 'VAL', 'C', 4, '', -14.470, 1.464, -43.338, 1.00, 41.45, 'C', '', 1],
             ['ATOM', 10572, 'C',  '', 'VAL', 'C', 4, '', -13.540, 1.977, -42.246, 1.00, 36.39, 'C', '', 1],
             ['ATOM', 10586, 'N',  '', 'GLN', 'd', 5, '', -13.811, 1.612, -40.993, 1.00, 36.76, 'N', '', 1],
             ['ATOM', 10587, 'CA', '', 'GLN', 'd', 5, '', -13.021, 2.105, -39.861, 1.00, 35.56, 'C', '', 1],
             ['ATOM', 10588, 'C',  '', 'GLN', 'd', 5, '', -12.114, 1.030, -39.260, 1.00, 30.89, 'C', '', 1],
             ['ATOM', 10589, 'O',  '', 'GLN', 'd', 5, '', -12.567, 0.013, -38.744, 1.00, 23.49, 'O', '', 1]],
            columns=['record_type', 'atom_number', 'atom_name', 'alt_loc', 'residue_name', 'chain_id',
                     'residue_seq_id', 'residue_insert_code', 'pos_x', 'pos_y', 'pos_z', 'occupancy', 'b_factor',
                     'element', 'charge', 'model_index'])

        structure = Structure.from_pandas(mock_df)

        self.assertEqual(len(structure[0]), 2)
        self.assertEqual(len(structure[0]), 2)
        self.assertEqual(structure[0]['C'][4].name, 'VAL')
        self.assertEqual(structure[0]['d'][5].name, 'GLN')

        self.assertEqual(len(structure[1]), 2)
        self.assertEqual(len(structure[1]), 2)
        self.assertEqual(structure[1]['C'][4].name, 'VAL')
        self.assertEqual(structure[1]['d'][5].name, 'GLN')

    def test_construct_from_records(self):
        records = [
            ModelRecord(0),
            AtomRecord(10570, 'N',  None, 'VAL', 'C', 4, None, -15.115, 2.602, -43.993, 1.00, 47.21, 'N', None),
            AtomRecord(10571, 'CA', None, 'VAL', 'C', 4, None, -14.470, 1.464, -43.338, 1.00, 41.45, 'C', None),
            AtomRecord(10572, 'C',  None, 'VAL', 'C', 4, None, -13.540, 1.977, -42.246, 1.00, 36.39, 'C', None),
            AtomRecord(10586, 'N',  None, 'GLN', 'd', 5, None, -13.811, 1.612, -40.993, 1.00, 36.76, 'N', None),
            AtomRecord(10587, 'CA', None, 'GLN', 'd', 5, None, -13.021, 2.105, -39.861, 1.00, 35.56, 'C', None),
            AtomRecord(10588, 'C',  None, 'GLN', 'd', 5, None, -12.114, 1.030, -39.260, 1.00, 30.89, 'C', None),
            AtomRecord(10589, 'O',  None, 'GLN', 'd', 5, None, -12.567, 0.013, -38.744, 1.00, 23.49, 'O', None),
            EndModelRecord(),
            ModelRecord(1),
            AtomRecord(10570, 'N',  None, 'VAL', 'C', 4, None, -15.115, 2.602, -43.993, 1.00, 47.21, 'N', None),
            AtomRecord(10571, 'CA', None, 'VAL', 'C', 4, None, -14.470, 1.464, -43.338, 1.00, 41.45, 'C', None),
            AtomRecord(10572, 'C',  None, 'VAL', 'C', 4, None, -13.540, 1.977, -42.246, 1.00, 36.39, 'C', None),
            AtomRecord(10586, 'N',  None, 'GLN', 'd', 5, None, -13.811, 1.612, -40.993, 1.00, 36.76, 'N', None),
            AtomRecord(10587, 'CA', None, 'GLN', 'd', 5, None, -13.021, 2.105, -39.861, 1.00, 35.56, 'C', None),
            AtomRecord(10588, 'C',  None, 'GLN', 'd', 5, None, -12.114, 1.030, -39.260, 1.00, 30.89, 'C', None),
            AtomRecord(10589, 'O',  None, 'GLN', 'd', 5, None, -12.567, 0.013, -38.744, 1.00, 23.49, 'O', None),
            EndModelRecord(),
        ]

        structure = Structure._construct_from_records(records)

        self.assertEqual(len(structure[0]), 2)
        self.assertEqual(len(structure[0]), 2)
        self.assertEqual(structure[0]['C'][4].name, 'VAL')
        self.assertEqual(structure[0]['d'][5].name, 'GLN')

        self.assertEqual(len(structure[1]), 2)
        self.assertEqual(len(structure[1]), 2)
        self.assertEqual(structure[1]['C'][4].name, 'VAL')
        self.assertEqual(structure[1]['d'][5].name, 'GLN')

    def test_split_states(self):
        mock_pdb = ('ATOM      1  N  AGLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      2  N  BGLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      3  CA AGLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      4  CA BGLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      5  C  AGLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      6  C  BGLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      7  O  AGLY A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      8  O  BGLY A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      10 N  AGLY A   4     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      11 N  BGLY A   4     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      12 CA AGLY A   4     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      13 CA BGLY A   4     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      14 C  AGLY A   4     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      15 C  BGLY A   4     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      16 O  AGLY A   4     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      17 O  BGLY A   4     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      18 N   GLY A   5     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      19 CA  GLY A   5     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      20 C   GLY A   5     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      21 O   GLY A   5     -13.721   1.351   1.280  1.00  0.00           O \n')

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=StructureConstructionWarning)
            structure = Structure.from_pdb(mock_pdb)

        structure.split_states()

        self.assertEqual(len(structure), 2)

    def test_split_states_all_combinations(self):
        mock_pdb = ('ATOM      1  N  AGLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      2  N  BGLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      3  CA AGLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      4  CA BGLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      5  C  AGLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      6  C  BGLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      7  O  AGLY A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      8  O  BGLY A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      10 N  AGLY A   4     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      11 N  BGLY A   4     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      12 CA AGLY A   4     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      13 CA BGLY A   4     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      14 C  AGLY A   4     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      15 C  BGLY A   4     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      16 O  AGLY A   4     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      17 O  BGLY A   4     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      18 N   GLY A   5     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      19 CA  GLY A   5     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      20 C   GLY A   5     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      21 O   GLY A   5     -13.721   1.351   1.280  1.00  0.00           O \n')

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=StructureConstructionWarning)
            structure = Structure.from_pdb(mock_pdb)

        structure.split_states(all_combinations=True)

        self.assertEqual(len(structure), 4)

    def test_split_states_different_chains(self):
        mock_pdb = ('ATOM      1  N  AGLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      2  N  BGLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      3  CA AGLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      4  CA BGLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      5  C  AGLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      6  C  BGLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      7  O  AGLY A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      8  O  BGLY A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      10 N  AGLY B   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      11 N  BGLY B   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      12 CA AGLY B   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      13 CA BGLY B   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      14 C  AGLY B   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      15 C  BGLY B   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      16 O  AGLY B   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      17 O  BGLY B   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      18 N   GLY B   4     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      19 CA  GLY B   4     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      20 C   GLY B   4     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      21 O   GLY B   4     -13.721   1.351   1.280  1.00  0.00           O \n')

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=StructureConstructionWarning)
            structure = Structure.from_pdb(mock_pdb)

        structure.split_states()

        self.assertEqual(len(structure), 2)

    def test_split_states_8ecq(self):
        '''pdb was causing issues splitting into three states instead of two as expected.'''
        with open('tests/data/8ecq.pdb', 'r') as fh:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', category=StructureConstructionWarning)
                structure = Structure.from_pdb(fh.read())

        self.assertEqual(len(structure), 1)
        structure.split_states()
        self.assertEqual(len(structure), 2)

    def test_dehydrate(self):
        mock_pdb = ('ATOM      1  N   ALA A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
                    'ATOM      2  CA  ALA A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
                    'ATOM      3  C   ALA A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
                    'ATOM      4  O   ALA A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
                    'ATOM      5  CB  ALA A   3     -11.858   2.289   4.034  1.00  0.00           C \n'
                    'HETATM    6  O   HOH A   4     -15.330  29.061 101.651  1.00 24.28           O \n')

        structure = Structure.from_pdb(mock_pdb)

        self.assertTrue(structure[0]['A'][4]['O'].het_atom)
        self.assertEqual(len(structure[0]['A']), 2)

        structure.dehydrate()

        self.assertEqual(len(structure[0]['A']), 1)

    def test_in(self):
        mock_structure = Structure()
        mock_model = Model()

        mock_structure.add_model(mock_model)

        self.assertTrue(mock_model in mock_structure)

        mock_model_2 = Model()
        mock_model_2.add_chain(Chain('A'))

        self.assertFalse(mock_model_2 in mock_structure)

    def test_get_coordinates(self):
        mock_atom_1 = Atom('N', 1, None, 1.00, 0.00, 0.00, 1.00, 7.9, 'N', None)
        mock_atom_2 = Atom('O', 1, None, 0.00, 1.00, 0.00, 1.00, 7.9, 'O', None)

        mock_res = Residue('ALA', 1, None)
        mock_res.add_atom(mock_atom_1.copy())
        mock_res.add_atom(mock_atom_2.copy())

        mock_chain = Chain('A')
        mock_chain.add_residue(mock_res.copy())
        mock_chain.add_residue(mock_res.copy())

        mock_model = Model()

        mock_model.add_chain(mock_chain.copy())
        mock_model.add_chain(mock_chain.copy())

        mock_structure = Structure()

        mock_structure.add_model(mock_model.copy())
        mock_structure.add_model(mock_model.copy())

        np.testing.assert_array_equal(mock_structure.get_coordinates(), np.array([[1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00],
                                                                                  [1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00],
                                                                                  [1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00],
                                                                                  [1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00],
                                                                                  [1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00],
                                                                                  [1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00],
                                                                                  [1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00],
                                                                                  [1.00, 0.00, 0.00],
                                                                                  [0.00, 1.00, 0.00]]))
