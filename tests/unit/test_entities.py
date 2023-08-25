from unittest import TestCase

import pandas as pd

from python_pdb.entities import Atom, Chain, Model, Residue, Structure


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

        self.assertEqual(len(structure[0].chains), 2)
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

        self.assertEqual(len(structure[0].chains), 2)
        self.assertEqual(len(structure[0].chains), 2)
        self.assertEqual(structure[0]['C'][4].name, 'VAL')
        self.assertEqual(structure[0]['d'][5].name, 'GLN')

        self.assertEqual(len(structure[1].chains), 2)
        self.assertEqual(len(structure[1].chains), 2)
        self.assertEqual(structure[1]['C'][4].name, 'VAL')
        self.assertEqual(structure[1]['d'][5].name, 'GLN')
