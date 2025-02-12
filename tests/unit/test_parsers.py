import os
import warnings
from unittest import TestCase

import pandas as pd

from python_pdb.entities import StructureConstructionWarning
from python_pdb.parsers import parse_pdb, parse_pdb_to_pandas, stringify_structure

TEST_DATA = 'tests/data'


class TestParsePDB(TestCase):
    def test(self):
        with open(os.path.join(TEST_DATA, '6zky_mock.pdb'), 'r') as fh:
            pdb_contents = fh.read()

        structure = parse_pdb(pdb_contents)

        self.assertEqual(len(structure), 1)

        self.assertEqual(len(structure[0]), 2)

        self.assertEqual(len(structure[0]['E']), 245)
        self.assertEqual(structure[0]['E'][1].name, 'ASN')
        self.assertEqual(structure[0]['E'][1].seq_id, 1)

        self.assertEqual(len(structure[0]['E'][1]), 8)
        self.assertEqual(structure[0]['E'][1]['N'].name, 'N')
        self.assertEqual(structure[0]['E'][1]['N'].number, 1431)
        self.assertEqual(structure[0]['E'][1]['N'].alt_loc, None)
        self.assertEqual(structure[0]['E'][1]['N'].position, (-41.307, -10.041, 5.372))
        self.assertEqual(structure[0]['E'][1]['N'].occupancy, 1.0)
        self.assertEqual(structure[0]['E'][1]['N'].b_factor, 96.09)
        self.assertEqual(structure[0]['E'][1]['N'].element, 'N')
        self.assertEqual(structure[0]['E'][1]['N'].charge, None)

    def test_silent(self):
        mock_pdb = (
            'ATOM      1  N  AGLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
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
            'ATOM      21 O   GLY A   5     -13.721   1.351   1.280  1.00  0.00           O \n'
        )

        with self.assertWarns(StructureConstructionWarning):
            parse_pdb(mock_pdb)

        with warnings.catch_warnings():
            warnings.filterwarnings('error', category=StructureConstructionWarning)
            parse_pdb(mock_pdb, silent=True)


class TestParsePDBToPandas(TestCase):
    def test(self):
        mock_pdb = (
            'ATOM      1  N   GLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
            'ATOM      2  CA  GLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
            'ATOM      3  C   GLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
            'ATOM      4  O   GLY A   3     -13.721   1.351   1.280  1.00  0.00           O '
        )

        test_df = parse_pdb_to_pandas(mock_pdb)

        mock_df = pd.DataFrame(
            [
                {
                    'record_type': 'ATOM',
                    'atom_number': 1,
                    'atom_name': 'N',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -14.239,
                    'pos_y': 2.261,
                    'pos_z': 3.769,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'N',
                    'charge': None,
                },
                {
                    'record_type': 'ATOM',
                    'atom_number': 2,
                    'atom_name': 'CA',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -12.988,
                    'pos_y': 2.726,
                    'pos_z': 3.102,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'C',
                    'charge': None,
                },
                {
                    'record_type': 'ATOM',
                    'atom_number': 3,
                    'atom_name': 'C',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -12.843,
                    'pos_y': 2.059,
                    'pos_z': 1.732,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'C',
                    'charge': None,
                },
                {
                    'record_type': 'ATOM',
                    'atom_number': 4,
                    'atom_name': 'O',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -13.721,
                    'pos_y': 1.351,
                    'pos_z': 1.280,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'O',
                    'charge': None,
                },
            ]
        )

        pd.testing.assert_frame_equal(test_df, mock_df)

    def test_with_models(self):
        mock_pdb = (
            'MODEL        1 \n'
            'ATOM      1  N   GLY A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
            'ATOM      2  CA  GLY A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
            'ATOM      3  C   GLY A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
            'ATOM      4  O   GLY A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
            'ENDMDL '
        )

        test_df = parse_pdb_to_pandas(mock_pdb)

        mock_df = pd.DataFrame(
            [
                {
                    'record_type': 'ATOM',
                    'atom_number': 1,
                    'atom_name': 'N',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -14.239,
                    'pos_y': 2.261,
                    'pos_z': 3.769,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'N',
                    'charge': None,
                    'model_index': 1,
                },
                {
                    'record_type': 'ATOM',
                    'atom_number': 2,
                    'atom_name': 'CA',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -12.988,
                    'pos_y': 2.726,
                    'pos_z': 3.102,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'C',
                    'charge': None,
                    'model_index': 1,
                },
                {
                    'record_type': 'ATOM',
                    'atom_number': 3,
                    'atom_name': 'C',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -12.843,
                    'pos_y': 2.059,
                    'pos_z': 1.732,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'C',
                    'charge': None,
                    'model_index': 1,
                },
                {
                    'record_type': 'ATOM',
                    'atom_number': 4,
                    'atom_name': 'O',
                    'alt_loc': None,
                    'residue_name': 'GLY',
                    'chain_id': 'A',
                    'residue_seq_id': 3,
                    'residue_insert_code': None,
                    'pos_x': -13.721,
                    'pos_y': 1.351,
                    'pos_z': 1.280,
                    'occupancy': 1.00,
                    'b_factor': 0.00,
                    'element': 'O',
                    'charge': None,
                    'model_index': 1,
                },
            ]
        )

        pd.testing.assert_frame_equal(test_df, mock_df)


class TestStringifyStructure(TestCase):
    def test_hydrogens(self):
        mock_contents = (
            'ATOM  10570  N   VAL C   4     -15.115   2.602 -43.993  1.00 47.21           N  \n'
            'ATOM  10571  CA  VAL C   4     -14.470   1.464 -43.338  1.00 41.45           C  \n'
            'ATOM  10572  C   VAL C   4     -13.540   1.977 -42.246  1.00 36.39           C  \n'
            'ATOM  10573  O   VAL C   4     -12.601   2.710 -42.537  1.00 34.14           O  \n'
            'ATOM  10574  CB  VAL C   4     -13.684   0.615 -44.374  1.00 32.54           C  \n'
            'ATOM  10575  CG1 VAL C   4     -12.700  -0.317 -43.694  1.00 30.75           C  \n'
            'ATOM  10576  CG2 VAL C   4     -14.657  -0.175 -45.232  1.00 30.79           C  \n'
            'ATOM  10577  H   VAL C   4     -14.563   3.180 -44.309  1.00 56.62           H  \n'
            'ATOM  10578  HA  VAL C   4     -15.153   0.895 -42.924  1.00 49.71           H  \n'
            'ATOM  10579  HB  VAL C   4     -13.179   1.215 -44.961  1.00 39.01           H  \n'
            'ATOM  10580 HG11 VAL C   4     -12.235  -0.823 -44.364  1.00 36.87           H  \n'
            'ATOM  10581 HG12 VAL C   4     -12.075   0.205 -43.185  1.00 36.87           H  \n'
            'ATOM  10582 HG13 VAL C   4     -13.182  -0.909 -43.112  1.00 36.87           H  \n'
            'ATOM  10583 HG21 VAL C   4     -14.161  -0.695 -45.868  1.00 36.92           H  \n'
            'ATOM  10584 HG22 VAL C   4     -15.174  -0.754 -44.667  1.00 36.92           H  \n'
            'ATOM  10585 HG23 VAL C   4     -15.237   0.438 -45.691  1.00 36.92           H  \n'
            'ATOM  10586  N   GLN C   5     -13.811   1.612 -40.993  1.00 36.76           N  \n'
            'ATOM  10587  CA  GLN C   5     -13.021   2.105 -39.861  1.00 35.56           C  \n'
            'ATOM  10588  C   GLN C   5     -12.114   1.030 -39.260  1.00 30.89           C  \n'
            'ATOM  10589  O   GLN C   5     -12.567   0.013 -38.744  1.00 23.49           O  \n'
            'ATOM  10590  CB  GLN C   5     -13.924   2.678 -38.765  1.00 37.26           C  \n'
            'ATOM  10591  CG  GLN C   5     -13.135   3.164 -37.558  1.00 58.53           C  \n'
            'ATOM  10592  CD  GLN C   5     -13.985   3.914 -36.548  1.00 69.58           C  \n'
            'ATOM  10593  OE1 GLN C   5     -14.976   3.392 -36.041  1.00 68.00           O  \n'
            'ATOM  10594  NE2 GLN C   5     -13.599   5.149 -36.253  1.00 65.40           N  \n'
            'ATOM  10595  H   GLN C   5     -14.447   1.078 -40.769  1.00 44.08           H  \n'
            'ATOM  10596  HA  GLN C   5     -12.446   2.832 -40.178  1.00 42.65           H  \n'
            'ATOM  10597  HB2 GLN C   5     -14.420   3.430 -39.124  1.00 44.68           H  \n'
            'ATOM  10598  HB3 GLN C   5     -14.537   1.988 -38.467  1.00 44.68           H  \n'
            'ATOM  10599  HG2 GLN C   5     -12.745   2.399 -37.108  1.00 70.20           H  \n'
            'ATOM  10600  HG3 GLN C   5     -12.435   3.763 -37.860  1.00 70.20           H  \n'
            'ATOM  10601 HE21 GLN C   5     -12.901   5.482 -36.628  1.00 78.45           H  \n'
            'ATOM  10602 HE22 GLN C   5     -14.047   5.615 -35.686  1.00 78.45           H  '
        )

        structure = parse_pdb(mock_contents)

        self.assertEqual(mock_contents, stringify_structure(structure))
