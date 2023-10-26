from unittest import TestCase

import numpy as np

from python_pdb.aligners import align_entities, align_pandas_structure, align_sequences
from python_pdb.parsers import parse_pdb, parse_pdb_to_pandas


class TestAlignSequence(TestCase):
    def test(self):
        alignment, score = align_sequences('GCATGCG', 'GATTACA')

        self.assertEqual(alignment, [('G', 'G'),
                                     ('C', '-'),
                                     ('A', 'A'),
                                     ('T', 'T'),
                                     ('G', 'T'),
                                     ('-', 'A'),
                                     ('C', 'C'),
                                     ('G', 'A')])

        self.assertEqual(score, 0.0)


class TestAlign(TestCase):

    def test(self):
        structure_1 = parse_pdb(
            ('ATOM      1  N   ALA A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
             'ATOM      2  CA  ALA A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
             'ATOM      3  C   ALA A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
             'ATOM      4  O   ALA A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
             'ATOM      5  CB  ALA A   3     -11.858   2.289   4.034  1.00  0.00           C \n'
             'ATOM      6  H1  ALA A   3     -14.406   2.827   4.625  1.00  0.00           H \n'
             'ATOM      7  H2  ALA A   3     -14.140   1.259   4.031  1.00  0.00           H \n'
             'ATOM      8  H3  ALA A   3     -15.041   2.374   3.119  1.00  0.00           H \n')
        )

        structure_2 = parse_pdb(
            ('ATOM      1  N   ALA A   3    -114.239  -8.261  3.769  1.00  0.00           N \n'
             'ATOM      2  CA  ALA A   3    -112.988  -8.726  3.102  1.00  0.00           C \n'
             'ATOM      3  C   ALA A   3    -112.843  -8.059  1.732  1.00  0.00           C \n'
             'ATOM      4  O   ALA A   3    -113.721  -9.351  1.280  1.00  0.00           O \n'
             'ATOM      5  CB  ALA A   3    -111.858  -8.289  4.034  1.00  0.00           C \n'
             'ATOM      6  H1  ALA A   3    -114.406  -8.827  4.625  1.00  0.00           H \n'
             'ATOM      7  H2  ALA A   3    -114.140  -9.259  4.031  1.00  0.00           H \n'
             'ATOM      8  H3  ALA A   3    -115.041  -8.374  3.119  1.00  0.00           H \n')
        )

        new_structure = align_entities(structure_1.get_coordinates(), structure_2.get_coordinates(), structure_1)

        np.testing.assert_array_almost_equal(new_structure.get_coordinates(),
                                             np.array([[-114.23727413,   -8.65418163,    3.78295986],
                                                       [-113.03039812,   -8.00613587,    3.19168835],
                                                       [-112.86117442,   -8.44017914,    1.73368157],
                                                       [-113.6961671,   -9.12785028,    1.18019476],
                                                       [-111.85988525,   -8.5041603,    4.03925208],
                                                       [-114.42926596,   -8.24116603,    4.71759254],
                                                       [-114.06723735,   -9.67554587,    3.88470255],
                                                       [-115.05459768,   -8.49678087,    3.16192828]]))


class TestAlignPandas(TestCase):

    def test(self):
        df_1 = parse_pdb_to_pandas(
            ('ATOM      1  N   ALA A   3     -14.239   2.261   3.769  1.00  0.00           N \n'
             'ATOM      2  CA  ALA A   3     -12.988   2.726   3.102  1.00  0.00           C \n'
             'ATOM      3  C   ALA A   3     -12.843   2.059   1.732  1.00  0.00           C \n'
             'ATOM      4  O   ALA A   3     -13.721   1.351   1.280  1.00  0.00           O \n'
             'ATOM      5  CB  ALA A   3     -11.858   2.289   4.034  1.00  0.00           C \n'
             'ATOM      6  H1  ALA A   3     -14.406   2.827   4.625  1.00  0.00           H \n'
             'ATOM      7  H2  ALA A   3     -14.140   1.259   4.031  1.00  0.00           H \n'
             'ATOM      8  H3  ALA A   3     -15.041   2.374   3.119  1.00  0.00           H \n')
        )

        df_2 = parse_pdb_to_pandas(
            ('ATOM      1  N   ALA A   3    -114.239  -8.261  3.769  1.00  0.00           N \n'
             'ATOM      2  CA  ALA A   3    -112.988  -8.726  3.102  1.00  0.00           C \n'
             'ATOM      3  C   ALA A   3    -112.843  -8.059  1.732  1.00  0.00           C \n'
             'ATOM      4  O   ALA A   3    -113.721  -9.351  1.280  1.00  0.00           O \n'
             'ATOM      5  CB  ALA A   3    -111.858  -8.289  4.034  1.00  0.00           C \n'
             'ATOM      6  H1  ALA A   3    -114.406  -8.827  4.625  1.00  0.00           H \n'
             'ATOM      7  H2  ALA A   3    -114.140  -9.259  4.031  1.00  0.00           H \n'
             'ATOM      8  H3  ALA A   3    -115.041  -8.374  3.119  1.00  0.00           H \n')
        )

        new_df = align_pandas_structure(df_1[['pos_x', 'pos_y', 'pos_z']].to_numpy(),
                                        df_2[['pos_x', 'pos_y', 'pos_z']].to_numpy(),
                                        df_1)

        np.testing.assert_array_almost_equal(new_df[['pos_x', 'pos_y', 'pos_z']].to_numpy(),
                                             np.array([[-114.23727413,   -8.65418163,    3.78295986],
                                                       [-113.03039812,   -8.00613587,    3.19168835],
                                                       [-112.86117442,   -8.44017914,    1.73368157],
                                                       [-113.6961671,   -9.12785028,    1.18019476],
                                                       [-111.85988525,   -8.5041603,    4.03925208],
                                                       [-114.42926596,   -8.24116603,    4.71759254],
                                                       [-114.06723735,   -9.67554587,    3.88470255],
                                                       [-115.05459768,   -8.49678087,    3.16192828]]))
