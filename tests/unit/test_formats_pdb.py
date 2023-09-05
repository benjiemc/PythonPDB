from unittest import TestCase

from python_pdb.formats.pdb import format_model_record


class TestFormatModelRecord(TestCase):
    def test_1(self):
        self.assertEqual(format_model_record(1), 'MODEL        1')

    def test_12(self):
        self.assertEqual(format_model_record(12), 'MODEL       12')
