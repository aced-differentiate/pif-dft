import unittest
from dfttopif.parsers.gpaw import GpawParser
from ..test_pif import unpack_example, delete_example
from pypif.obj.common.value import Value
import os
import shutil


class TestGpawParser(unittest.TestCase):

    def get_parser(self,name):
        '''Get a GpawParser for a certain test'''
        unpack_example(os.path.join('examples', 'gpaw', name+'.tar.gz'))
        return GpawParser.generate_from_directory(name)

    def test_H2(self):
        """Test that a static calculation is parseable"""

        # Parse the results
        parser = self.get_parser("H2")

        # Test the settings
        self.assertEquals("GPAW", parser.get_name())
        self.assertEquals("19.8.1", parser.get_version_number())

        grid_spacing = parser.get_grid_spacing()
        self.assertEquals(0.18, grid_spacing.scalars[0].value)
        self.assertEquals("angstrom", grid_spacing.units)

        xc = parser.get_xc_functional()
        self.assertEquals("BEEF-vdW",xc.scalars[0].value)

        calc_mode = parser.get_calc_mode()
        self.assertEquals("Finite Difference",calc_mode.scalars[0].value)

        smear_func = parser.get_smearing_function()
        self.assertEquals("fermi-dirac",smear_func.scalars[0].value)

        smear_width = parser.get_smearing_width()
        self.assertEquals(0.05,smear_width.scalars[0].value)
        self.assertEquals("eV",smear_width.units)

        # Test the results
        energy = parser.get_total_energy()
        self.assertAlmostEquals(-8.010749238310888,energy.scalars[0].value)
        self.assertEquals("eV", energy.units)

        # Test the structure
        strc = parser.get_output_structure()
        self.assertEquals(["H","H"], strc.get_chemical_symbols())
        self.assertEquals("H2",parser.get_composition())
        self.assertEquals(10.0, strc.cell[0][0])

        delete_example("H2")

if __name__ == '__main__':
    unittest.main()
