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

    def test_H2_traj_and_txt(self):
        """Test parsing when both traj and txt present"""

        # Parse the results
        parser = self.get_parser("H2")

        # Test the settings
        self.assertEquals("GPAW", parser.get_name())
        self.assertEquals("19.8.1", parser.get_version_number())

        self.assertEquals("3.19.1", parser.get_ase_version())
        self.assertEquals("3.0.0", parser.get_libxc_version())
        self.assertEquals("3.7.7", parser.get_python_version())
        self.assertEquals("1.18.1", parser.get_numpy_version())
        self.assertEquals("1.4.1", parser.get_scipy_version())

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

        ncores = parser.get_ncores()
        self.assertEquals(8,ncores.scalars[0].value)

        kppra = parser.get_KPPRA()
        self.assertEquals(2.,kppra.scalars[0].value)

        # Test the results
        energy = parser.get_total_energy()
        self.assertAlmostEquals(-8.010749238310888,energy.scalars[0].value)
        self.assertEquals("eV", energy.units)

        # Test the timing
        total_timing = parser.get_total_timing()
        self.assertAlmostEquals(501.283,total_timing.scalars[0].value)
        wf_time = parser.get_wf_initialization_timing()
        self.assertIsNone(wf_time)
        force_time = parser.get_forces_timing()
        self.assertIsNone(force_time)
        lcao_time = parser.get_lcao_initialization_timing()
        self.assertAlmostEquals(1.018,lcao_time.scalars[0].value)
        scf_time = parser.get_scf_cycle_timing()
        self.assertAlmostEquals(464.621,scf_time.scalars[0].value)
        other_time = parser.get_other_timing()
        self.assertAlmostEquals(7.459,other_time.scalars[0].value)

        # Test getting memory usage
        mem = parser.get_memory_usage()
        self.assertAlmostEquals(415.17,mem.scalars[0].value)
        self.assertEquals(mem.units, "MiB")

        # Test the run date
        date = parser.get_run_date()
        self.assertEquals("Fri Jan 8 21:28:12 2021",date.scalars[0].value)

        # Test the structure
        strc = parser.get_output_structure()
        self.assertEquals(["H","H"], strc.get_chemical_symbols())
        self.assertEquals("H2",parser.get_composition())
        self.assertEquals(10.0, strc.cell[0][0])

        # Test getting the setup files
        pp = parser.get_pp_name()
        self.assertEquals("/home/azeeshan/software/gpaw-setups-0.9.20000/H.PBE.gz",pp.scalars[0])

        delete_example("H2")

    def test_H2_txt_only(self):
        """Test parsing when only the txt is present"""

        # Parse the results
        parser = self.get_parser("H2_txt_only")

        # Test the settings
        self.assertEquals("GPAW", parser.get_name())
        self.assertEquals("19.8.1", parser.get_version_number())

        self.assertEquals("3.19.1", parser.get_ase_version())
        self.assertEquals("3.0.0", parser.get_libxc_version())
        self.assertEquals("3.7.7", parser.get_python_version())
        self.assertEquals("1.18.1", parser.get_numpy_version())
        self.assertEquals("1.4.1", parser.get_scipy_version())

#        grid_spacing = parser.get_grid_spacing()
#        self.assertEquals(0.18, grid_spacing.scalars[0].value)
#        self.assertEquals("angstrom", grid_spacing.units)

        xc = parser.get_xc_functional()
        self.assertEquals("BEEF-vdW",xc.scalars[0].value)

#        calc_mode = parser.get_calc_mode()
#        self.assertEquals("Finite Difference",calc_mode.scalars[0].value)

#        smear_func = parser.get_smearing_function()
#        self.assertEquals("fermi-dirac",smear_func.scalars[0].value)

#        smear_width = parser.get_smearing_width()
#        self.assertEquals(0.05,smear_width.scalars[0].value)
#        self.assertEquals("eV",smear_width.units)

        ncores = parser.get_ncores()
        self.assertEquals(8,ncores.scalars[0].value)

        kppra = parser.get_KPPRA()
        self.assertEquals(2.,kppra.scalars[0].value)

        # Test the results
        energy = parser.get_total_energy()
        self.assertAlmostEquals(-8.010749,energy.scalars[0].value)
        self.assertEquals("eV", energy.units)

        # Test the timing
        total_timing = parser.get_total_timing()
        self.assertAlmostEquals(501.283,total_timing.scalars[0].value)
        wf_time = parser.get_wf_initialization_timing()
        self.assertIsNone(wf_time)
        force_time = parser.get_forces_timing()
        self.assertIsNone(force_time)
        lcao_time = parser.get_lcao_initialization_timing()
        self.assertAlmostEquals(1.018,lcao_time.scalars[0].value)
        scf_time = parser.get_scf_cycle_timing()
        self.assertAlmostEquals(464.621,scf_time.scalars[0].value)
        other_time = parser.get_other_timing()
        self.assertAlmostEquals(7.459,other_time.scalars[0].value)

        # Test getting memory usage
        mem = parser.get_memory_usage()
        self.assertAlmostEquals(415.17,mem.scalars[0].value)
        self.assertEquals(mem.units, "MiB")

        # Test the run date
        date = parser.get_run_date()
        self.assertEquals("Fri Jan 8 21:28:12 2021",date.scalars[0].value)

        # Test the structure
        strc = parser.get_output_structure()
        self.assertEquals(["H","H"], strc.get_chemical_symbols())
        self.assertEquals("H2",parser.get_composition())
        self.assertEquals(10.0, strc.cell[0][0])

        # Test getting the setup files
        pp = parser.get_pp_name()
        self.assertEquals("/home/azeeshan/software/gpaw-setups-0.9.20000/H.PBE.gz",pp.scalars[0])

        delete_example("H2_txt_only")

if __name__ == '__main__':
    unittest.main()
