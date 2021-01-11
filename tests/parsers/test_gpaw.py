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

        delete_example("H2")

if __name__ == '__main__':
    unittest.main()
