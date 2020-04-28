from pypif.obj.common import Property, Scalar

from .base import DFTParser, Value_if_true, InvalidIngesterException
import os
from pypif.obj.common.value import Value

from ase import Atoms

class GpawParser(DFTParser):
    '''
    Parser for GPAW calculations
    '''
    def __init__(self, files):
        super(GpawParser, self).__init__(files)
        self.settings = {}
        self.all_parsed_data = {}

    # Look for appropriate files
        self.outputf = None
        for f in self._files:
            try:
                if self._get_line('gpaw', f, return_string=False, case_sens=False):
                    if self.outputf is not None:
                        raise InvalidIngesterException('More than one output file!')
                    self.outputf = f
            except UnicodeDecodeError as e:
                pass

        if self.outputf is None:
            raise InvalidIngesterException('Failed to find output file')

     # Reads each output file into a text file and writes it to a ase database

    def get_name(self): return "GPAW"

    def _get_line(self, search_string, search_file, return_string=True, case_sens=True):
        '''Return the first line containing a set of strings in a file.

        If return_string is False, we just return whether such a line
        was found. If case_sens is False, the search is case
        insensitive.

        '''
        if os.path.isfile(search_file):
            # if single search string
            if type(search_string) == type(''): search_string = [search_string]
            # if case insensitive, convert everything to lowercase
            if not case_sens: search_string = [i.lower() for i in search_string]
            with open(search_file) as fp:
                # search for the strings line by line
                for line in fp:
                    query_line = line if case_sens else line.lower()
                    if all([i in query_line for i in search_string]):
                        return line if return_string else True
                if return_string:
                    raise Exception('%s not found in %s'%(' & '.join(search_string), search_file))
                else: return False
        else: raise Exception('%s file does not exist'%search_file)
