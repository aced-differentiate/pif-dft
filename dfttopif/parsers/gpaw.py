from pypif.obj.common import Property, Scalar

from .base import DFTParser, Value_if_true, InvalidIngesterException
import os
from pypif.obj.common.value import Value

from ase import Atoms
from ase.io import read
from ase.io.ulm import InvalidULMFileError
import ase.db

class GpawParser(DFTParser):
    '''
    Parser for GPAW calculations
    '''
    def __init__(self, files):
        super(GpawParser, self).__init__(files)
        self.settings = {}
        self.all_parsed_data = {}

    # Look for ase traj files
        def _find_traj():
            traj_file = None
            for f in self._files:
                if os.path.basename(f).split('.')[-1] == 'traj':
                    try:
                        test_traj_file = read(f, format='traj')
                        if traj_file is not None:
                            raise InvalidIngesterException('Found more than one valid traj file')
                        traj_file = f
                    except InvalidULMFileError:
                        pass
                    except OSError:
                        pass
            return traj_file

        self.outputf = _find_traj()



    # Look for appropriate txt if no traj files
        if self.outputf is None:
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

    # Use ase db functionality to read in data
    atoms = read(self.outputf)
    temp_db = ase.db.connect('temp_db.db')
    temp_db.write(atoms)

    def get_total_energy(self):
        ener = temp_db.get(id=1).energy
        return Property(scalars=[Scalar(value=ener)], units='eV')


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
