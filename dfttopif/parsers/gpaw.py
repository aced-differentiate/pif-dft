import numpy as np

from pypif.obj.common import Property, Scalar

from .base import DFTParser, Value_if_true, InvalidIngesterException
import os
from pypif.obj.common.value import Value

from dftparse.gpaw.stdout_parser import GpawStdOutputParser
from ase import Atoms
from ase.io import read
from ase.io.ulm import InvalidULMFileError
from ase.calculators.calculator import PropertyNotImplementedError
import ase.db

class GpawParser(DFTParser):
    '''
    Parser for GPAW calculations
    '''
    def __init__(self, files):
        super(GpawParser, self).__init__(files)
        self.computational_data = {}
        self.output_traj = self._find_output("traj")
        self.output_txt = self._find_output("txt")
        parser = GpawStdOutputParser() 


        if self.output_traj is not None:
            self.atoms = read(self.output_traj)
        elif self.output_txt is not None:
            self.atoms = read(self.output_txt)
        else:
            raise InvalidIngesterException('Failed to find output file')


        if self.output_txt is not None:
            with open(self.output_txt,"r") as f:
                for line in parser.parse(f.readlines()):
                    self.computational_data.update(line)


        # Use ase db functionality to read in data to temporary ase db
        def _write_temp_asedb():
            '''

            Reads designated output file and writes it to a temporary ase db.
            Checks to see if cit_temp_db.db already exists, and overwrites it if it does.

            '''
            if os.path.exists('cit_temp_db.db'):
                os.remove('cit_temp_db.db')
            temp_db = ase.db.connect('cit_temp_db.db')
            temp_db.write(self.atoms)
            return temp_db

        self.temp_db = _write_temp_asedb()
        self.settings = self.temp_db.get(id=1).calculator_parameters


        def _get_mode():
            '''Determine calculation mode used.

            Default value pulled from https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html on May 18, 2020


            Returns:
               String, Possibilities are 'fd' (finite difference real space grid), 'lcao', 'pw' (plane-wave)

            '''
            try:
                mode=self.settings['mode']['name']
            except KeyError:
                mode = 'fd'
            return mode


        # Get mode for calculation
        self.calc_mode = _get_mode()

        # Look for ase traj files
    def _find_output(self,fmt):
        '''Searches for GPAW readable traj file and returns name '''
        out_file = None
        for f in self._files:
            if os.path.basename(f) == "output." + fmt:
                out_file = f
            else:
                if os.path.basename(f).split('.')[-1] == fmt:
                    try:
                        test_out_file = read(f)
                        if out_file is not None:
                            raise InvalidIngesterException('Found more than one valid {} file'.format(fmt))
                        out_file = f
                    except InvalidULMFileError:
                        pass
                    except OSError:
                        pass
        return out_file

    def get_name(self): return "GPAW"

    def get_version_number(self):
        '''Determine the GPAW version number from the output'''
        if self.output_txt is not None:
            return self.computational_data["gpaw_version"]
        else:
            return None


    def get_setting_functions(self):
        base_settings = super(GpawParser, self).get_setting_functions()
        base_settings["Grid Spacing"] = "get_grid_spacing"
        base_settings["Calculation Mode"] = "get_calc_mode"
        return base_settings

    def get_result_functions(self):
        base_results = super(GpawParser, self).get_result_functions()
        return base_results


    def get_total_energy(self):
        '''Determine total energy from the temporary ase db'''
        ener = self.temp_db.get(id=1).energy
        return Property(scalars=[Scalar(value=ener)], units='eV')


    def get_calc_mode(self):
        '''Takes string of calc mode and converts it to a pif Value'''
        m = self.calc_mode
        if m == 'fd':
            return Value(scalars=[Scalar(value='Finite Difference')])
        elif m == 'lcao':
            return Value(scalars=[Scalar(value='LCAO')])
        elif m == 'pw':
            return Value(scalars=[Scalar(value='Plane-Wave')])

    def get_grid_spacing(self):
        '''Determine grid spacing from the temporary ase db

        Default value pulled from https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html on May 18, 2020

        Returns: Value, None if the calculator mode was not 'fd'

        '''
        if self.calc_mode == 'fd':
            try:
                h = self.settings['h']
            except KeyError:
                h = 0.2
            return Value(scalars=[Scalar(value=h)],units='angstrom')
        return None

    def get_xc_functional(self):
        '''Determine the xc functional from the temporary ase db

        Default value pulled from https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html on May 18, 2020

        '''
        try:
            xc = self.settings['xc']
        except KeyError:
            xc = 'LDA'
        return Value(scalars=[Scalar(value=xc)])

    def get_cutoff_energy(self):
        '''Determine cutoff energy if calc mode is plane wave, else returns None

        Default value pulled from https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html on May 18, 2020

        '''
        if self.calc_mode == 'pw':
            try:
                ecut = self.settings['mode']['ecut']
            except KeyError:
                ecut = 340.0
            return Value(scalars=[Scalar(value=ecut)],units='eV')
        return None

    def get_KPPRA(self):
        ''' Determine the no. of k-points in the BZ times the no. of atoms

        Default value pulled from https://wiki.fysik.dtu.dk/gpaw/documentation/manual.html on May 18, 2020

        '''
        try:
            kp = self.settings['kpts']
        except KeyError:
            kp = np.array([0.])
        natoms = len(self.atoms)
        return Value(scalars=[Scalar(value=kp[0]*natoms)])

    def get_total_magnetization(self):
        try:
            spin = self.settings['spinpol']
        except KeyError:
            return None
        if spin:
            try:
                tot_mag = self.atoms.get_magnetic_moment()
            except KeyError:
                return None
            return Property(scalars=[Scalar(value=tot_mag)],units='Bohr')
        return None

    def get_output_structure(self):
        ''' Returns ase Atoms object'''
        return self.atoms


    def get_final_volume(self):
        vol = self.atoms.get_cell().volume
        return Property(scalars=[Scalar(value=vol)],units='angstrom^3')


    def get_forces(self):
        try:
            all_forces = self.atoms.get_forces()
            wrapped = [[Scalar(value=x) for x in y] for y in all_forces]
            return Property(vectors=wrapped,units='eV/angstrom')
        except PropertyNotImplementedError:
            return None


    def get_max_force(self):
        ''' Returns maximum force in the structure'''
        max_f = self.temp_db.get(id=1).fmax
        return Property(scalars=[Scalar(value=max_f)],units='eV/angstrom')

    @Value_if_true
    def is_relaxed(self):
        '''Determine if relaxation run by checking if more than one image present'''
        if self.output_traj is not None:
            return len(read(self.output_traj,index=':')) > 1
        else:
            return len(read(self.output_txt,index=':')) > 1

# Begin function placeholders

    @Value_if_true
    def uses_SOC(self): return None

    def get_pp_name(self): return None

    def get_U_settings(self): return None

    def get_vdW_settings(self): return None

    def is_converged(self): return None

    def get_band_gap(self): return None

    def get_pressure(self): return None

    def get_dos(self): return None

    def get_density(self): return None

    def get_stresses(self): return None

    def get_initial_volume(self): return None


# End function placeholders

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
