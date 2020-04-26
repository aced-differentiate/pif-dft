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
