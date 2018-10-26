#!/usr/bin/python
import sys, os
from optparse import OptionParser

# set sys.path as shown for the package folder with __init__.py file (modformPro) from the current location
sys.path.append('.')

from mainsolver import MainSolver
from specifications import ModformProSpecifications

# Command line options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parser = OptionParser()

parser.add_option("-f", "--fname_spec", dest="fname_spec", default="", help="File having all the specifications for the package")

(cmd_options,args) = parser.parse_args()

spec = ModformProSpecifications()
spec.read_specifications_from_file(cmd_options.fname_spec)

main_solver                 = MainSolver(spec)
main_solver.run()

'''
Assumptions:
1. Ionization efficiency for peptides are identical
2. Peak area measurments have error bound of 2%, we set NOISE_STDEV=2 (being generous about the error) for intact mass spectrometry as well for bottom-up/MRM experiments.

'''
