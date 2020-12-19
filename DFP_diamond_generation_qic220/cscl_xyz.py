#!/usr/bin/env python3
# 02_multicomponent.py: Free energy of multicomponent Lennard-Jones system.

import inspect
import os
import sys
import numpy

# Locate and import PyPhase.

#sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
#    os.path.realpath(inspect.getfile(inspect.currentframe()))))))
sys.path.insert(0,'/raid/codes/colloid/rum216/Package/pyphase')

import pyphase
import pyphase_tools._test_utilities

pwd = '/raid/codes/colloid/rum216/Package/pyphase_back/doc/examples/structures/lammps/'

pyphase.lammps.set_lammps_global("/raid/codes/colloid/rum216/Package/lammps/src/lmp_serial")

# Make a 250-atom CuAu cuboid with density 1
#cell = pyphase.crystal.CellTools.scale_by(pyphase.crystal.CellTools.tile(pyphase.crystal.CellCodecs.read_lammps(pwd+'Cu3Au.dat'), (5, 5, 5)), 2.0)
cell = pyphase.crystal.CellTools.scale_by(pyphase.crystal.CellTools.tile(pyphase.structure.Primitive.CsCl(), (3, 3, 3)), 1.0)

pyphase.crystal.CellCodecs.write_xyz(cell,'cscl.xyz')
