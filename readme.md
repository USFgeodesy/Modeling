<h1>README</h1>

<h2>okadaPY</h2>

collection of python scripts for using Okada's DC3D0 fortran code
for modeling of fault dislocations

Description of programs:

**/okadawrapper1 :** directory with okada wrappers from tbenthompson: https://github.com/tbenthompson/okada_wrapper.git

**okadaPY.py :** main program, takes in fault file, a gps file, and a slip file

**readData.py:** functions for reading in the data files including the data files specified on the command line

**okada_functions.py:** functions for doing the dislocation calculations including the DCD0 fortran code, and the rotations

**outputData.py:** writes the computed surface displacements to files

**okada.py:** wrapper for the DC3D0

**plotData.py:** simple plotting routines for the displacments at gps locations and on a grid

**plotPatch.py:** higher quality plot of just GPS displacments, and slip patterns on fault



