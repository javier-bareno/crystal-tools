# crystal-tools

A [crystal structure](https://en.wikipedia.org/wiki/Crystal_structure) is defined by a [Bravais lattice](https://en.wikipedia.org/wiki/Bravais_lattice) (three vectors that describe the translation simmetry) and a basis (the set of repeating atoms). This project provides a series of Python2 modules to work with crystal structures.

## src/Base.py

This is an older version of the [transform](https://github.com/javier-bareno/transform) module. It defines the Base class, which provides the basic functionality to build Bravais lattices, changes coordinates, and calculate distances and angles between directions and crystal planes.

## src/Structure.py

This module builds on src/Base.py to define the Structure class. Structure implements a crystal structure: a basis (an instance of Base) and a set of atoms specified by symbol and coordinates.

### I/O

* **Structure.read** Create Structure instance from a text file representation.
* **Structure.read_cell** Create Structure instance from a [Carine](http://carine.crystallography.pagesperso-orange.fr) cell definition ascii file.
* **Structure.write** Write a text file representation of Structure instance.
* **Structure.write_xyz** Write an xyz file to view with e.g., [Jmol](http://jmol.sourceforge.net)
* **Structure.write_pov** Write an #include file for a [POV-ray](http://www.povray.org) scene.

### Transformations
* **Structure.new_basis** Creates an equivalent structure with new basis vectors (need to be compatible).
* **Structure.extend** Extends unit cell repeating it along unit vector, creating a supercell; e.g., 2x2x4.
* **Structure.translate** Displaces origin of basis set within unit cell; e.g., to change a surface termination.
* **Structure.sort** Sorts list of basis atoms by distance along a given plane normal.
* **Strucure.cluster** Returns a list of structures, each containing a subset of the atoms in the original structure. Each subset is composed of coplanar atoms, of specified plane normal. Useful to find plane stacking of arbitrary structures.

## src/XRD.py

Defines the XRD_pattern class, which provides methods to calculate XRD peak positions, store a list of experimental hkl indices and positions, and obtain lattice parametters (using scipy.optimize.fmin_powell or scipy.optimize.fmin).

## examples/

Contains a few examples of code to create crsytal strcutures, change bases, and split into plane stacks.
