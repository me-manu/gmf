	#################################
	#	gmf			#
	#################################
	# python package to compute 	#
	# the galactic magnetic field	#
	# from the Jansson & Farrar	#
	# model (2012)			#
	#				#
	# Author:			#
	# -------			#
	# Manuel Meyer			#
	# Version: 0.02 07/19/2013	#
	# manuel.meyer@desy.de		#
	#################################

1. Introduction
---------------

The python scripts included in this package can be used
to compute the the regular component of Galactic magnetic field model of Jansson and Farrar (2012).

If you use any of the packages, please give a reference to the 
following papers:

- Jansson & Farrar, 2012, The Astrophysical Journal, Volume 757, Issue 1, article id. 14, 13 pp.
- Horns et al., 2012, Physical Review D, vol. 86, Issue 7
- Meyer et al., 2013, Physical Review D, vol. 87, Issue 3, id. 035027

2. Prerequisites
----------------

The scripts require further packages written available at https://github.com/me-manu/
namely the package eblstud, and a running version of the modified
NE2001 code (Cordes & Lazio, 2001), also available at the above repository.
Download the packages and the add the paths to your PYTHONPATH variable.
Set the NE2001_PATH variable to the file path of the NE2001 installation,
e.g.,
export NE2001_PATH="/path/to/ne2001/implementation/NE2001/bin.NE2001/"

You will also need the following standard python packages:
numpy
scipy
matplotlib


3. Package contents:
--------------------
README.md		-	the file you are currently reading
__init__.py		-	init the python packages
gmf.py			-	functions to calculate the B-field components
trafo.py		-	functions to compute the transformation from galactocentric cylindrical 
				coordinates to heliocentric spherical (l,b) coordinates
ne2001.py		-	functions to compute electron density from fortran library
example.py		- 	example file to calculate the conversion (not yet added)
