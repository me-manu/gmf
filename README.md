	#################################
	#       gmf                     #
	#################################
	# python package to compute     #
	# the galactic magnetic field   #
	#                               #
	# Author:                       #
	# -------                       #
	# Manuel Meyer                  #
	# manuel.meyer@fysik.su.se      #
	#################################

1. Introduction
---------------

The python scripts included in this package can be used
to compute the the regular component of Galactic magnetic field model of Jansson and Farrar (2012).

If you use any of the packages, please give a reference to the 
following papers:

- Jansson & Farrar, 2012, The Astrophysical Journal, Volume 757, Issue 1, article id. 14, 13 pp.
- Pshirkov et al., 2011, The Astrophysical Journal, Volume 738, Issue 2, article id. 192, 14 pp.
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
- numpy
- scipy


3. Package contents:
--------------------
- README.md: the file you are currently reading
- __init__.py: init the python packages
- gmf.py	: functions to calculate the B-field components
- trafo.py: functions to compute the transformation from galactocentric cylindrical 
- coordinates to heliocentric spherical (l,b) coordinates
- ne2001.py: functions to compute electron density from fortran library

4. License
----------
gmf is distributed under the modified BSD License.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.
- Neither the name of the gmf developers  nor the
names of its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE GMF DEVELOPERS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
