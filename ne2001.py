import numpy as np
import warnings
from subprocess import PIPE,Popen
from trafo import x_HC2GC, y_HC2GC, z_HC2GC
import logging
import pickle
import os

def density_2001(x,y,z):
    """
    Return thermal electron density in cm^-3 of NE2001 model

    Compiled versions of NE2001_dens (my code) and NE2001 are mandatory,
    they are called via subprocess.

    The paths to the defined variables are given in environment variable NE2001_PATH

    Parameters
    ----------
    Coordinates are galactocentric with sun at x=-8.5, y=0, z=0
    x: x coordinate in kpc (N-dim array or scalar)
    y: y coordinate in kpc (N-dim array or scalar)
    z: z coordinate in kpc (N-dim array or scalar)

    Returns
    -------
    thermal electron density at (x,y,z) in cm^-3
    as N-dim numpy array (for each (x,y,z)-tuple)

    Notes
    -----
    Model is described in Cordes & Lazio (2002), see http://adsabs.harvard.edu/abs/2002astro.ph..7156C
    The source code is available at http://www.astro.cornell.edu/~cordes/NE2001/
    and was compiled succesfully using g77 (3.4.6)


    Description of density_2001 fortran subroutine in original code
    ---------------------------------------------------------------
    Returns seven components of the free electron density of the
    interstellar medium at Galactic location (x,y,z).
    Calling arguments:
    input:
	 x, y, z = galactocentric location (kpc)
	 Right-handed coordinate system
	 x is in l=90 direction
	 y is in l=180 direction
	 The sun is at (x,y,z) = (0,R0,0)
    output:
      electron densities in cm^{-3}:
	 ne1:    outer, thick disk
	 ne2:    inner, thin disk (annular in form)
	 nea:    spiral arms
	 negc:   galactic center component
	 nelism: local ISM component
	 necN:   contribution from discrete 'clumps'
	 nevN:   contribution from voids
      fluctuation parameters (one for each ne component):
	 F1, F2, Fa, Fgc, Flism, FcN, FvN
      flags:
	 whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
	    wlism: 1 if x,y,z is in any of the four LISM components
	     wLDR: 1 if in LDR, 0 if not
	     wLHB: 1 if in LHB, 0 if not
	     wLSB: 1 if in LSB, 0 if not
	   wLOOPI: 1 if in LoopI, 0 if not
	 (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
	 hitclump: clump number that x,y,z is in (0 if none)
	  hitvoid: void number that x,y,z is in (0 if none)
    25 May 2002
    based on routines from TC93 and test routines from 1999-2002 by JMC.
    """

# Rotate x and y: in NE2001 x is in l=90 and y in l=180 direction
# I prefer x to be in the l = 0 and y in the l = 90 direction
# Hence: (x_NE2001, y_NE2001) = R(3\pi/2) (x,y) and z_NE2001 = z
# where R(theta) is the standard 2x2 rotation matrix of SO(2)
    x_NE2001 = x * np.cos(np.pi * 1.5) - y * np.sin(np.pi * 1.5)
    y_NE2001 = x * np.sin(np.pi * 1.5) + y * np.cos(np.pi * 1.5)
    z_NE2001 = z

    if np.isscalar(x):
	x_NE2001 = np.array([x_NE2001])
	y_NE2001 = np.array([y_NE2001])
	z_NE2001 = np.array([z_NE2001])

    ne = np.zeros(x_NE2001.shape[0])

    try:
	ne2001_path = os.environ['NE2001_PATH']
    except KeyError:
	ne2001_path = '/afs/desy.de/user/m/meyerm/projects/NE2001/bin.NE2001/'
	logging.warning("Did not find environmental variable for NE2001_PATH, using {0} instead.".format(ne2001_path))

    for i,xN in enumerate(x_NE2001):
	if np.sqrt(xN**2. + y_NE2001[i]**2.) > 15.:
	    ne[i] = 0.
	else:
	    try:
		p = Popen([ne2001_path + 'NE2001_dens','{0}'.format(xN),'{0}'.format(y_NE2001[i]),'{0}'.format(z_NE2001[i])], stdout = PIPE)
	    except OSError:
		warnings.warn('Could not run NE2001_dens, is it compiled correctly?',RuntimeWarning)
		return -1
	    out, err = p.communicate()
	    #logging.debug('NE2001 code output: {0},{1}'.format(out,err))
	    if not len(out):
		warnings.warn('No output of NE2001_dens, is it compiled correctly?',RuntimeWarning)
		return -1
# n is list with 
#      electron densities in cm^{-3}:
#	 ne1:    outer, thick disk
#	 ne2:    inner, thin disk (annular in form)
#	 nea:    spiral arms
#	 negc:   galactic center component
#	 nelism: local ISM component
#	 necN:   contribution from discrete 'clumps'
#	 nevN:   contribution from voids
	    n = np.array(map(lambda x: float(x),out.split()))
	    if n[-1] > 0.:	# we are in a void: supersedes all but clumps
		ne[i] = n[-2] + n[-1]
	    else:
		ne[i] = np.sum(n)
    return ne


def density_2001_los(s,l,b, saveoutput,d = -8.5):
    """
    Save thermal electron density in cm^-3 of NE2001 model
    for galactic coordinates l and b and an array of distances to a pickle file

    The paths to the defined variables are given in environment variable NE2001_PATH

    Parameters
    ----------
    Coordinates are galactocentric with sun at x=-8.5, y=0, z=0
    s: N-dim array with distances from the sun in kpc
    l: galactic longitude
    b: galactic latitude
    saveoutput: string with full filename to save output
    d: distance of sun along x-axis in kpc

    Returns
    -------
    thermal electron density at along line of sight in cm^-3

    Notes
    -----
    Model is described in Cordes & Lazio (2002), see http://adsabs.harvard.edu/abs/2002astro.ph..7156C
    The source code is available at http://www.astro.cornell.edu/~cordes/NE2001/
    and was compiled succesfully using g77 (3.4.6)


    Description of density_2001 fortran subroutine in original code
    ---------------------------------------------------------------
    Returns seven components of the free electron density of the
    interstellar medium at Galactic location (x,y,z).
    Calling arguments:
    input:
	 x, y, z = galactocentric location (kpc)
	 Right-handed coordinate system
	 x is in l=90 direction
	 y is in l=180 direction
	 The sun is at (x,y,z) = (0,R0,0)
    output:
      electron densities in cm^{-3}:
	 ne1:    outer, thick disk
	 ne2:    inner, thin disk (annular in form)
	 nea:    spiral arms
	 negc:   galactic center component
	 nelism: local ISM component
	 necN:   contribution from discrete 'clumps'
	 nevN:   contribution from voids
      fluctuation parameters (one for each ne component):
	 F1, F2, Fa, Fgc, Flism, FcN, FvN
      flags:
	 whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
	    wlism: 1 if x,y,z is in any of the four LISM components
	     wLDR: 1 if in LDR, 0 if not
	     wLHB: 1 if in LHB, 0 if not
	     wLSB: 1 if in LSB, 0 if not
	   wLOOPI: 1 if in LoopI, 0 if not
	 (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
	 hitclump: clump number that x,y,z is in (0 if none)
	  hitvoid: void number that x,y,z is in (0 if none)
    25 May 2002
    based on routines from TC93 and test routines from 1999-2002 by JMC.
    """

    x = x_HC2GC(s,l,b,d)
    y = y_HC2GC(s,l,b,d)
    z = z_HC2GC(s,l,b,d)

# Rotate x and y: in NE2001 x is in l=90 and y in l=180 direction
# I prefer x to be in the l = 0 and y in the l = 90 direction
# Hence: (x_NE2001, y_NE2001) = R(3\pi/2) (x,y) and z_NE2001 = z
# where R(theta) is the standard 2x2 rotation matrix of SO(2)
    x_NE2001 = x * np.cos(np.pi * 1.5) - y * np.sin(np.pi * 1.5)
    y_NE2001 = x * np.sin(np.pi * 1.5) + y * np.cos(np.pi * 1.5)
    z_NE2001 = z

    ne = np.zeros(x_NE2001.shape[0])

    try:
	ne2001_path = os.environ['NE2001_PATH']
    except KeyError:
	ne2001_path = '/afs/desy.de/user/m/meyerm/projects/NE2001/bin.NE2001/'
	logging.warning("Did not find environmental variable for NE2001_PATH, using {0} instead.".format(ne2001_path))

    m = np.sqrt(x_NE2001**2. + y_NE2001**2.) < 15.	# mask, if x^2 + y^2 > 15 model is zero anyways
    for i,xN in enumerate(x_NE2001[m]):
	try:
	    p = Popen([ ne2001_path + 'NE2001_dens','{0}'.format(xN),'{0}'.format(y_NE2001[i]),'{0}'.format(z_NE2001[i])], stdout = PIPE)
	except OSError:
	    warnings.warn('Could not run NE2001_dens, is it compiled correctly?',RuntimeWarning)
	    return -1
	out, err = p.communicate()
	#logging.debug('NE2001 code output: {0},{1}'.format(out,err))
	if not len(out):
	    warnings.warn('No output of NE2001_dens, is it compiled correctly?',RuntimeWarning)
	    return -1
# n is list with 
#      electron densities in cm^{-3}:
#	 ne1:    outer, thick disk
#	 ne2:    inner, thin disk (annular in form)
#	 nea:    spiral arms
#	 negc:   galactic center component
#	 nelism: local ISM component
#	 necN:   contribution from discrete 'clumps'
#	 nevN:   contribution from voids
	n = np.array(map(lambda x: float(x),out.split()))
	if n[-1] > 0.:	# we are in a void: supersedes all but clumps
	    ne[i] = n[-2] + n[-1]
	else:
	    ne[i] = np.sum(n)
    f = open(saveoutput,'w')
    pickle.dump(ne,f)
    f.close()

    return ne
