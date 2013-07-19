"""
Class for calculating the galactic magnetic field according to Jansson & Farrar (2012)
see  http://adsabs.harvard.edu/abs/2012arXiv1204.3662J

History:
- 05/07/2012 created
"""
__version__ = '0.1'
__author__ ="M. Meyer // manuel.meyer@physik.uni-hamburg.de"

import numpy as np
import logging
import warnings

signum = lambda x: (x < 0.) * -1. + (x >= 0) * 1.
pi = np.pi

class GMF(object):
    """
    Class with analytical functions that describe the 
    galactic magnetic field according to the model of Jannson & Farrar (2012)

    Only the regular field components are implemented. The striated field component is missing as well.

    Attributes
    ----------
    Disk:
	bring, bring_unc	: floats, field strength in ring at 3 kpc < rho < 5 kpc
	hdisk, hdisk_unc	: float, disk/halo transition height
	wdisk, wdisk_unc	: floats, transition width
	b, b_unc		: (8,1)-dim np.arrays, field strengths of spiral arms at 5 kpc
	rx			: (8,1)-dim np.array, dividing lines of spiral arms, coordinates of neg. x-axes that intersect with arm 
	idisk			: float, spiral arms opening angle
    Halo:
	Bn, Bn_unc		: floats, field strength northern halo
	Bs, Bs_unc		: floats, field strength southern halo
	rhon, rhon_unc		: floats, transition radius north
	rhos, rhos_unc		: floats, transition radius south, lower limit
	whalo, whalo_unc	: floats, transition width
	z0, z0_unc		: floats, vertical scale height
    Out of plaxe or "X" component:
	BX0, BX_unc		: floats, field strength at origin
	ThetaX0, ThetaX0_unc	: floats, elev. angle at z = 0, rho > rhoXc
	rhoXc, rhoXc_unc	: floats, radius where thetaX = thetaX0
	rhoX, rhoX_unca		: floats, exponential scale length

    striated field:
	gamma, self.gamma_unc	: floats striation and / or rel. elec. number dens. rescaling

    Notes
    -----
    Paper:
    http://adsabs.harvard.edu/abs/2012arXiv1204.3662J
    Jansson & Farrar (2012)
    """

    def __init__(self):
	"""
	Init the GMF class,
	all B-field values are in muG

	Parameters
	----------
	None

	Returns
	-------
	Nothing

	"""
	# Best fit values, see Table 1 of Jansson & Farrar --------#
	# Disk
	self.bring, self.bring_unc	= 0.1,0.1	# ring at 3 kpc < rho < 5 kpc
	self.hdisk, self.hdisk_unc	= 0.4, 0.03	# disk/halo transition
	self.wdisk, self.wdisk_unc	= 0.27,0.08	# transition width
	self.b		= np.array([0.1,3.,-0.9,-0.8,-2.0,-4.2,0.,2.7])	# field strength of spiral arms at 5 kpc
	self.b_unc	= np.array([1.8,0.6,0.8,0.3,0.1,0.5,1.8,1.8])
	self.rx		= np.array([5.1,6.3,7.1,8.3,9.8,11.4,12.7,15.5])# dividing lines of spiral lines
	self.idisk	= 11.5 * pi/180.		# spiral arms opening angle
	# Halo
	self.Bn, self.Bn_unc		= 1.4,0.1	# northern halo
	self.Bs, self.Bs_unc		= -1.1,0.1	# southern halo
	self.rhon, self.rhon_unc	= 9.22,0.08	# transition radius north
	self.rhos, self.rhos_unc	= 16.7,0.	# transition radius south, lower limit
	self.whalo, self.whalo_unc	= 0.2,0.12	# transition width
	self.z0, self.z0_unc		= 5.3, 1.6	# vertical scale height
	# Out of plaxe or "X" component
	self.BX0, self.BX_unc		= 4.6,0.3	# field strength at origin
	self.ThetaX0, self.ThetaX0_unc	= 49. * pi/180., pi/180. # elev. angle at z = 0, rho > rhoXc
	self.rhoXc, self.rhoXc_unc	= 4.8, 0.2	# radius where thetaX = thetaX0
	self.rhoX, self.rhoX_unc	= 2.9, 0.1	# exponential scale length
	# striated field
	self.gamma, self.gamma_unc	= 2.92,0.14	# striation and / or rel. elec. number dens. rescaling
	return

    def L(self,z,h,w):
	"""
	Transition function, see Ronnie & Farrar Eq. 5

	Parameters:
	-----------
	z: scalar or numpy array with positions (height above disk, z; distance from center, rho)
	h: scaler: height parameter
	w: scalar: width parameter

	Returns:
	--------
	numpy array with transition function values
	"""
	if np.isscalar(z):
	    z = np.array([z])
	ones = np.ones(z.shape[0])
	return 1./(ones + np.exp(-2. * (np.abs(z) - h) / w))

    def r_log_spiral(self,phi):
	"""
	return distance from center for angle phi of logarithmic spiral

	Parameters
	----------
	phi: scalar or np.array with polar angle values

	Returns
	-------
	r(phi) = rx * exp(b * phi) as np.array

	Notes
	-----
	see http://en.wikipedia.org/wiki/Logarithmic_spiral
	"""
	if np.isscalar(phi):
	    phi = np.array([phi])
	ones = np.ones(phi.shape[0])

	# self.rx.shape = 8
	# phi.shape = p
	# then result is given as (8,p)-dim array, each row stands for one rx

	result = np.tensordot(self.rx , np.exp((phi - 3.*pi*ones) / np.tan(pi/2. - self.idisk)),axes = 0)
	result = np.vstack((result, np.tensordot(self.rx , np.exp((phi - pi*ones) / np.tan(pi/2. - self.idisk)),axes = 0) ))
	result = np.vstack((result, np.tensordot(self.rx , np.exp((phi + pi*ones) / np.tan(pi/2. - self.idisk)),axes = 0) ))
	return np.vstack((result, np.tensordot(self.rx , np.exp((phi + 3.*pi*ones) / np.tan(pi/2. - self.idisk)),axes = 0) ))

    def Bdisk(self,rho,phi,z):
	"""
	Disk component of galactic magnetic field 
	in galactocentric cylindrical coordinates (rho,phi,z)

	Parameters
	----------
	rho:	N-dim np.array,	distance from origin in GC cylindrical coordinates, is in kpc
	z:	N-dim np.array, height in kpc in GC cylindrical coordinates
	phi:	N-dim np.array, polar angle in GC cylindircal coordinates, in radian

	Returns
	-------
	Bdisk:	(3,N)-dim np.array with (rho,phi,z) components of disk field for each coordinate tuple
	|Bdisk|: N-dim np.array, absolute value of Bdisk for each coordinate tuple
	"""
	if (not rho.shape[0] == phi.shape[0]) and (not z.shape[0] == phi.shape[0]):
	    warnings.warn("List do not have equal shape! returning -1", RuntimeWarning)
	    return -1

	Bdisk = np.zeros((3,rho.shape[0]))	# Bdisk vector in rho, phi, z
						# rows: rho, phi and z component

	ones		= np.ones(rho.shape[0])
	m_center	= (rho >= 3.) & (rho < 5.1)
	m_disk		= (rho >= 5.1) & (rho <= 20.)

	Bdisk[1,m_center] = self.bring

	# Determine in which arm we are
	# this is done for each coordinate individually, possible to convert into array task?
	if np.sum(m_disk):
	    rls = self.r_log_spiral(phi[m_disk])

	    rls = np.abs(rls - rho[m_disk])
	    narm = np.argmin(rls, axis = 0) % 8

	    Bdisk[0,m_disk] = np.sin(self.idisk)* self.b[narm] * (5. / rho[m_disk])
	    Bdisk[1,m_disk] = np.cos(self.idisk)* self.b[narm] * (5. / rho[m_disk])

	Bdisk  *= (ones - self.L(z,self.hdisk,self.wdisk))

	return Bdisk, np.sqrt(np.sum(Bdisk**2.,axis = 0))

    def Bhalo(self,rho,z):
	"""
	Halo component of galactic magnetic field 
	in galactocentric cylindrical coordinates (rho,phi,z)

	Bhalo is purely azimuthal (toroidal), i.e. has only a phi component

	Parameters
	----------
	rho:	N-dim np.array,	distance from origin in GC cylindrical coordinates, is in kpc
	z:	N-dim np.array, height in kpc in GC cylindrical coordinates

	Returns
	-------
	Bhalo:	(3,N)-dim np.array with (rho,phi,z) components of halo field for each coordinate tuple
	|Bhalo|: N-dim np.array, absolute value of Bdisk for each coordinate tuple
	"""

	if (not rho.shape[0] == z.shape[0]):
	    warnings.warn("List do not have equal shape! returning -1", RuntimeWarning)
	    return -1

	Bhalo = np.zeros((3,rho.shape[0]))	# Bhalo vector in rho, phi, z
						# rows: rho, phi and z component

	ones		= np.ones(rho.shape[0])
	m = ( z != 0. )

	Bhalo[1,m] = np.exp(-np.abs(z[m])/self.z0) * self.L(z[m], self.hdisk, self.wdisk) * \
			( self.Bn * (ones[m] - self.L(rho[m], self.rhon, self.whalo)) * (z[m] > 0.) \
			+ self.Bs * (ones[m] - self.L(rho[m], self.rhos, self.whalo)) * (z[m] < 0.) )
	return Bhalo, np.sqrt(np.sum(Bhalo**2.,axis = 0))

    def BX(self,rho,z):
	"""
	X (out of plane) component of galactic magnetic field 
	in galactocentric cylindrical coordinates (rho,phi,z)

	BX is purely poloidal, i.e. phi component = 0

	Parameters
	----------
	rho:	N-dim np.array,	distance from origin in GC cylindrical coordinates, is in kpc
	z:	N-dim np.array, height in kpc in GC cylindrical coordinates

	Returns
	-------
	BX:	(3,N)-dim np.array with (rho,phi,z) components of halo field for each coordinate tuple
	|BX|: N-dim np.array, absolute value of Bdisk for each coordinate tuple
	"""

	if (not rho.shape[0] == z.shape[0]):
	    warnings.warn("List do not have equal shape! returning -1", RuntimeWarning)
	    return -1

	BX= np.zeros((3,rho.shape[0]))	# BX vector in rho, phi, z
					# rows: rho, phi and z component
	m = np.sqrt(rho**2. + z**2.) >= 1.

	bx = lambda rho_p: self.BX0 * np.exp(-rho_p / self.rhoX)
	tx = lambda rho,z,rho_p: np.arctan(np.abs(z)/(rho - rho_p))

	rho_p	= rho[m] *self.rhoXc/(self.rhoXc + np.abs(z[m] ) / np.tan(self.ThetaX0))

	m_rho_b = rho_p > self.rhoXc	# region with constant elevation angle
	m_rho_l = rho_p <= self.rhoXc	# region with varying elevation angle

	theta	= np.zeros(z[m].shape[0])
	b	= np.zeros(z[m].shape[0])

	rho_p0	= (rho[m])[m_rho_b]  - np.abs( (z[m])[m_rho_b] ) / np.tan(self.ThetaX0)
	b[m_rho_b]	= bx(rho_p0) * rho_p0/ (rho[m])[m_rho_b]
	theta[m_rho_b]	= self.ThetaX0 * np.ones(m_rho_b.shape[0])

	b[m_rho_l]	= bx(rho_p[m_rho_l]) * (rho_p[m_rho_l]/(rho[m])[m_rho_l] )**2.
	theta[m_rho_l]	= tx((rho[m])[m_rho_l] ,(z[m])[m_rho_l] ,rho_p[m_rho_l])
	mz = (z[m] == 0.)
	theta[mz]	= np.pi/2.
	#logging.debug('rho,z,rho_p, theta: {0:.3f}  {1:.3f}  {2:.3f}  {3:.3f}'.format(rho,z,rho_p, theta))

	BX[0,m] = b * (np.cos(theta) * (z[m] >= 0) + np.cos(pi*np.ones(theta.shape[0]) - theta) * (z[m] < 0))
	BX[2,m] = b * (np.sin(theta) * (z[m] >= 0) + np.sin(pi*np.ones(theta.shape[0]) - theta) * (z[m] < 0))

	return BX, np.sqrt(np.sum(BX**2.,axis=0))
