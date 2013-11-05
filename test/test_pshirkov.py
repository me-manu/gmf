import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize


class Test_Pshirkov(object):
    """
    class to test implementation of Pshirkov GMF model
    """

    def test_disk(self):
	"""
	test the b field model by pshirkov
	"""
	import gmf.gmf as GMF

	gmfm = GMF.GMF_Pshirkov(mode = 'BSS')
	#gmfm = GMF.GMF_Pshirkov(mode = 'ASS')
	#gmfm = GMF.GMF()

	x = np.linspace(-20.,20.,500)
	y = np.linspace(-20.,20.,x.shape[0])
	z = np.ones(x.shape[0]) * 1.1

	xx,yy = np.meshgrid(x,y)

        rr  = np.sqrt(xx**2. + yy**2.)
        pp  = np.arctan2(yy,xx)

	B  = np.zeros((x.shape[0],y.shape[0]))

	idx0 = x.shape[0] / 2	# index where x is closest to zero
	u,v  = np.zeros(x.shape[0]/3), np.zeros(x.shape[0]/3)

	for i,r in enumerate(rr[:]):
	    Bd,Babs_d	= gmfm.Bdisk(r,pp[i],z)
	    Bh,Babs_h	= gmfm.Bhalo(r,z)
	    Btot	= Bd + Bh
	    B[i] = np.sqrt(np.sum((Btot)**2.,axis = 0)) * GMF.signum(Btot[1])


	    if i == idx0:	# calculate length of arrows for quiver function, only do it for a third of the points 
				# and at x close to zero
		u = Btot[0,::3] * np.cos(pp[i,::3]) - Btot[1,::3] * np.sin(pp[i,::3])	# this comes from coordinate transformation
		v = Btot[0,::3] * np.sin(pp[i,::3]) + Btot[1,::3] * np.cos(pp[i,::3])

	### --- Plot the field
	###

	fig	= plt.figure(figsize=(12,9))
	ax	= fig.add_subplot(1,1,1)
	norm	= Normalize(vmin = -3., vmax = 3.)
	im	= plt.imshow(B, interpolation='bilinear', origin='lower',
	cmap	= cm.RdBu, aspect=1.,norm=norm,extent=(x[0],x[-1],y[0],y[-1]))
	#cmap	=cm.RdBu, aspect=1.,norm=norm,extent=(x[0],x[-1],z[0],z[-1]))

	plt.quiver(x[::3],y[idx0],u,v,width=0.002, headwidth = 3)	# draw the vector field

	CBI	= plt.colorbar(im, orientation='vertical', shrink=0.9,format = "%.1f")

	CBI.set_label(r"$B\,(\mu\mathrm{G})$")
	plt.xlabel('$x$ (kpc)')
	plt.ylabel('$y$ (kpc)')

	#plt.ylabel('$z$ (kpc)')
	plt.axis([np.min(x), np.max(x), np.min(y), np.max(y)])
	#plt.axis([np.min(x), np.max(x), np.min(z), np.max(z)])
	plt.show()

	return
