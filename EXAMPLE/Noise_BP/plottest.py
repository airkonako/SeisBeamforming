import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#-- Generate Data -----------------------------------------
# Using linspace so that the endpoint of 360 is included...
azimuths = np.radians(np.linspace(0, 360, 20))
zeniths = np.arange(0, 70, 10)

r, theta = np.meshgrid(zeniths, azimuths)
values = np.random.random((azimuths.size, zeniths.size))

colormap = cm.rainbow
#-- Plot... ------------------------------------------------
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
ax.set_xticklabels(['E', '', 'N', '', 'W', '', 'S', ''])
ax.set_rticks([10, 40]) 
cs = ax.contourf(theta, r, values, cmap=colormap)

cbar = fig.colorbar(cs, orientation="horizontal", fraction=0.05, pad=0.1)
cbar.set_label('Beam power [dB]')


plt.show()