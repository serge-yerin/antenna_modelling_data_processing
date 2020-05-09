# Python3
#*************************************************************
#                        PARAMETERS                          *
#*************************************************************
max = 8000

#*************************************************************
#                   IMPORT LIBRARIES                         *
#*************************************************************
import pylab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#from matplotlib.pyplot import figure, show, grid, tight_layout

# Original pattern
theta = np.linspace(0, 2. * np.pi, 360)
r = 2 *  (1. + np.cos(theta+np.pi/6))**12

# 'Horizontal mirror'
theta_1 = np.zeros_like(theta)
for i in range(360):
    theta_1[i] = theta[359-i]

# 'Vertical mirror'
theta_2 = np.zeros_like(theta)
for i in range(180):
    theta_2[i] = theta[179-i]
    theta_2[180+i] = theta[359-i]

#'Cardinal point mirror'
theta_3 = np.zeros_like(theta)
theta_3a = np.zeros_like(theta)
for i in range(360):
    theta_3a[i] = theta[359-i]
for i in range(180):
    theta_3[i] = theta_3a[179-i]
    theta_3[180+i] = theta_3a[359-i]



radii = [max/4, max/2, 3*max/4, max]
man_pad = -4.5
man_angle = 90
rc('font', size=4, weight='bold')
fig = plt.figure(figsize=(9, 5))
ax1 = fig.add_subplot(221, projection='polar')
ax1.plot(theta, r, color = 'C0')
ax1.set_title('Original pattern', fontsize=6)
ax1.tick_params(pad=man_pad)
ax1.set_rgrids(radii, angle = man_angle)
ax2 = fig.add_subplot(222, projection='polar')
ax2.plot(theta_1, r, color = 'C1')
ax2.tick_params(pad=man_pad)
ax2.set_rgrids(radii, angle = man_angle)
ax2.set_title('Horizontal mirror', fontsize=6)
ax3 = fig.add_subplot(223, projection='polar')
ax3.plot(theta_2, r, color = 'C3')
ax3.tick_params(pad=man_pad)
ax3.set_rgrids(radii, angle = man_angle)
ax3.set_title('Vertical mirror', fontsize=6)
ax4 = fig.add_subplot(224, projection='polar')
ax4.plot(theta_3, r, color = 'C4')
ax4.tick_params(pad=man_pad)
ax4.set_rgrids(radii, angle = man_angle)
ax4.set_title('Cardinal point mirror', fontsize=6)
fig.subplots_adjust(wspace = -0.5, hspace=0.25, top=0.90) #
fig.suptitle('Mirror of antenna patterns check', fontsize=8, fontweight='bold')
pylab.savefig('Fig.1.png', bbox_inches='tight', dpi=300)
plt.close('all')

