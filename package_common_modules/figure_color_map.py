'''
'''
import matplotlib.pyplot as plt
from matplotlib import rc
import pylab

# Plotting color map figures
def figure_color_map(data, XLabel, YLabel, SupTitle, Title, FileName):
    fig = plt.figure(figsize = (10, 6))
    ax1 = fig.add_subplot(111)
    rc('font', weight='normal', size = 8)
    plot = ax1.pcolor(data, cmap = 'jet')
    ax1.set_aspect('equal', 'box') #axes().set_aspect('equal')
    ax1.set_xlabel(XLabel)
    ax1.set_ylabel(YLabel)
    fig.colorbar(plot)
    fig.suptitle(SupTitle, fontsize=8, fontweight='bold')
    ax1.set_title(Title, fontsize = 8)
    # plt.subplots_adjust(left=0.2, right=0.3)
    pylab.savefig(FileName, bbox_inches='tight', dpi = 160)
    plt.close('all')
