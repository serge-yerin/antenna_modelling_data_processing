'''
'''
import matplotlib.pyplot as plt
from matplotlib import rc
import pylab

# Plotting color map figures
def figure_color_map(data, XLabel, YLabel, SupTitle, Title, FileName):
    plt.figure(figsize=(10, 6))
    rc('font', weight='normal', size = 8)
    plt.pcolor(data, cmap = 'jet')
    plt.axes().set_aspect('equal')
    plt.xlabel(XLabel)
    plt.ylabel(YLabel)
    plt.colorbar()
    plt.suptitle(SupTitle, fontsize=8, fontweight='bold')
    plt.title(Title, fontsize = 8)
    # plt.subplots_adjust(left=0.2, right=0.3)
    pylab.savefig(FileName, bbox_inches='tight', dpi = 160)
    plt.close('all')
