import matplotlib.pyplot as plt
import numpy as np
from numpy import arange,sqrt,ma
from numpy.ma import masked_array
from pylab import figure,arange,colorbar,ma,sqrt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as clr
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.signal import savgol_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

SMALL_SIZE = 12
MEDIUM_SIZE = 16
BIGGER_SIZE = 18

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE) 

k_filter = 4
win_length = 53
# Make an output directory
dirname='./2FigComp/'
if not os.path.exists(dirname):
    os.makedirs(dirname)

trPR = [] # tracer for peak ring material < 25 GPa
for line in open('./C30/I30_V15/RingIndex.txt',"r"):
    trPR.append(int(line.rstrip('\n')))
tr25GPa = [] # tracer for all material < 25 GPa
for line in open('./C30/I30_V15/25Index.txt',"r"):
    tr25GPa.append(int(line.rstrip('\n')))
trPR_all = [] # tracer for the entire peak ring material
for line in open('./C30/I30_V15/RingIndexAll.txt',"r"):
    trPR_all.append(int(line.rstrip('\n')))
# tracer mass
mass = np.load('./C30/I30_V15/mass_C30_V15.npy')
# material x, y position
mat_x = np.load('./C30/I30_V15/modelx_C30_V15.npy')
mat_y = np.load('./C30/I30_V15/modely_C30_V15.npy')

trPR2 = [] 
for line in open('./C45/I30_V15/RingIndex.txt',"r"):
    trPR2.append(int(line.rstrip('\n')))
tr25GPa2 = [] 
for line in open('./C45/I30_V15/25Index.txt',"r"):
    tr25GPa2.append(int(line.rstrip('\n')))
trPR2_all = [] 
for line in open('./C45/I30_V15/RingIndexAll.txt',"r"):
    trPR2_all.append(int(line.rstrip('\n')))
mass2 = np.load('./C45/I30_V15/mass_C45_V15.npy')
mat2_x = np.load('./C45/I30_V15/modelx_C45_V15.npy')
mat2_y = np.load('./C45/I30_V15/modely_C45_V15.npy')

# Set up a figure
fig, axs = plt.subplots(1, 2, sharey=True, figsize=(12, 5))
# Remove horizontal space between axes
fig.subplots_adjust(wspace=0)
axs[0].spines['bottom'].set_linewidth(1.6)
axs[0].spines['left'].set_linewidth(1.6)
axs[0].spines['right'].set_linewidth(1.6)
axs[0].spines['top'].set_linewidth(1.6)
axs[1].spines['bottom'].set_linewidth(1.6)
axs[1].spines['left'].set_linewidth(1.6)
axs[1].spines['right'].set_linewidth(1.6)
axs[1].spines['top'].set_linewidth(1.6)
cmap2 = LinearSegmentedColormap.from_list('mycmap2', [ '#ffe5e5', '#e50000'])
cmap = LinearSegmentedColormap.from_list('mycmap', [ '#80755b', '#f0f0f0'])
xmajorLocator = MultipleLocator(50)
ymajorLocator = MultipleLocator(30)

time = [0, 5, 75, 110, 180]
for i in time:
    trp1 = np.load('./C30/I30_V15/TrP_C30_V15.npy')
    trp2 = np.load('./C45/I30_V15/TrP_C45_V15.npy')

    with np.load('./C30/I30_V15/mat_C30_V15-T{:03d}.npz'.format(i)) as npz:
        mat = np.ma.MaskedArray(**npz)
    with np.load('./C45/I30_V15/mat_C45_V15-T{:03d}.npz'.format(i)) as npz:
        mat2 = np.ma.MaskedArray(**npz)

    xmark = np.load('./C30/I30_V15/xmark_C30_V15-T{:03d}.npy'.format(i))
    ymark = np.load('./C30/I30_V15/ymark_C30_V15-T{:03d}.npy'.format(i))
    xmark2 = np.load('./C45/I30_V15/xmark_C45_V15-T{:03d}.npy'.format(i))
    ymark2 = np.load('./C45/I30_V15/ymark_C45_V15-T{:03d}.npy'.format(i))

    #plot materials
    axs[0].yaxis.set_major_locator(ymajorLocator)
    axs[0].set_xticks([-250, -200, -150, -100, -50])
    axs[0].set_yticks([-100, -50, 0, 50])
    axs[0].pcolormesh(-mat_x,mat_y,mat,cmap=cmap,vmin=1,vmax=3)
    axs[0].set_xlabel('Distance (km)')
    axs[0].set_xlim([-260,0])
    axs[0].set_ylim([-135, 75])
    axs[0].set_ylabel('Depth (km)')
    axs[1].xaxis.set_major_locator(xmajorLocator)
    axs[1].pcolormesh(mat2_x,mat2_y,mat2,cmap=cmap,vmin=1,vmax=3)   
    axs[1].set_xlabel('Distance (km)')
    axs[1].set_xlim(0,260)
    axs[1].text(.95, .95, "t ={: 2.0f} s".format(i*10), transform=axs[1].transAxes, ha="right", va="top")
    fig.savefig('{}/30v45Mat-{:03d}.jpg'.format(dirname,i),dpi=300, bbox_inches = 'tight', pad_inches = 0)
    axs[0].cla()
    axs[1].cla()

    rim30 = -478/2
    rim45 = 472/2
    if i > 75:
        axs[0].annotate('peak ring', xy=(-101, 3), xytext=(-101, 25),
            arrowprops=dict(facecolor='black', shrink=0.05, width = 3, headwidth = 5, headlength = 5),va = 'center', ha='center')
        axs[0].annotate('rim', xy=(rim30, 3), xytext=(rim30, 25),
            arrowprops=dict(facecolor='black', shrink=0.05, width = 3, headwidth = 5, headlength = 5),va = 'center', ha='center')
        axs[1].annotate('peak ring', xy=(95, 3), xytext=(95, 25),
            arrowprops=dict(facecolor='black', shrink=0.05, width = 3, headwidth = 5, headlength = 5),va = 'center', ha='center')
        axs[1].annotate('rim', xy=(rim45, 3), xytext=(rim45, 25),
            arrowprops=dict(facecolor='black', shrink=0.05, width = 3, headwidth = 5, headlength = 5),va = 'center', ha='center')
    
    if i == 180:
        axs[0].set_xlabel('Distance (km)')
        axs[1].set_xlabel('Distance (km)')
        axs[0].set_ylabel('Depth (km)')
        plt.rc('xtick', labelsize=8) 
        values = [1, 5, 20, 25, 30, 50, 100]
        norm = clr.BoundaryNorm(values, ncolors = 256)
        cmap_trp = matplotlib.cm.RdYlBu_r
        axs[0].yaxis.set_major_locator(ymajorLocator)
        axs[0].set_xticks([-250, -200, -150, -100, -50])
        axs[0].set_yticks([-100, -50, 0, 50])
        axs[0].set_xlim([-260,0])
        axs[0].set_ylim([-135, 75])
        axs[1].set_xlim(0,260)
        im1 = axs[0].scatter(-xmark, ymark, 
        c=trp1, norm = norm, cmap=cmap_trp, s=30, linewidths=0, alpha = 0.7)
        axs[1].scatter(xmark2, ymark2, 
        c=trp2, norm = norm, cmap=cmap_trp, s=30, linewidths=0, alpha = 0.7)
        axs[1].text(.95, .95, "t ={: 2.0f} s".format(i*10), transform=axs[1].transAxes, ha="right", va="top")
        axins = inset_axes(axs[0], width = '50%', height = '3%', loc = 'upper left', borderpad = 0.6)
        cb1 = plt.colorbar(im1, cax = axins, orientation = 'horizontal')
        cb1.set_label(r'$P_{peak}$ (GPa)', fontsize = 12)
        fig.savefig('{}/30v45MatTrP-{:03d}.jpg'.format(dirname,i),dpi=300, bbox_inches = 'tight', pad_inches = 0)
        axs[0].cla()
        axs[1].cla()
