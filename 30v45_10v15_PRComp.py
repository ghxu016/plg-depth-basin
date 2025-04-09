import matplotlib.pyplot as plt
import numpy as np
from numpy import arange,sqrt,ma
from numpy.ma import masked_array
from pylab import figure,arange,colorbar,ma,sqrt
# Need this for the colorbars we will make on the mirrored plot
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 18

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)     # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE) 

# Make an output directory
dirname='./4FigComp/'
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

trPR3 = [] 
for line in open('./C30/I30_V10/RingIndex.txt',"r"):
    trPR3.append(int(line.rstrip('\n')))
tr25GPa3 = [] 
for line in open('./C30/I30_V10/25Index.txt',"r"):
    tr25GPa3.append(int(line.rstrip('\n')))
trPR3_all = [] 
for line in open('./C30/I30_V10/RingIndexAll.txt',"r"):
    trPR3_all.append(int(line.rstrip('\n')))
mass3 = np.load('./C30/I30_V10/mass_C30_V10.npy')
mat3_x = np.load('./C30/I30_V10/modelx_C30_V10.npy')
mat3_y = np.load('./C30/I30_V10/modely_C30_V10.npy')

trPR4 = [] 
for line in open('./C45/I30_V10/RingIndex.txt',"r"):
    trPR4.append(int(line.rstrip('\n')))
tr25GPa4 = [] 
for line in open('./C45/I30_V10/25Index.txt',"r"):
    tr25GPa4.append(int(line.rstrip('\n')))
trPR4_all = [] 
for line in open('./C45/I30_V10/RingIndexAll.txt',"r"):
    trPR4_all.append(int(line.rstrip('\n')))
mass4 = np.load('./C45/I30_V10/mass_C45_V10.npy')
mat4_x = np.load('./C45/I30_V10/modelx_C45_V10.npy')
mat4_y = np.load('./C45/I30_V10/modely_C45_V10.npy')

# Set up a figure
fig, axs = plt.subplots(2, 2, sharey=True, figsize=(10, 8))
# Remove horizontal space between axes
fig.subplots_adjust(wspace=0, hspace=0.03)
axs[1][0].spines['bottom'].set_linewidth(1.5)
axs[1][1].spines['bottom'].set_linewidth(1.5)
axs[0][0].spines['bottom'].set_linewidth(1.5)
axs[0][1].spines['bottom'].set_linewidth(1.5)
axs[0][0].spines['left'].set_linewidth(1.5)
axs[1][0].spines['left'].set_linewidth(1.5)
axs[0][1].spines['right'].set_linewidth(1.5)
axs[1][1].spines['right'].set_linewidth(1.5)
axs[0][0].spines['right'].set_linewidth(1.5)
axs[1][0].spines['right'].set_linewidth(1.5)
axs[0][0].spines['top'].set_linewidth(1.5)
axs[0][1].spines['top'].set_linewidth(1.5)
axs[1][0].spines['top'].set_linewidth(1.5)
axs[1][1].spines['top'].set_linewidth(1.5)
cmap2 = LinearSegmentedColormap.from_list('mycmap2', [ '#ffe5e5', '#e50000'])
cmap = LinearSegmentedColormap.from_list('mycmap', [ '#80755b', '#f0f0f0'])
xmajorLocator = MultipleLocator(50)
ymajorLocator = MultipleLocator(30)

time = [0, 5, 75, 110, 180]
for i in time:
    # load material type in each cells
    with np.load('./C30/I30_V15/mat_C30_V15-T{:03d}.npz'.format(i)) as npz:
        mat = np.ma.MaskedArray(**npz)
    with np.load('./C45/I30_V15/mat_C45_V15-T{:03d}.npz'.format(i)) as npz:
        mat2 = np.ma.MaskedArray(**npz)
    with np.load('./C30/I30_V10/mat_C30_V10-T{:03d}.npz'.format(i)) as npz:
        mat3 = np.ma.MaskedArray(**npz)
    with np.load('./C45/I30_V10/mat_C45_V10-T{:03d}.npz'.format(i)) as npz:
        mat4 = np.ma.MaskedArray(**npz)

    #load tracer position
    xmark = np.load('./C30/I30_V15/xmark_C30_V15-T{:03d}.npy'.format(i))
    ymark = np.load('./C30/I30_V15/ymark_C30_V15-T{:03d}.npy'.format(i))
    xmark2 = np.load('./C45/I30_V15/xmark_C45_V15-T{:03d}.npy'.format(i))
    ymark2 = np.load('./C45/I30_V15/ymark_C45_V15-T{:03d}.npy'.format(i))
    xmark3 = np.load('./C30/I30_V10/xmark_C30_V10-T{:03d}.npy'.format(i))
    ymark3 = np.load('./C30/I30_V10/ymark_C30_V10-T{:03d}.npy'.format(i))
    xmark4 = np.load('./C45/I30_V10/xmark_C45_V10-T{:03d}.npy'.format(i))
    ymark4 = np.load('./C45/I30_V10/ymark_C45_V10-T{:03d}.npy'.format(i))

    #plot materials
    axs[0][0].yaxis.set_major_locator(ymajorLocator)
    axs[0][0].set_xticks([-150, -100, -50])
    axs[0][0].pcolormesh(-mat_x,mat_y,mat,cmap=cmap,vmin=1,vmax=3)
    axs[0][0].set_xlim([-180,0])
    axs[0][0].set_ylim([-105,45])
    axs[0][1].xaxis.set_major_locator(xmajorLocator)
    axs[0][1].pcolormesh(mat2_x,mat2_y,mat2,cmap=cmap,vmin=1,vmax=3)
    axs[0][1].set_xlim(0,180)
    axs[1][0].yaxis.set_major_locator(ymajorLocator)
    axs[1][0].set_xticks([-150, -100, -50])
    axs[1][0].pcolormesh(-mat3_x,mat3_y,mat3,cmap=cmap,vmin=1,vmax=3)
    axs[1][0].set_xlim([-180,0])
    axs[1][0].set_ylim([-105,45])
    axs[1][1].xaxis.set_major_locator(xmajorLocator)
    axs[1][1].pcolormesh(mat4_x,mat4_y,mat4,cmap=cmap,vmin=1,vmax=3)
    axs[1][1].set_xlim(0,180)
    
    fig.text(0.5, 0.04, 'Distance to Center (km)', ha='center', fontsize = 16)
    fig.text(0.05, 0.5, 'Depth (km)', va='center', rotation='vertical', fontsize = 16)

    # plot the location of the entire peak ring material 
    axs[0][0].scatter(-xmark[trPR_all], ymark[trPR_all], c=mass[trPR_all]/7e14, vmin = 0, vmax = 1, cmap=cmap2, s=30, linewidths=0)
    im1 = axs[0][1].scatter(xmark2[trPR2_all], ymark2[trPR2_all], c=mass2[trPR2_all]/7e14, vmin = 0, vmax = 1, cmap=cmap2, s=30, linewidths=0)
    im1 = axs[1][0].scatter(-xmark3[trPR3_all], ymark3[trPR3_all], c=mass3[trPR3_all]/5e14, vmin = 0, vmax = 1, cmap=cmap2, s=30, linewidths=0)
    axs[1][1].scatter(xmark4[trPR4_all], ymark4[trPR4_all], c=mass4[trPR4_all]/5e14, vmin = 0, vmax = 1, cmap=cmap2, s=30, linewidths=0)
    
    axs[0][0].text(.05, .95, "t ={:2.0f} s".format(i*10), transform=axs[0][0].transAxes, ha="left", va="top")
    axs[0][0].text(.55, .85, "$D_{proj}$ = 30 km", transform=axs[0][0].transAxes, ha="left", va="top")
    # axs[0][1].text(.95, .95, "t={:2.0f} s".format(time), transform=axs[1].transAxes, ha="right", va="top")
    axs[0][0].text(.05, .6, "$T_{c}$ = 30 km", transform=axs[0][0].transAxes, ha="left", va="top")
    axs[0][1].text(.95, .6, "45 km", transform=axs[0][1].transAxes, ha="right", va="top")
    axs[1][0].text(.05, .6, "$T_{c}$ = 30 km", transform=axs[1][0].transAxes, ha="left", va="top")
    axs[1][1].text(.95, .6, "45 km", transform=axs[1][1].transAxes, ha="right", va="top")

    axs[0][0].text(.05, .45, "$v$ = 15 km/s", transform=axs[0][0].transAxes, ha="left", va="top")
    axs[1][0].text(.05, .45, "$v$ = 10 km/s", transform=axs[1][0].transAxes, ha="left", va="top")

    if i == 0:
        # plot the shock zone of > 25 GPa
        # which will be ignored in calculation 
        x30 = -xmark[tr25GPa]
        y30 = ymark[tr25GPa]
        x30 = x30.astype(int)
        y30 = y30.astype(int)
        x45 = xmark2[tr25GPa2]
        y45 = ymark2[tr25GPa2]
        x45 = x45.astype(int)
        y45 = y45.astype(int)
        axs[0][0].scatter(x30, y30, c='k', s = 0.8, linewidths=0)
        axs[0][1].scatter(x45, y45, c='k', s = 0.8, linewidths=0)

        xx30 = -xmark3[tr25GPa3]
        yy30 = ymark3[tr25GPa3]
        xx30 = xx30.astype(int)
        yy30 = yy30.astype(int)
        xx45 = xmark4[tr25GPa4]
        yy45 = ymark4[tr25GPa4]
        xx45 = xx45.astype(int)
        yy45 = yy45.astype(int)
        axs[1][0].scatter(xx30, yy30, c='k', s = 0.8, linewidths=0)
        axs[1][1].scatter(xx45, yy45, c='k', s = 0.8, linewidths=0)
        
    
    axins = inset_axes(axs[0][1], width = '50%', height = '8%', loc = 'upper right', borderpad = 0.8)
    cb1 = plt.colorbar(im1, cax = axins, orientation = 'horizontal')
    cb1.ax.tick_params(labelsize=10)
    cb1.set_ticks(ticks=[0, 0.5, 1])
    cb1.set_label('Normalized tracer mass', fontsize = 10)
    fig.savefig('{}/30v45_10v15Ring-{:03d}.jpg'.format(dirname,i),dpi=300, bbox_inches = 'tight', pad_inches = 0)
    axs[0][0].cla()
    axs[0][1].cla()
    axs[1][0].cla()
    axs[1][1].cla()