import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange,sqrt,ma
from numpy.ma import masked_array
from pylab import figure,arange,colorbar,ma,sqrt
# Need this for the colorbars we will make on the mirrored plot
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap

from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.signal import savgol_filter
# This example plotting script designed to plot 
# material and temperature in the Chicxulub example

# If viridis colormap is available, use it here
try:
    plt.set_cmap('viridis')
except:
    plt.set_cmap('YlGnBu_r')

# distances between tracers
def get_distances(s,line):
    x=s.xmark[line]
    y=s.ymark[line]
    return sqrt((x[:-1]-x[1:])**2+(y[:-1]-y[1:])**2)

# Define the maximum separation allowed when plotting lines
maxsep=3.

def make_colorbar(ax,p,f):
    # Create axes either side of the plot to place the colorbars
    divider = make_axes_locatable(ax)
    cx=divider.append_axes("left", size="5%", pad=0.7)
    cb=fig.colorbar(p,cax=cx)
    cb.set_label(psp.longFieldName(f))
    # Need to set the labels on the left for this colorbar
    cx.yaxis.tick_left()
    cx.yaxis.set_label_position('left')


SMALL_SIZE = 12
MEDIUM_SIZE = 14
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
dirname='I30/30v45_10v15/0117'
# dirname='I30/30v45_10v15/shade'
psp.mkdir_p(dirname)

# Open the datafile
tempdir = '/home/xgh/isale_runs_d/output/peak_ring/'
model=psp.opendatfile(tempdir+'30/v15/I30_V15/jdata.dat')
model2=psp.opendatfile(tempdir+'45/v15/I30_V15/jdata.dat')
model3=psp.opendatfile(tempdir+'30/v10/I30_V10/jdata.dat')
model4=psp.opendatfile(tempdir+'45/v10/I30_V10/jdata.dat')
[model.modelInfo()]
[model2.modelInfo()]
# Set the distance units to km
model.setScale('km')
model2.setScale('km')
model3.setScale('km')
model4.setScale('km')
# Set up a pylab figure
# fig, axs = plt.subplots(1, 2, sharey=True, figsize=(10,4.5))
fig, axs = plt.subplots(2, 2, sharey=True, figsize=(10, 8))
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

# ax=fig.add_subplot(111)
# Remove horizontal space between axes
fig.subplots_adjust(wspace=0, hspace=0.03)
cmap = LinearSegmentedColormap.from_list('mycmap', [ 'gray', 'sandybrown'])
cmap2 = LinearSegmentedColormap.from_list('mycmap2', [ '#ffe5e5', '#e50000'])
cmap = LinearSegmentedColormap.from_list('mycmap', [ '#80755b', '#f0f0f0'])
lb = '#aaaaaa'
trPR = [] 
for line in open(tempdir+'30/v15/I30_V15/RingIndex.txt',"r"):
    trPR.append(int(line.rstrip('\n')))

trPR2 = [] # Front facing PR tracers
for line in open(tempdir+'45/v15/I30_V15/RingIndex.txt',"r"):
    trPR2.append(int(line.rstrip('\n')))

trPR3 = [] # Front facing PR tracers
for line in open(tempdir+'30/v10/I30_V10/RingIndex.txt',"r"):
    trPR3.append(int(line.rstrip('\n')))

trPR4 = [] # Front facing PR tracers
for line in open(tempdir+'45/v10/I30_V10/RingIndex.txt',"r"):
    trPR4.append(int(line.rstrip('\n')))

tr25PR = []
for line in open(tempdir+'30/v15/I30_V15/25gpaRingIndex.txt',"r"):
    tr25PR.append(int(line.rstrip('\n')))

tr25PR2 = [] # Front facing PR tracers
for line in open(tempdir+'45/v15/I30_V15/25gpaRingIndex.txt',"r"):
    tr25PR2.append(int(line.rstrip('\n')))

tr25PR3 = [] # Front facing PR tracers
for line in open(tempdir+'30/v10/I30_V10/25gpaRingIndex.txt',"r"):
    tr25PR3.append(int(line.rstrip('\n')))

tr25PR4 = [] # Front facing PR tracers
for line in open(tempdir+'45/v10/I30_V10/25gpaRingIndex.txt',"r"):
    tr25PR4.append(int(line.rstrip('\n')))

trPR_full = [] 
for line in open(tempdir+'30/v15/I30_V15/RingIndexFULL.txt',"r"):
    trPR_full.append(int(line.rstrip('\n')))

trPR2_full = [] # Front facing PR tracers
for line in open(tempdir+'45/v15/I30_V15/RingIndexFULL.txt',"r"):
    trPR2_full.append(int(line.rstrip('\n')))

trPR3_full = [] # Front facing PR tracers
for line in open(tempdir+'30/v10/I30_V10/RingIndexFULL.txt',"r"):
    trPR3_full.append(int(line.rstrip('\n')))

trPR4_full = [] # Front facing PR tracers
for line in open(tempdir+'45/v10/I30_V10/RingIndexFULL.txt',"r"):
    trPR4_full.append(int(line.rstrip('\n')))


model.tracerMassVol()
mass = model.tracerMass
model2.tracerMassVol()
mass2 = model2.tracerMass
model3.tracerMassVol()
mass3 = model3.tracerMass
model4.tracerMassVol()
mass4 = model4.tracerMass

xmajorLocator = MultipleLocator(50)
ymajorLocator = MultipleLocator(30)
for i in arange(0,200,5):
    step=model.readStep(['Den'],i)
    step2=model2.readStep(['Den'],i)
    step3=model3.readStep(['Den'],i)
    step4=model4.readStep(['Den'],i)
    # Plot each graph, and manually set the y tick values
    #plot materials
    axs[0][0].yaxis.set_major_locator(ymajorLocator)
    axs[0][0].set_xticks([-150, -100, -50])
    # axs[0][0].set_xticks([-100, -50])
    axs[0][0].pcolormesh(-model.x,model.y,step.mat,cmap=cmap,vmin=1,vmax=model.nmat+1)
    # axs[0][0].set_xlabel('Distance to Center (km)')
    axs[0][0].set_xlim([-180,0])
    # axs[0][0].set_ylim([-110,60])
    # axs[0][0].set_xlim([-110,0])
    axs[0][0].set_ylim([-105,45])
    # axs[0][0].set_ylabel('Depth (km)')
    axs[0][1].xaxis.set_major_locator(xmajorLocator)
    axs[0][1].pcolormesh(model2.x,model2.y,step2.mat,
            cmap=cmap,vmin=1,vmax=model2.nmat+1)
    
    # axs[0][1].set_xlabel('Distance to Center (km)')
    axs[0][1].set_xlim(0,180)
    # axs[0][1].set_xlim(0,110)

    axs[1][0].yaxis.set_major_locator(ymajorLocator)
    axs[1][0].set_xticks([-150, -100, -50])
    # axs[1][0].set_xticks([-100, -50])
    axs[1][0].pcolormesh(-model3.x,model3.y,step3.mat,cmap=cmap,vmin=1,vmax=model3.nmat+1)
    # axs[1][0].set_xlabel('Distance to Center (km)')
    axs[1][0].set_xlim([-180,0])
    # axs[1][0].set_ylim([-110,60])
    # axs[1][0].set_xlim([-110,0])
    axs[1][0].set_ylim([-105,45])
    # axs[1][0].set_ylabel('Depth (km)')
    axs[1][1].xaxis.set_major_locator(xmajorLocator)
    axs[1][1].pcolormesh(model4.x,model4.y,step4.mat,
            cmap=cmap,vmin=1,vmax=model4.nmat+1)
    
    # axs[1][1].set_xlabel('Distance to Center (km)')
    fig.text(0.5, 0.04, 'Distance to Center (km)', ha='center', fontsize = 16)
    fig.text(0.05, 0.5, 'Depth (km)', va='center', rotation='vertical', fontsize = 16)
    axs[1][1].set_xlim(0,180)
    # axs[1][1].set_xlim(0,110)
    # plot peak rings

    axs[0][0].scatter(-step.xmark[trPR_full],
                    step.ymark[trPR_full], c=mass[trPR_full],
                    cmap=cmap2, s=30, linewidths=0)
    axs[0][1].scatter(step2.xmark[trPR2_full],
                    step2.ymark[trPR2_full], c=mass2[trPR2_full],
                    cmap=cmap2, s=30, linewidths=0)
    axs[1][0].scatter(-step3.xmark[trPR3_full],
                    step3.ymark[trPR3_full], c=mass3[trPR3_full],
                    cmap=cmap2, s=30, linewidths=0)
    axs[1][1].scatter(step4.xmark[trPR4_full],
                    step4.ymark[trPR4_full], c=mass4[trPR4_full],
                    cmap=cmap2, s=30, linewidths=0)
    # Material boundaries
    [axs[0][0].contour(-model.xc,model.yc,step.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
    [axs[0][1].contour(model2.xc,model2.yc,step2.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
    [axs[1][0].contour(-model3.xc,model3.yc,step3.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
    [axs[1][1].contour(model4.xc,model4.yc,step4.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
    
    axs[0][0].text(.05, .95, "t ={:2.0f} s".format(step.time), transform=axs[0][0].transAxes, ha="left", va="top")
    axs[0][0].text(.55, .85, "$D_{proj}$ = 30 km", transform=axs[0][0].transAxes, ha="left", va="top")
    # axs[0][1].text(.95, .95, "t={:2.0f} s".format(step.time), transform=axs[1].transAxes, ha="right", va="top")
    axs[0][0].text(.05, .6, "$T_{c}$ = 30 km", transform=axs[0][0].transAxes, ha="left", va="top")
    axs[0][1].text(.95, .6, "45 km", transform=axs[0][1].transAxes, ha="right", va="top")
    axs[1][0].text(.05, .6, "$T_{c}$ = 30 km", transform=axs[1][0].transAxes, ha="left", va="top")
    axs[1][1].text(.95, .6, "45 km", transform=axs[1][1].transAxes, ha="right", va="top")

    axs[0][0].text(.05, .45, "$v$ = 15 km/s", transform=axs[0][0].transAxes, ha="left", va="top")
    axs[1][0].text(.05, .45, "$v$ = 10 km/s", transform=axs[1][0].transAxes, ha="left", va="top")

    if i == 0:
        x30 = -step.xmark[tr25PR]
        y30 = step.ymark[tr25PR]
        x30 = x30.astype(int)
        y30 = y30.astype(int)
        x45 = step2.xmark[tr25PR2]
        y45 = step2.ymark[tr25PR2]
        x45 = x45.astype(int)
        y45 = y45.astype(int)
        # axs[0][0].scatter(-step.xmark[trPR_full],
        #             step.ymark[trPR_full], c=lb,
        #             s=15, linewidths=0)
        # axs[0][1].scatter(step2.xmark[trPR2_full],
        #             step2.ymark[trPR2_full], c=lb,
        #              s=15, linewidths=0)
        # axs[1][0].scatter(-step3.xmark[trPR3_full],
        #             step3.ymark[trPR3_full], c=lb,
        #             s=15, linewidths=0)
        # axs[1][1].scatter(step4.xmark[trPR4_full],
        #             step4.ymark[trPR4_full], c=lb,
        #             s=15, linewidths=0)
        axs[0][0].scatter(x30, y30, c='k', s = 0.8, linewidths=0)
        axs[0][1].scatter(x45, y45, c='k', s = 0.8, linewidths=0)

        xx30 = -step.xmark[tr25PR3]
        yy30 = step.ymark[tr25PR3]
        xx30 = xx30.astype(int)
        yy30 = yy30.astype(int)
        xx45 = step2.xmark[tr25PR4]
        yy45 = step2.ymark[tr25PR4]
        xx45 = xx45.astype(int)
        yy45 = yy45.astype(int)
        axs[1][0].scatter(xx30, yy30, c='k', s = 0.8, linewidths=0)
        axs[1][1].scatter(xx45, yy45, c='k', s = 0.8, linewidths=0)
        
    
    
    fig.savefig('{}/30v45_10v15Ring-{:05d}.jpg'.format(dirname,i),dpi=300, bbox_inches = 'tight', pad = 0)
    axs[0][0].cla()
    axs[0][1].cla()
    axs[1][0].cla()
    axs[1][1].cla()

#Draw Choosing Time
step=model.readStep(['Den'],65)
step2=model2.readStep(['Den'],75)
step3=model3.readStep(['Den'],62)
step4=model4.readStep(['Den'],70)
# Plot each graph, and manually set the y tick values
#plot materials
axs[0][0].yaxis.set_major_locator(ymajorLocator)
axs[0][0].set_xticks([-150, -100, -50])
# axs[0][0].set_xticks([-100, -50])
axs[0][0].pcolormesh(-model.x,model.y,step.mat,cmap=cmap,vmin=1,vmax=model.nmat+1)
# axs[0][0].set_xlabel('Distance to Center (km)')
axs[0][0].set_xlim([-180,0])
# axs[0][0].set_ylim([-110,60])
# axs[0][0].set_xlim([-110,0])
axs[0][0].set_ylim([-105,45])
axs[0][0].set_ylabel('Depth (km)')
axs[0][1].xaxis.set_major_locator(xmajorLocator)
axs[0][1].pcolormesh(model2.x,model2.y,step2.mat,
        cmap=cmap,vmin=1,vmax=model2.nmat+1)

# axs[0][1].set_xlabel('Distance to Center (km)')
axs[0][1].set_xlim(0,180)
# axs[0][1].set_xlim(0,110)

axs[1][0].yaxis.set_major_locator(ymajorLocator)
axs[1][0].set_xticks([-150, -100, -50])
# axs[1][0].set_xticks([-100, -50])
axs[1][0].pcolormesh(-model3.x,model3.y,step3.mat,cmap=cmap,vmin=1,vmax=model3.nmat+1)
axs[1][0].set_xlabel('Distance to Center (km)')
axs[1][0].set_xlim([-180,0])
# axs[1][0].set_ylim([-110,60])
# axs[1][0].set_xlim([-110,0])
axs[1][0].set_ylim([-105,45])
axs[1][0].set_ylabel('Depth (km)')
axs[1][1].xaxis.set_major_locator(xmajorLocator)
axs[1][1].pcolormesh(model4.x,model4.y,step4.mat,
        cmap=cmap,vmin=1,vmax=model4.nmat+1)

axs[1][1].set_xlabel('Distance to Center (km)')
axs[1][1].set_xlim(0,180)
# axs[1][1].set_xlim(0,110)
    

axs[0][0].text(.05, .95, "t = 650s", transform=axs[0][0].transAxes, ha="left", va="top")
axs[0][1].text(.95, .95, "t = 750s", transform=axs[0][1].transAxes, ha="right", va="top")
axs[1][0].text(.05, .95, "t = 620s", transform=axs[1][0].transAxes, ha="left", va="top")
axs[1][1].text(.95, .95, "t = 700s", transform=axs[1][1].transAxes, ha="right", va="top")
axs[0][0].text(.05, .6, "$T_{c}$ = 30km", transform=axs[0][0].transAxes, ha="left", va="top")
axs[0][1].text(.95, .6, "45km", transform=axs[0][1].transAxes, ha="right", va="top")
axs[1][0].text(.05, .6, "$T_{c}$ = 30km", transform=axs[1][0].transAxes, ha="left", va="top")
axs[1][1].text(.95, .6, "45km", transform=axs[1][1].transAxes, ha="right", va="top")

axs[0][0].text(.05, .1, "$v$ = 15km/s", transform=axs[0][0].transAxes, ha="left", va="top")
axs[1][0].text(.05, .1, "$v$ = 10km/s", transform=axs[1][0].transAxes, ha="left", va="top")

#plot peak rings
axs[0][0].scatter(-step.xmark[trPR_full],
                    step.ymark[trPR_full], c=mass[trPR_full],
                    cmap=cmap2, s=15, linewidths=0)
axs[0][1].scatter(step2.xmark[trPR2_full],
                step2.ymark[trPR2_full], c=mass2[trPR2_full],
                cmap=cmap2, s=15, linewidths=0)
axs[1][0].scatter(-step3.xmark[trPR3_full],
                step3.ymark[trPR3_full], c=mass3[trPR3_full],
                cmap=cmap2, s=15, linewidths=0)
axs[1][1].scatter(step4.xmark[trPR4_full],
                step4.ymark[trPR4_full], c=mass4[trPR4_full],
                cmap=cmap2, s=15, linewidths=0)
# Material boundaries
[axs[0][0].contour(-model.xc,model.yc,step.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
[axs[0][1].contour(model2.xc,model2.yc,step2.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
[axs[1][0].contour(-model3.xc,model3.yc,step3.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
[axs[1][1].contour(model4.xc,model4.yc,step4.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
# axs[0].set_title('30 km')
# axs[1].set_title('45 km')
#plt.show()

fig.savefig('{}/30v45_10v15RingMax.jpg'.format(dirname),dpi=300, bbox_inches = 'tight', pad = 0)
axs[0][0].cla()
axs[0][1].cla()
axs[1][0].cla()
axs[1][1].cla()
