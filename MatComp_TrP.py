import pySALEPlot as psp
import matplotlib.pyplot as plt
import numpy as np
from numpy import arange,sqrt,ma
from numpy.ma import masked_array
from pylab import figure,arange,colorbar,ma,sqrt
# Need this for the colorbars we will make on the mirrored plot
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors as clr
from pylab import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.signal import savgol_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
maxsep=6.

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
dirname='I30/Mat'
psp.mkdir_p(dirname)

# Open the datafile
tempdir = '/home/xgh/isale_runs_d/output/peak_ring/'
model=psp.opendatfile(tempdir+'30/v15/I30_V15/jdata.dat')
model2=psp.opendatfile(tempdir+'45/v15/I30_V15/jdata.dat')
[model.modelInfo()]
[model2.modelInfo()]
# Set the distance units to km
model.setScale('km')
model2.setScale('km')
# Set up a pylab figure
# fig, axs = plt.subplots(1, 2, sharey=True, figsize=(8,4.5))
fig, axs = plt.subplots(1, 2, sharey=True, figsize=(12, 5))
# Remove horizontal space between axes
fig.subplots_adjust(wspace=0)
# ax = fig.gca()
axs[0].spines['bottom'].set_linewidth(1.6)
axs[0].spines['left'].set_linewidth(1.6)
axs[0].spines['right'].set_linewidth(1.6)
axs[0].spines['top'].set_linewidth(1.6)
axs[1].spines['bottom'].set_linewidth(1.6)
axs[1].spines['left'].set_linewidth(1.6)
axs[1].spines['right'].set_linewidth(1.6)
axs[1].spines['top'].set_linewidth(1.6)
cmap = LinearSegmentedColormap.from_list('mycmap', [ 'gray', 'sandybrown'])
cmap2 = LinearSegmentedColormap.from_list('mycmap2', [ '#ffe5e5', '#e50000'])
cmap = LinearSegmentedColormap.from_list('mycmap', [ '#80755b', '#f0f0f0'])

trPR = [] 
for line in open(tempdir+'30/v15/I30_V15/RingIndex.txt',"r"):
    trPR.append(int(line.rstrip('\n')))

trPR2 = [] # Front facing PR tracers
for line in open(tempdir+'45/v15/I30_V15/RingIndex.txt',"r"):
    trPR2.append(int(line.rstrip('\n')))

tr25PR = []
for line in open(tempdir+'30/v15/I30_V15/25gpaRingIndex.txt',"r"):
    tr25PR.append(int(line.rstrip('\n')))

tr25PR2 = [] # Front facing PR tracers
for line in open(tempdir+'45/v15/I30_V15/25gpaRingIndex.txt',"r"):
    tr25PR2.append(int(line.rstrip('\n')))

model.tracerMassVol()
mass = model.tracerMass
model2.tracerMassVol()
mass2 = model2.tracerMass

xmajorLocator = MultipleLocator(50)
ymajorLocator = MultipleLocator(30)
for i in arange(0,256,5):
    step1=model.readStep('TrP',i)
    trp1 = step1.data[0]/1e9
    step2=model2.readStep('TrP',i)
    trp2 = step2.data[0]/1e9
    step=model.readStep(['Den'],i)
    step2=model2.readStep(['Den'],i)
    # Plot each graph, and manually set the y tick values
    #plot materials
    axs[0].yaxis.set_major_locator(ymajorLocator)
    axs[0].set_xticks([-250, -200, -150, -100, -50])
    axs[0].set_yticks([-100, -50, 0, 50])
    # axs[0].set_xticks([ -100, -50])
    axs[0].pcolormesh(-model.x,model.y,step.mat,cmap=cmap,vmin=1,vmax=model.nmat+1)
    axs[0].set_xlabel('Distance (km)')
    axs[0].set_xlim([-260,0])
    axs[0].set_ylim([-135, 75])
    # axs[0].set_xlim([-175,0])
    # axs[0].set_ylim([-100,50])
    axs[0].set_ylabel('Depth (km)')
    axs[1].xaxis.set_major_locator(xmajorLocator)
    axs[1].pcolormesh(model2.x,model2.y,step2.mat,
            cmap=cmap,vmin=1,vmax=model2.nmat+1)
    
    axs[1].set_xlabel('Distance (km)')
    axs[1].set_xlim(0,260)
    # axs[1].set_xlim(0,175)
    #plot peak rings
    # p = axs[0].scatter(-step.xmark[trPR],
    #                 step.ymark[trPR], c=mass[trPR],
    #                 cmap=cmap2, s=15, linewidths=0)
    # axs[1].scatter(step2.xmark[trPR2],
    #                 step2.ymark[trPR2], c=mass2[trPR2],
    #                 cmap=cmap2, s=15, linewidths=0)
    # # Material boundaries
    [axs[0].contour(-model.xc,model.yc,step.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
    [axs[1].contour(model2.xc,model2.yc,step2.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
    
    # axs[0].text(.05, .95, "t ={: 2.0f} s".format(step.time), transform=axs[0].transAxes, ha="left", va="top")
    axs[1].text(.95, .95, "t ={: 2.0f} s".format(step.time), transform=axs[1].transAxes, ha="right", va="top")

    # if i == 0:
    #     x30 = -step.xmark[tr25PR]
    #     y30 = step.ymark[tr25PR]
    #     x30 = x30.astype(int)
    #     y30 = y30.astype(int)
    #     x45 = step2.xmark[tr25PR2]
    #     y45 = step2.ymark[tr25PR2]
    #     x45 = x45.astype(int)
    #     y45 = y45.astype(int)
    #     axs[0].scatter(x30, y30, c='k', s = 0.8, linewidths=0)
    #     axs[1].scatter(x45, y45, c='k', s = 0.8, linewidths=0)

        # cb=fig.colorbar(p)
        # cb.set_label('Density [kg/m$^3$]')
    #plot grid
    
    for u in range(1,model.tracer_numu):
        tru=model.tru[u]
        
        # Plot the tracers in horizontal lines, every 5 lines
        for l in arange(0,len(tru.xlines),5):
    
            # Get the distances between pairs of tracers in xlines
            dist=get_distances(step,tru.xlines[l])
            # Mask the xmark values if separation too big... means the line won't be connected here
            axs[0].plot(ma.masked_array(-step.xmark[tru.xlines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                    step.ymark[tru.xlines[l]][:-1],
                    c="k",marker='None',linestyle='-',linewidth=0.8)
        
        for l in arange(0,len(tru.ylines),5):
    
            # Get the distances between pairs of tracers in xlines
            dist=get_distances(step,tru.ylines[l])
            # Mask the xmark values if separation too big... means the line won't be connected here
            axs[0].plot(ma.masked_array(-step.xmark[tru.ylines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                    step.ymark[tru.ylines[l]][:-1],
                    c="k",marker='None',linestyle='-',linewidth=0.8)
    
    for u in range(1,model2.tracer_numu):
        tru=model2.tru[u]
        
        # Plot the tracers in horizontal lines, every 5 lines
        for l in arange(0,len(tru.xlines),5):
    
            # Get the distances between pairs of tracers in xlines
            dist=get_distances(step2,tru.xlines[l])
            # Mask the xmark values if separation too big... means the line won't be connected here
            axs[1].plot(ma.masked_array(step2.xmark[tru.xlines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                    step2.ymark[tru.xlines[l]][:-1],
                    c="k",marker='None',linestyle='-',linewidth=0.8)
        
        for l in arange(0,len(tru.ylines),5):
    
            # Get the distances between pairs of tracers in xlines
            dist=get_distances(step2,tru.ylines[l])
            # Mask the xmark values if separation too big... means the line won't be connected here
            axs[1].plot(ma.masked_array(step2.xmark[tru.ylines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                    step2.ymark[tru.ylines[l]][:-1],
                    c="k",marker='None',linestyle='-',linewidth=0.8)


        
    # fig.suptitle('{: 5.2f} s'.format(step.time))
    # axs[0].set_title('30 km, 15 km/s')
    # axs[1].set_title('45 km, 15 km/s')
    #plt.show()
    fig.savefig('{}/30v45Mat-{:05d}.jpg'.format(dirname,i),dpi=300, bbox_inches = 'tight', pad = 0)
    axs[0].cla()
    axs[1].cla()

    rim30 = -478/2
    rim45 = 472/2
    if i > 95:
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
        plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
        values = [1, 5, 20, 25, 30, 50, 100]
        norm = clr.BoundaryNorm(values, ncolors = 256)
        cmap_trp = matplotlib.cm.RdYlBu_r
        axs[0].yaxis.set_major_locator(ymajorLocator)
        axs[0].set_xticks([-250, -200, -150, -100, -50])
        axs[0].set_yticks([-100, -50, 0, 50])
        axs[0].set_xlim([-260,0])
        axs[0].set_ylim([-135, 75])
        axs[1].set_xlim(0,260)
        # ax.pcolormesh(model.x,model.y,step.mat,cmap=cmap_trp,vmin=1,vmax=model.nmat+1)
        [axs[1].contour(model2.xc,model2.yc,step2.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1]]      
        [axs[0].contour(-model.xc,model.yc,step1.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1]]      
        im1 = axs[0].scatter(-step1.xmark, step1.ymark, 
        c=trp1, norm = norm, cmap=cmap_trp, s=30, linewidths=0, alpha = 0.7)
        axs[1].scatter(step2.xmark, step2.ymark, 
        c=trp2, norm = norm, cmap=cmap_trp, s=30, linewidths=0, alpha = 0.7)
        axs[1].text(.95, .95, "t ={: 2.0f} s".format(step.time), transform=axs[1].transAxes, ha="right", va="top")

        for u in range(1,model.tracer_numu):
            tru=model.tru[u]
        
            # Plot the tracers in horizontal lines, every 5 lines
            for l in arange(0,len(tru.xlines),5):
        
                # Get the distances between pairs of tracers in xlines
                dist=get_distances(step,tru.xlines[l])
                # Mask the xmark values if separation too big... means the line won't be connected here
                axs[0].plot(ma.masked_array(-step.xmark[tru.xlines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                        step.ymark[tru.xlines[l]][:-1],
                        c="k",marker='None',linestyle='-',linewidth=0.8)
            
            for l in arange(0,len(tru.ylines),5):
        
                # Get the distances between pairs of tracers in xlines
                dist=get_distances(step,tru.ylines[l])
                # Mask the xmark values if separation too big... means the line won't be connected here
                axs[0].plot(ma.masked_array(-step.xmark[tru.ylines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                        step.ymark[tru.ylines[l]][:-1],
                        c="k",marker='None',linestyle='-',linewidth=0.8)
    
        for u in range(1,model2.tracer_numu):
            tru=model2.tru[u]
            
            # Plot the tracers in horizontal lines, every 5 lines
            for l in arange(0,len(tru.xlines),5):
        
                # Get the distances between pairs of tracers in xlines
                dist=get_distances(step2,tru.xlines[l])
                # Mask the xmark values if separation too big... means the line won't be connected here
                axs[1].plot(ma.masked_array(step2.xmark[tru.xlines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                        step2.ymark[tru.xlines[l]][:-1],
                        c="k",marker='None',linestyle='-',linewidth=0.8)
            
            for l in arange(0,len(tru.ylines),5):
        
                # Get the distances between pairs of tracers in xlines
                dist=get_distances(step2,tru.ylines[l])
                # Mask the xmark values if separation too big... means the line won't be connected here
                axs[1].plot(ma.masked_array(step2.xmark[tru.ylines[l]][:-1],mask=dist > maxsep*tru.d[0]),
                        step2.ymark[tru.ylines[l]][:-1],
                        c="k",marker='None',linestyle='-',linewidth=0.8)

        # axins = axs.inset_axes([0.5, 0.5, 0.47, 0.47], xlim = (-95, -85), ylim = (-5, -10), xticklabels = [], yticklabels = [])
        # axes.indicate_inset_zoom(axins, edgecolor = 'k')
        axins = inset_axes(axs[0], width = '50%', height = '3%', loc = 'upper left', borderpad = 0.6)
        cb1 = plt.colorbar(im1, cax = axins, orientation = 'horizontal')
        # cb1 = fig.colorbar(im1, ax=axs[1], pad = 0.02)
        cb1.set_label(r'$P_{peak}$ (GPa)', fontsize = 12)
        fig.savefig('{}/30v45MatTrP-{:05d}.jpg'.format(dirname,i),dpi=300, bbox_inches = 'tight', pad = 0)
        axs[0].cla()
        axs[1].cla()
#Draw Choosing Time
step=model.readStep(['Den'],70)
step2=model2.readStep(['Den'],85)
# Plot each graph, and manually set the y tick values
#plot materials
axs[0].yaxis.set_major_locator(ymajorLocator)
axs[0].set_xticks([-150, -100, -50])
axs[0].pcolormesh(-model.x,model.y,step.mat,cmap=cmap,vmin=1,vmax=model.nmat+1)
axs[0].set_xlabel('Distance (km)')
axs[0].set_ylim([-95,40])
axs[0].set_xlim([-175,0])
axs[0].set_ylabel('Depth (km)')
axs[1].xaxis.set_major_locator(xmajorLocator)
axs[1].pcolormesh(model2.x,model2.y,step2.mat,
        cmap=cmap,vmin=1,vmax=model2.nmat+1)

axs[1].set_xlabel('Distance (km)')
axs[1].set_xlim(0,175)

axs[0].text(.05, .95, "t = 700s", transform=axs[0].transAxes, ha="left", va="top")
axs[1].text(.95, .95, "850s", transform=axs[1].transAxes, ha="right", va="top")

#plot peak rings
axs[0].scatter(-step.xmark[trPR],
                step.ymark[trPR], c=mass[trPR],
                cmap=cmap2, s=30, linewidths=0)
axs[1].scatter(step2.xmark[trPR2],
                step2.ymark[trPR2], c=mass2[trPR2],
                cmap=cmap2, s=30, linewidths=0)
# Material boundaries
[axs[0].contour(-model.xc,model.yc,step.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
[axs[1].contour(model2.xc,model2.yc,step2.cmc[mat],1,colors='k',linewidths=1) for mat in [0,1,2]]
# axs[0].set_title('30 km')
# axs[1].set_title('45 km')
#plt.show()
fig.savefig('{}/30v45MatMax.jpg'.format(dirname),dpi=300, bbox_inches = 'tight', pad = 0)
axs[0].cla()
axs[1].cla()