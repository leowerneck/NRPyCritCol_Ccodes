import numpy as np
# import matplotlib
# matplotlib.rcParams['text.usetex'] = True
# matplotlib.rcParams['font.family'] = 'serif'
# matplotlib.rcParams['font.serif' ] = 'Computer Modern Roman'
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle

def draw_rectangle_and_connecting_lines(ax1_in,ax2_in,rec_left,rec_right,rec_bot,rec_top,connect):

    # Draw rectangle
    con = Rectangle((rec_left,rec_bot),rec_right-rec_left,rec_top-rec_bot,edgecolor='black',facecolor='none',lw=0.5,ls='--')
    ax1_in.add_artist(con)

    # Draw connecting lines from plot 1 to plot 2
    if connect == 'bottom_to_top':
        xys_1 = [(rec_left,rec_top),(rec_right,rec_top)]
        xys_2 = [(rec_left,rec_bot),(rec_right,rec_bot)]
        for i in range(2):
            con = ConnectionPatch(xyA=xys_1[i], xyB=xys_2[i], coordsA="data", coordsB="data",
                                  axesA=ax2_in, axesB=ax1_in,color="black",ls='--',lw=0.5)
            ax1_in.add_artist(con)
    elif connect == 'left_to_right':
        xys_1 = [(rec_left,rec_top),(rec_left,rec_bot)]
        xys_2 = [(rec_right,rec_top),(rec_right,rec_bot)]
        for i in range(2):
            con = ConnectionPatch(xyA=xys_1[i], xyB=xys_2[i], coordsA="data", coordsB="data",
                                  axesA=ax1_in, axesB=ax2_in,color="black",ls='--',lw=0.5)
            ax1_in.add_artist(con)

# Output file name
outfile = "lapse_self_similarity.png"

# Load weak data
t_w,alp_w,sf_w = np.loadtxt("out_weak.dat").T

# Load strong data
t_s,alp_s,sf_s = np.loadtxt("out_strong.dat").T

fig = plt.figure()

linewidth = 1
color_w   = 'blue'
color_s   = 'orange'

# Top panel
ax1 = plt.subplot(2, 1, 1)
plt.grid(ls=':')
plt.plot(t_w,alp_w,lw=linewidth,c=color_w,label=r"$\eta_{\rm weak} = 0.30332394090$")
plt.plot(t_s,alp_s,lw=linewidth,c=color_s,ls='--',label=r"$\eta_{\rm strong} = 0.30332394095$")
plt.legend(loc=2,markerfirst=False)
plt.xlim(-0.4,7.4)
plt.ylim(-0.2,1.2)
plt.xticks([0,1,2,3,4,5,6,7],['0','1','2','3','4','5','6','7'])
plt.yticks([0,0.2,0.4,0.6,0.8,1],['0.0','0.2','0.4','0.6','0.8','1.0'])

# Bottom right panel
xl2 = [+5.6,+6.6]
yl2 = [-0.1,+0.9]
ax2 = plt.subplot(2, 2, 4)
plt.grid(ls=':')
plt.plot(t_w,alp_w,lw=linewidth,c=color_w)
plt.plot(t_s,alp_s,lw=linewidth,c=color_s,ls='--')
plt.xlim(xl2[0],xl2[1])
plt.ylim(yl2[0],yl2[1])
plt.xticks([5.9,6.2,6.5],['5.9','6.2','6.5'])
plt.yticks([0,0.4,0.8],['0.0','0.4','0.8'])

# Bottom left panel
xl3 = [+6.53,+6.59]
yl3 = [-0.05,+0.55]
ax3 = plt.subplot(2, 2, 3)
plt.grid(ls=':')
plt.plot(t_w,alp_w,lw=linewidth,c=color_w)
plt.plot(t_s,alp_s,lw=linewidth,c=color_s,ls='--')
plt.xlim(xl3[0],xl3[1])
plt.ylim(yl3[0],yl3[1])
plt.xticks([6.54,6.56,6.58],['6.54','6.56','6.58'])
plt.yticks([0,0.25,0.5],['0.00','0.25','0.50'])

draw_rectangle_and_connecting_lines(ax1,ax2,xl2[0],xl2[1],yl2[0],yl2[1],'bottom_to_top')
draw_rectangle_and_connecting_lines(ax2,ax3,xl3[0],xl3[1],yl3[0],yl3[1],'left_to_right')

# Set labels
fig.text(0.5, 0.03, r"$t$", ha='center', va='center')
fig.text(0.03, 0.5, r"$\alpha_{\rm central}$", ha='center', va='center', rotation='vertical')

plt.savefig(outfile,dpi=300,bbox_inches='tight',facecolor='white')
plt.close(fig)
