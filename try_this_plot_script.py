import numpy as np
import matplotlib.pyplot as plt

def set_bounds(x,y,margin):
    width = max(x)-min(x)
    height = max(y)-min(y)
    hpad = margin*width
    vpad = margin*height
    x0 = min(x) - hpad
    x1 = max(x) + hpad
    y0 = min(y) - vpad
    y1 = max(y) + vpad
    return x0,x1,y0,y1

# here is some fake data
Teff = np.linspace(3000.0,3500.0,6)
lgL = np.linspace(-1.0,-0.4,6)

# some sizing to make the plot look nice
charsize = 14   # points
major_ticklength = 0.6*charsize
major_tickwidth = 0.9   # points
minor_ticklength = 0.3*charsize
minor_tickwidth = 0.7   # points

plt.tick_params(axis='both',\
    length=major_ticklength, width=major_tickwidth,which='major')
plt.tick_params(axis='both',\
    length=minor_ticklength,width=minor_tickwidth,which='minor')

# set the limits: we'll add a margin of 10% of the plot width around the data 
# points
margin = 0.1
x0,x1,y0,y1 = set_bounds(Teff,lgL,margin)

# notice the reversal on the x-axis
plt.xlim(x1,x0)
plt.ylim(y0,y1)

# now for some nice text formatting:
#   Text inside the $...$ is set in math mode:  a different font is used, and 
#   all kinds of mathematical symbols are available.  Variables, like T and L, 
#   are set in italic font.  To make 'log' set in upright font, we use the \log 
#   directive (there are also commands \ln, \sin, \cos, \tan...)  Likewise, we 
#   want the subscript 'eff' to be in upright font, so 
#   we use the \mathrm{...} directive, meaning 'typeset in math roman'.  The 
#   command \odot makes a dot with a circle around it, which is the symbol for 
#   the sun.  Because of the use of {} and \, we *MUST* precede the string with 
#   an r, which tells python the string contains special characters.
#
plt.xlabel(r'$T_{\mathrm{eff}}$ (K)',fontsize=charsize)
plt.ylabel(r'$\log(L/L_{\odot})$',fontsize=charsize)

plt.plot(Teff,lgL,\
    color='k',linestyle='-',linewidth=2.0,solid_capstyle='round',\
    solid_joinstyle='miter',marker='o',markerfacecolor='r',markersize=10,\
    markeredgecolor='k',)
plt.show()

