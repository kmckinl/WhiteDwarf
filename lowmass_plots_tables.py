import numpy as np
import matplotlib.pyplot as plt
import astro_constants as ac


def plot_Tr_and_Lr(r,T,L):
    margin = 0.1
    x0,x1,y0,y1 = set_bounds(r,T,margin)
    plt.plot(r,T) #Temperature and radius plot
    plt.xlim(x0,x1)
    plt.ylim(y0,y1)
    plt.title('Temperature as a Function of the Radius')
    plt.xlabel('Radius (R/Rsun)')
    plt.ylabel('Temperature (Kelvin)')
    plt.savefig('T(r).pdf',format='pdf')
    plt.close()
    
    x0,x1,y0,y1 = set_bounds(r,L,margin)
    plt.plot(r,L) #Luminosity and radius plot    
    plt.xlim(x0,x1)
    plt.ylim(y0,y1)
    plt.title('Luminosity as a Function of the Radius')
    plt.xlabel('Radius (R/Rsun)')
    plt.ylabel('Luminosity (L/Lsun)')
    plt.savefig('L(r).pdf',format='pdf')
    plt.close()


def Lum_Teff_Central_plots(L,T,rhoc,Tc):
    margin = 0.1
    x0,x1,y0,y1 = set_bounds(T,L,margin)        
    plt.loglog(T,L)
    plt.xlim(x0,x1)
    plt.ylim(y0,y1)    
    plt.title('Log-log of Luminosity and Effective Temperature')
    plt.gca().invert_xaxis()
    plt.xlabel('Log of Temperature (T/K)')
    plt.ylabel('Log of Luminosity (L/Lsun)')
    plt.savefig('loglog.pdf',format='pdf')
    plt.close()

    x0,x1,y0,y1 = set_bounds(rhoc,Tc,margin)    
    plt.loglog(rhoc,Tc)
    plt.xlim(x1,x0)
    plt.ylim(y0,y1)
    plt.title('Log-Log of Central Density and Central Temperature')
    plt.xlabel('Log of Central Density (rhoc/gc**-3)')
    plt.ylabel('Log of Central Temperature (Tc/K)')
    plt.savefig('loglogrhocTc.pdf',format='pdf')
    plt.close()


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