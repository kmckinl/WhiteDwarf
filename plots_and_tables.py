##########################################################################
#   Team Fuzzbots: Colin Gordon, Eric Kopins, Kyle McKinley, Jiahong Xue
#
#   File to create list of Masses in units solar mass, Radii in units of 
#   solar radius, Central Pressure in cgs units, and Central Pressure in 
#   units of (GM^2)/R^4
#   
#   File creates plot of calculated Masses versus Radii in solar units vs
#   observed masses and radii from Provencal
#
##########################################################################

import numpy as np
import matplotlib.pyplot as plt
import astro_constants as ac
from EOS import density

#m = np.array([0.1,0.2,0.3,0.4,0.5])
#r = np.array([1,2,3,4,5])
#pc = np.array([1,2,3,4,5])

def print_table(m,r,pc):
    """
    prints to screen a formatted version of table 1 in the instructions.
    
    Arguments:
        m := array of masses in units of Msun
        r :=array of radii in units of Rsun
        pc := array of central pressures in cgs units
    """
    header1 = 'M', 'R', 'Pc', 'Pc'
    header2 = '(Msun)', '(Rsun)', '(CGS)', '(GM^2/R^4)'
    header_format_string = '{:>3}{:>10}{:>10}{:>11}'
    header_format_string2 = '{:>4}{:>9}{:>9}{:>14}'
    values_format_string = '{:4.2f}{:10.3f}{:12.3e}{:12.3e}'
    print(header_format_string.format(*header1))
    print(header_format_string2.format(*header2))
    pc2list = []
    Rsol = []
    Msol = []
    for j in range(len(m)):
        Msol.append(m[j]/ac.Msun)
        Rsol.append(r[j]/ac.Rsun)
        pc2 = pc[j]/((ac.G*(ac.Msun**2.0))/(ac.Rsun**4.0))
        pc2list.append(pc2)
    if len(m) == len(pc) == len(r):
        for j in range(len(m)):
            print(values_format_string.format(Msol[j],Rsol[j],pc[j],pc2list[j]))
    else:
        print("FAIL: Radius,Mass, and Pc Arrays aren't of equal length!")
        
def plot_mr(m,r):
    """
    Plots the M-R relation along with the observed masses and radii from
    Provencal et al. (1998). The file is written to 'MR.pdf'.
    
    Arguments:
        m := array of computed masses in units of Msun
        r := array of computed radii in units of Rsun
        filename := if supplied, will produce a pdf
    """
    #Store the observed white dwarfs from provencal
    mobs = np.array([1.03,0.48,0.43,0.65,0.594])
    merr = np.array([0.015,0.045,0.02,0.15,0.012])    
    robs = np.array([0.0074,0.0111,0.0124,0.0127,0.0096])
    rerr = np.array([0.0007,0.0015,0.0005,0.002,0.0096])    
    
    #Convert the Masses and Radii to solar units
    for j in range(len(m)):
        m[j] = m[j]/ac.Msun
        r[j] = r[j]/ac.Rsun
    
    #Plot the data with errorbars
    plt.errorbar(mobs,robs,yerr=merr,xerr=rerr,linestyle='None')
    plt.plot(m,r)
    plt.title('Mass-Radius Relation for White Dwarf with HIPPARCUS Data')
    plt.xlabel('Mass (M/Msun)')
    plt.ylabel('Radius (R/Rsun)')
    plt.savefig('MR.pdf',format='pdf')
    plt.close()
