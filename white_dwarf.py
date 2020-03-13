######################################################################
#   Team Fuzzbots: Colin Gordon, Eric Kopinser, Kyle McKinley, Jiahong Xue
#   Michigan State University
#   
#   Primary code for computation of White Dwarf Properties
#   Computes Pressures, Masses, and Radii of White Dwarf
#
########################################################################

import numpy as np
import astro_constants as ac
from EOS import density, K
from ode import rk4
from bisection import bisect
from plots_and_tables import print_table, plot_mr

def derivatives(r,z,mue):
    """
    Creates derivatives for Pressure (cgs) and Mass (Msun)
    
    Arguments:
        r := current value of the radius
        z := array containing current values of Mass and Pressure
        
    returns:
        dz := array of dm and dP (with respect to r)
    """
    dz = np.zeros(len(z))
    rho = density(z[1],mue) #From the EOS file in terms of pressure and mue
    
    #Fill the array with dm and dP
    dz[0] = 4*np.pi*(r**2.0)*rho
    dz[1] = -rho*(ac.G*z[0])/(r**2.0)
    return dz

def get_stepsize(r,z,dz,xi):
    """
    determines step size based on which scale size is smaller 
    between pressure and mass
    
    Arguments:
        r := array of radii in units of Rsun
        z := array of mass in Msun and pressure in cgs
        dP := derivative of Pressure with respect to r
        dm := derivative of mass with respect to r
        xi := some constant much smaller than 1
    """
    
    #Set Pressure and Mass Scale Heights
    HP = z[1]/abs(dz[1])
    Hm = z[0]/abs(dz[0])
    return [xi*min(HP,Hm),xi*HP,xi*Hm]

def initial_cond(Pc,delta,mue):
    """
    Sets initial z a distance delta off-center for integration 
    Arguments:
        Pc := Central Pressure guess
        delta := some number much smaller than radius
    returns r in cm, m in solar mass units, and Pc in cgs
    """
    z = np.zeros(2)
    r = delta
    rhoc = density(Pc,mue)
    z[0] = ((4.0*np.pi/3.0)*rhoc*delta**3.0)
    z[1] = Pc
    return r,z

def integrate(Pc,xi,delta,eta,mue):
    """
    Takes Central pressure and factors to integrate the pressure and mass to 
    compute radius and mass.
    Arguments:
        Pc := Central Pressure
        xi := scaling value for min of Pressure/Mass Scalings
        delta := initial radius (very small compared to expected radius)
        eta := check value for Pressure to stop integration
        mue := average mass of composition
    returns tuple of Mass, Pressure, and Radius lists
    """
    r,z = initial_cond(Pc,delta,mue)
    
    while z[1] > eta*Pc:
        dz = derivatives(r,z,mue)
        h = get_stepsize(r,z,dz,xi)        
        z4 = rk4(derivatives,r,z,h[0],args=(mue,))
        z = z4
        mlist.append(z[0])
        Plist.append(z[1])
        rlist.append(r)
        Hpplot.append(h[1])
        Hmplot.append(h[2])
        r += h[0]
    return [mlist,Plist,rlist,Hpplot,Hmplot]


def Pscale(m):
    """
    Takes m in units of a solar mass and returns Pscale, a ballpark
    estimate for the central pressure.
    """
    Ps = np.zeros(2)
    fudge_factor = 1000000
    
    #Calibrate the relation between P and M and K
    P0 = ((ac.G**5)/((K(mue)**4)))*(((4*np.pi)/3)**(20/3))*(ac.Msun**(10/3))
    Mw = integrate(P0,xi,delta,eta,mue)[0][-1]
    C = (ac.Msun/Mw)**(10/3)
    
    #Using m, calculate the central Pressure and Low/High bounds
    P = C*P0*((m/ac.Msun)**(10/3))
    Ps[0] = P/fudge_factor
    Ps[1] = P*fudge_factor
    return [P,Ps]


def check_mass(Pc,Mwant):
    """
    Takes guess for central pressure and Desired Mass in solar mass units
    Arguments:
        Pc := Central Pressure
        Mwant := Desired Mass in solar mass units
    returns difference between mass of Pc and desired Mass
    """
    params = integrate(Pc,xi,delta,eta,mue)
    m = params[0][-1]
    return m-Mwant
    
    
################################################################################  

##### Run the program #####
mue = 2.0
M = 0.1*ac.Msun
xi = 0.05
delta = 0.0001*ac.Rsun
eta = 1E-5

#Create lists for radius, mass, pressure, for plotting, and for errors  
rlist = []
mlist = []
Plist = []

#Lists for plotting radius mass and pressure (for table)
rplot = []
mplot = []
Pplot = []

#Lists for the scaling heights of mass and pressure
Hpplot = []
Hmplot = []

#Lists for the calculated errors in mass and radius
rerr = []
merr = []

#Set initial conditions, mass in gs, delta in cm
Mwant = np.array([0.1,0.2,0.3,0.4,0.5])*ac.Msun
xi = 0.05
delta = 0.0001*ac.Rsun
eta = 1E-5
mue = 2.0

#Code to create the white dwarf
for j in range(len(Mwant)):
    #Calculate The Central Pressure for Desired Mass and Integrate for 
    #Mass and Radius
    P, Ps = Pscale(Mwant[j])
    
    #Use bisection method with bounds of pressure to find Central Pressure 
    #and calculate mass and radius once again
    Pc = bisect(check_mass,Ps[0],Ps[1],args=(Mwant[j],))
    Inter = integrate(Pc,xi,delta,eta,mue)
    
    #Append the last values of mass and radius as well as scale heights and
    #calculate the errors from scale heights
    m = Inter[0][-1]
    r = Inter[2][-1]
  
    mplot.append(m)
    rplot.append(r)
    Pplot.append(Pc)

#Make a table of the wanted Masses and calculated Radii (solar units) 
#and Central Pressure in cgs and (GM^2/R^4)
print_table(Mwant,rplot,Pplot)

#Create new lists to plot a wider range of masses
Mwant = np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.9,1.0])*ac.Msun
Mplot = []
Rplot = []
Pcplot = []

#Determine Pc,M, and R for wider range of Masses
for j in range(len(Mwant)):
    P,Ps = Pscale(Mwant[j])
    Pc = bisect(check_mass,Ps[0],Ps[1],args=(Mwant[j],))
    Inter = integrate(Pc,xi,delta,eta,mue)
    M = Inter[0][-1]
    R = Inter[2][-1]
    Mplot.append(M)
    Rplot.append(R)
    Pcplot.append(Pc)
#Plot Masses and Radii (solar units) with data from Provencal
plot_mr(Mplot,Rplot)