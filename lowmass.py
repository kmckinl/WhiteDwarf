###############################################################################
#   Team Fuzzbots: Colin Gordon, Eric Kopins, Kyle McKinley, Jiahong Xue
#
#   Code to calculate the luminosity of an adiabatic star of some mass
#
#   Routines to calculate central values of star, including integration of 
#   stellar equations for mass, pressure, and luminosity
#
###############################################################################
import numpy as np
import astro_constants as ac
from ode import rk4
from rootfind import bisect
from adiabat import get_rho_and_T
from chabrier_Teff import CBAH
from lowmass_plots_tables import plot_Tr_and_Lr, Lum_Teff_Central_plots

def initial_cond(m,R,delta,mu):
    """
    Sets initial z a distance delta off-center for integration 
    Arguments:
        m := desired mass in g
        R := radius guess in cm
        delta := some number much smaller than radius
        mu := mean molecular weight
    returns r in cm, m in solar mass units, Pc in cgs, and Luminosity in ergs
    and central density, pressure, and luminosity
    """
    #Set the radius for integration equal to delta 
    z = []
    r = delta
    
    #Calculate the central pressure, density, and temp using m and R
    zc = central(m,R,mu)
    #rhocplot.append(zc[1])
    #Tcplot.append(zc[2])
    
    #Calculate initial mass, pressure, and luminosity
    z.append(((4.0*np.pi/3.0)*zc[1]*delta**3.0))
    z.append(zc[0])
    z.append((4.0*np.pi/3.0)*(delta**3.0)*zc[1]*epsnuc(zc[1],zc[2]))

    return r,z,zc
    
def epsnuc(rho,T):
    T_9 = T*1.0E-9
    eps = ((2.4E4)*rho*(0.706**2.0)/(T_9**(2.0/3.0)))*np.exp(-3.380/(T_9**(1.0/3.0)))
    return eps
    
def central(m,R,mu):
    """
    Takes mass and radius and computes central Pressure and 
    central density.
    
    Arguments:
        m := mass in g
        R := radius in cm
        mu := mean molecular weight of composition
    returns Pc, rhoc, and Tc
    """
        
    Pc = 0.77*(ac.G*m**2.0)/(R**4.0)
    rhoc = 5.99*(3*m)/(4*np.pi*(R**3.0))
    Tc = ((mu*ac.mu)/ac.kb)*(ac.G*m/R)
    zc = np.array([Pc,rhoc,Tc])
    return zc
    
    
def composition(q,x,a):
    """
    Takes arrays of charge (electron count), mass fractions,
    and number of nucleons in nucleus
    Arguments:
        q := array of electron count
        x := array of mass fractions
        a := array of nucleon count
    
    Returns:
        Mean Molecular Weight mu and mue
    """
    muilist = np.zeros(len(q))
    muelist = np.zeros(len(q))
    for i in range(len(q)):
        muilist[i] = (x[i]/a[i])
        muelist[i] = ((q[i]*x[i])/a[i])

    muion = sum(muilist)
    mue = sum(muelist)
    mu = 1/(muion+mue)
    return mu
    
def check_lum(R,m,Teff):
    """
    Takes mass (g), radius (cm), and effective temp for radius (K)
        m := mass in g
        R := radius in cm
        Teff := effective temperature for mass m in K
    returns difference between luminosity of R and desired luminosity
    """
    m_arry,P_arry,r_arry,L_arry, Pc, rhoc, Tc = integrate(m,R,xi,delta,eta,mu)
    L = L_arry[-1]
    return L-Lwant(R,Teff)
    
########### Stellar Structure Equations ######################################

def derivatives(r,z,Pc,rhoc,Tc):
    """
    Creates derivatives for Pressure (cgs) and Mass (Msun)
    
    Arguments:
        r := current value of the radius
        z := array containing current values of Mass, Pressure, and Luminosity
        rho := density in g/cm**3.0
    returns:
        dz := array of dm, dP, and dL (with respect to r)
    """
    
    dz = np.zeros(len(z))
    P = z[1]
    rho,T = get_rho_and_T(P,Pc,rhoc,Tc)

    #Fill the array with dm, dP, and dL
    dz[0] = 4*np.pi*(r**2.0)*rho
    dz[1] = -rho*(ac.G*z[0])/(r**2.0)
    dz[2] = dz[0]*epsnuc(rho,T)
    return dz


######## Integration of Stellar Functions #####################################

def get_stepsize(r,z,dz,xi):
    """
    determines step size based on which scale size is smaller 
    between pressure and mass
    
    Arguments:
        r := array of radii in units of Rsun
        z := array of mass in Msun, pressure in cgs, temp in kelvin
        dP := derivative of Pressure with respect to r
        dm := derivative of mass with respect to r
        dL := derivative of luminosity with respect to density
        xi := some constant much smaller than 1
    """
    
    #Set Pressure and Mass Scale Heights
    HP = z[1]/abs(dz[1])
    Hm = z[0]/abs(dz[0])
    HL = z[2]/abs(dz[2])
    return [xi*min(HP,Hm,HL),xi*HP,xi*Hm,xi*HL]
    
def integrate(m,R,xi,delta,eta,mu):
    """
    Takes Central pressure and factors to integrate the pressure and mass to 
    compute radius and mass.
    Arguments:
        m := mass (g)
        R := desired radius (cm)
        xi := scaling value for min of Pressure/Mass Scalings
        delta := initial radius (very small compared to expected radius)
        eta := check value for Pressure to stop integration
        mue := average mass of composition
    returns tuple of Mass, Pressure, and Radius lists
    """
    r,z,zc = initial_cond(m,R,delta,mu)
    Pc = zc[0]
    rhoc = zc[1]
    Tc = zc[2]
    
    mlist = [z[0]]
    Plist = [z[1]]
    rlist = [r]
    Llist = [z[2]] 
    
    Pclist = [zc[0]]
    rhoclist = [zc[1]]
    Tclist = [zc[2]]    

    while z[1] > eta*Pc:
        dz = derivatives(r,z,Pc,rhoc,Tc)
        h = get_stepsize(r,z,dz,xi)
        z4 = rk4(derivatives,r,z,h[0],args=(Pc,rhoc,Tc))
        z = z4
        r += h[0]        
        mlist.append(z[0])
        Plist.append(z[1])
        rlist.append(r)
        Llist.append(z[2])

    return np.array(mlist),np.array(Plist),np.array(rlist),np.array(Llist),Pc,rhoc,Tc   
    
def Lwant(R,Teff):
    Lwant = (4*np.pi)*(R**2.0)*(ac.sigmaSB)*(Teff**4.0)
    return Lwant
    
############## Routine for Computation of Luminosity ###########################

#Lists for plotting
mplot = []
rplot = []
Lplot = []
rhocplot = []
Tcplot = []
#Command to allow for extraction of Teff from chabrier_Teff.py
tbl = CBAH()

#Create arrays for composition
q = np.array([1,2,7])
x = np.array([0.706,0.275,0.019])
a = np.array([1,4,14])
mu = composition(q,x,a)

#Set initial conditions
M = np.array([0.1,0.15,0.2,0.3])            #Mass in solar units for Teff extraction
Mg = M*ac.Msun                              #Convert mass to g for calculation
xi = 0.1
delta = 0.0001*ac.Rsun
eta = 1E-5

#Set bounds for radius for bisection method
Rhigh = 2.0*ac.Rsun
Rlow = 0.01*ac.Rsun

Tefflist = []
for j in range(len(M)):
    Teff = tbl.get_Teff(M[j])
    Tefflist.append(Teff)


for j in range(len(M)): 
    print('M = ',M[j],' Msun')
    Rb = bisect(check_lum,Rlow,Rhigh,args=(Mg[j],Tefflist[j],))
    m_arry,P_arry,r_arry,L_arry, Pc, rhoc, Tc =\
     integrate(Mg[j],Rb,xi,delta,eta,mu)

    m = m_arry[-1]
    r = r_arry[-1]
    L = L_arry[-1]
    P = P_arry[-1]
    rhocplot.append(rhoc)
    Tcplot.append(Tc)
    print('m,r,l,p = ',m/ac.Msun, r/ac.Rsun, L/ac.Lsun, P/Pc)
    mplot.append(m)
    rplot.append(r)
    Lplot.append(L)

#Calculate the Temperature array as a function of r and plot T and L vs R    
m_arry,P_arry,r_arry,L_arry, Pc, rhoc, Tc =\
     integrate(mplot[-1],rplot[-1],xi,delta,eta,mu)
    
T_arry = np.zeros(len(P_arry))
for i in range(len(P_arry)):
    rho,T = get_rho_and_T(P_arry[i],Pc,rhoc,Tc)
    T_arry[i] = T
    
plot_Tr_and_Lr(r_arry/ac.Rsun,T_arry,L_arry/ac.Lsun)

#Plot log-log of Luminosity and Teff
for i in range(len(Lplot)):
    Lplot[i] = Lplot[i]/ac.Lsun
    
Lum_Teff_Central_plots(Lplot,Tefflist,rhocplot,Tcplot)