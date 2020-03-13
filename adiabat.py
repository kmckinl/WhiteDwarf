########################################################################
#   Team <your team>: <team members>
#   Michigan State University
#   
#   Computes the equation of state for an adiabatic ideal gas.
#   Provides a function returning density and temperature 
#   for a given pressure and reference values Pc, rhoc, Tc, as well as
#   the adiabatic index gamma.
#
#   The adiabatic index has a default value gamma = 5.0/3.0.
#
########################################################################

import numpy as np

default_gamma = 5.0/3.0

def get_rho_and_T(P,Pc,rhoc,Tc,gamma=default_gamma):
    """
    Compute density and temperature along an adiabat given a pressure and a
    reference point (Pc,rhoc,Tc).
    
    Arguments:
    P       := pressure
    Pc      := reference pressure
    rhoc    := reference density
    Tc      := reference temperature
    gamma   := adiabatic index, defaults to default_gamma = 5.0/3.0
    
    Returns:
    density, temperature
    """

    # fill these in
    rho = rhoc*(P/Pc)**(1/gamma)
    T = Tc*(P/Pc)**(1.0-(1.0/gamma))

    return rho,T


if __name__ == '__main__':
    """
    Does a quick sanity check on get_rho_and_T: constructs consistent ideal gas
    values for Pc, rhoc, and Tc, calls get_rho_and_T, and then checks that the
    new values also obey the ideal gas EOS.
    """

    import astro_constants as cgs
    
    # set the tolerance: eps is the smallest number such that 1.0 + eps != 1.0
    tol = 4.0*np.finfo(1.0).eps
    
    # set mean molecular weight and fiducial values for Pc, rhoc, and Tc
    mu = 0.615
    Pc = 0.77*cgs.G*cgs.Msun**2/cgs.Rsun**4
    rhoc = 5.99*(3.0*cgs.Msun/4.0/np.pi/cgs.Rsun**3)
    Tc  = Pc/rhoc * mu*cgs.mu/cgs.kb
   
    P_arry = Pc*np.logspace(0.0,-3.0,15)
    print('testing get_rho_and_T')
    print('\n'+'-'*70)
    print('{0:>14}{1:>14}{2:>14}{3:>14}{4:>14}\n'.format(
        'P','rho','T','Pcheck','diff'))
    print('='*70)
    for P in P_arry:
       rho,T = get_rho_and_T(P,Pc,rhoc,Tc)
       Pcheck = rho*T*(cgs.kb/mu/cgs.mu)
       resid = 1.0 - Pcheck/P
       if (abs(resid) > tol):
           print('******** found get_rho_and_T inconsistency ********')
       print('{0:14.6e}{1:14.6e}{2:14.6e}{3:14.6e}{4:14.6e}'.format(
           P, rho, T, Pcheck, resid))
    print('-'*70)
