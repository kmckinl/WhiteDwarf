########################################################################
#   Team Fuzzbots: Colin Gordon, Eric Kopins, Kyle McKinley, Jiahong Xue
#   Michigan State University
#   
#   Computes the equation of state for a degenerate electron gas.
#   Functions for returning pressure(density, mu_e) and 
#   density(pressure, mu_e) are provided.
#
########################################################################
import numpy as np
import astro_constants as ac

def K(mue):
    """
    Calculate constant K from Equation of State for Degenerate Gas
    Converts every gram to Solar Mass Units
    """
    K = (0.5)*(1.0/5.0)*((3.0/(8.0*np.pi))**(2.0/3.0))*((ac.h**2.0)/(ac.me))*((1/(mue*ac.mu))**(5.0/3.0))
    return K
    
def pressure(rho,mue):
    """
    Arguments:
        rho     :=  mass density (<units>)
        mue     :=  baryon/electron ratio
    
    Returns:
        electron degeneracy pressure (<units>)
    """
    P = K(mue)*(rho**(5.0/3.0))
    return P
    
def density(P,mue):
    """
    Arguments:
        P       :=  electron degeneracy pressure (<units>)
        mue     :=  baryon/electron ratio
        
    Returns:
        mass density (<units>)
    """
    rho = ((1.0/K(mue))*P)**(3.0/5.0)
    return rho
    
def temperature(P,mue):
    """
    Takes Pressure and calculates temperature
    Arguments:
        P       := pressure (<units>)
        mue     := baryon/electron ratio
        
    Returns:
        temperature (Kelvin)
    """
    T = ((mue*ac.mu)/ac.kb)*(P/density(P,mue))
    return T
  
if __name__ == "__main__":
    from test_eos import *
    mue = 2.0
    print('EOS compared to tests/sample_EOS_table.txt...')
    if (test_eos(pressure,density,mue,tolerance=1.0e-10)):
        print(' all values within tolerance')
    else:
        print('FAIL')
