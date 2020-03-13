########################################################################
#   Team <your team>: <team members>
#   Michigan State University
#   
#   Contains a list of constants used in astronomy.  
#   All are in CGS units unless specified otherwise.
#
########################################################################

import scipy.constants as sc

# dictionary for storing the names, values, and symbols
cgs_constants = {}

# common fractions
onethird = 1.0/3.0
twothird = 2.0*onethird
fourthird = 4.0*onethird
fivethird = 5.0*onethird

# conversion factors
m_to_cm = 1.0/sc.centi
kg_to_g = sc.kilo
J_to_erg = kg_to_g*m_to_cm**2

# physical constants, converted from scipy.constants
c = sc.c * m_to_cm
cgs_constants['speed of light'] = (c,'c')
h = sc.h * J_to_erg  # Planck
cgs_constants['Planck constant'] = (h,'h')
hbar = h/2.0/sc.pi
cgs_constants['Planck constant reduced'] = (hbar,'hbar')
me = sc.m_e * kg_to_g
cgs_constants['electron mass'] = (me,'me')
mu = sc.m_u * kg_to_g
cgs_constants['atomic mass unit'] = (mu,'mu')
kb = sc.k * J_to_erg
cgs_constants['Boltzmann constant'] = (kb,'kb')
G = sc.G * m_to_cm**3 / kg_to_g
cgs_constants['Newton gravitational constant'] = (G,'G')
sigmaSB = sc.constants.Stefan_Boltzmann * J_to_erg / m_to_cm**2
cgs_constants['Stefan Boltzmann constant'] = (sigmaSB,'sigmaSB')

# astrophysical constants, cgs units.  Obtained from the mesa code
Msun = 1.9892e33
cgs_constants['solar mass'] = (Msun,'Msun')
Rsun = 6.9598e10
cgs_constants['solar radius'] = (Rsun,'Rsun')
Lsun = 3.8418e33
cgs_constants['solar luminosity'] = (Lsun,'Lsun')
AU = 1.495978921e13
cgs_constants['astronomical unit'] = (AU,'AU')
pc = 3.261633*9.460528e17
cgs_constants['parsec'] = (pc,'pc')
yr = 3.1558149984e7
cgs_constants['year'] = (yr,'yr')


if __name__ == '__main__':
    """
    prints out a sorted list of the cgs constants we use
    """
    result = [ key for key in cgs_constants ]
    result.sort()
    print('-'*60)
    print('{0:<32}{1:<8}{2:>20}'.format('constant','symbol','cgs value'))
    print('='*60)
    for key in result:
        print('{0:<32}{1[1]:<8}{1[0]:20.12e}'.format(key,cgs_constants[key]))
    print('-'*60)
