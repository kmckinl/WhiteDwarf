########################################################################
#   Team <your team>: <team members>
#   Michigan State University
#   
#   Stores the Teff(M) data from
#       Chabrier, Baraffe, Allard, and Hauschildt. Evolutionary Models 
#       for Very Low-Mass Stars and Brown Dwarfs with Dusty Atmospheres.
#       Astrophys. Jour. 542:464--472, Oct. 2000.
#
########################################################################

from numpy import array, interp

class CBAH:
    """
    Storage of the mass-Teff relation from Chabrier et al. (2000). 
    Provides a simple linear interpolation for intermediate values via the
    function get_Teff.  Units for mass and Teff are solar mass and Kelvin.

    Example usage:
        In [1]: from chabrier_Teff import CBAH

        In [2]: tbl = CBAH()

        In [3]: tbl.M
        Out[3]: array([ 0.1 ,  0.15,  0.2 ,  0.3 ])

        In [4]: tbl.Teff
        Out[4]: array([ 2800.,  3150.,  3300.,  3400.])

        In [5]: tbl.get_Teff(0.1)
        Out[5]: 2800.0

        In [6]: tbl.get_Teff(0.13)
        Out[6]: 3010.0
        
        In [7]: tbl.get_Teff(0.08)
        Out[7]: 2800.0
        
    Notice for In/Out [7] how the value is 'clipped' to the table.
    """
    
    def __init__(self):
        self.M = array([0.1,0.15,0.2,0.3])
        self.Teff = array([2800.0,3150.0,3300.0,3400.0])

    def get_Teff(self, Mwant):
        """
        Convenient lookup for effective temperatures given M; interpolates     
        table using piecewise linear interpolation. Clips to the range of the 
        table.
        
        Example: if tbl = CBAH(), then Teff = tbl.get_Teff(0.13) returns 3010.
        """
        return interp(Mwant,self.M,self.Teff)
