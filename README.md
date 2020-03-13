Project 4: The low-mass main-sequence
=====================================

A full set of instructions are in the docs folder; click on [this link](./docs/instructions.pdf) for a pdf copy. If you just see a blank square with a link `view raw`, then click the link to download the pdf.

In the directory
----------------

*   `README.md`:  this file

*   `adiabat.py`: The equation of state for an adiabatically stratifed gas.  Low-mass stars are fully convective, so their pressures and temperatures lie anlong an adiabat.  Please complete the function `get_rho_and_T`.  To test your formula, run `python adiabat.py`. You will see

        testing get_rho_and_T

        ----------------------------------------------------------------------
                     P           rho             T        Pcheck          diff

        ======================================================================
        
    followed by a table of numbers.  If your formula is inconsistent, the text
    
        ******** found get_rho_and_T inconsistency ********
        
    will appear.  If there are no warning messages, then your EOS passes the test.

*   `chabrier_Teff.py`: Stores the effective temperatures for the low-mass zero-age main sequence, as computed by Chabrier et al. (2000).  This file is usuable as is: you do not need to edit it, and you should not modify it. Here is an example of how the routines in the file are used to get the effective temperature of a low mass star.

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

    Calling `tbl.get_Teff(m)`, where `m` is a star mass in units of _M_<sub>sun</sub>, returns that star's main-sequence effective temperature.
        
*   `rootfind.py`: contains `bisect`; we'll do a rootfind for where nuclear reaction supply the stars' luminosity, which defines where the star joins the main-sequence. **NB:** Use this version: in addition to being named to avoid a clash with a anaconda package, it also contains a better check for convergence.

*   `astro_constants.py`: contains many useful constants

*   `ode.py`: the solvers for ordinary differential equations.

*   `try_this_plot_script.py`: Some of you might have wondered how I get subscrips and mathematical symbols in the plots I've shown.  Or you might wonder how to tweak the plots to make them look nicer. (And yes, having nice looking plots is _very_ important if you are ever in a position to give a talk when it counts, e.g., as part of a job interview.) If you are curious, have a look at this file, and then run it and examine the output.

