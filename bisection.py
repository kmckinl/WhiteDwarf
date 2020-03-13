########################################################################
#
#   Edward Brown, Michigan State University
#    
#   Routine to find root of a function via bisection, and a test script.
#
########################################################################

def bisect(f,a,b,tol=1.0e-4,rtol=1.0e-4,itmax=80,args=()):
    """
    Returns the root of the function f(x,*args) == 0 by bisection.
    
    Arguments:
        f   := user defined function
        a,b := the root r is bounded by a, b
        tol := (optional) desired tolerance; bisection stops when |a-b|/|a+b| < tol
        rtol:= (optional) desired tolerance; bisection stops if |f(x)| < rtol
        itmax:= (optional) maximum number of iterations
        args:= tuple of optional arguments to pass to f
    """
    from numpy import abs
        
    # set the boundary
    fa = f(a,*args)
    fb = f(b,*args)
    
    # check that the root is indeed bracketed
    if (fa*fb > 0):
        raise ValueError('root is not bracketed')
        return
    
    # orient interval so that fa < 0, fb > 0; swap if necessary
    if (fa > fb):
       a, b = b, a

    # main loop: c is the midpoint
    c = 0.5*(a+b)
    fc = f(c,*args)
    for i in range(itmax):
        
        # # for testing
        # print '\nround {0}'.format(i)
        # print 'f({0:11.8f}) = {1:11.8f}'.format(a,f(a,*args))
        # print 'f({0:11.8f}) = {1:11.8f}'.format(b,f(b,*args))
                
        # are we done?
        if (abs(b-a) < tol): return c
        if (abs(fc) < rtol): return c
        
        # reset the interval bounds
        if fc > 0:
            b = c
            fb = f(b,*args)
        else:
            a = c
            fa = f(a,*args)
        
        # and set up for the next round
        c = 0.5*(a+b)
        fc = f(c,*args)
    
    # if we get here, we made more than itmax intervals, so raise an error
    raise Exception('too many iterations')


# The following script performs a test of the routine
#
if (__name__ == "__main__"):
    from numpy import cos, pi, sqrt, abs
    
    # test routines
    def f1(x):
        return x**2 - 2.0
    
    def f2(x,r,p):
        return cos(2.0*pi*x/p) - r
    
    r1 = bisect(f1,1.0,2.0)
    print('r1 = {0:11.8f}; compare with {1:11.8f}: diff = {2:11.8f}'.format(r1,sqrt(2.0),abs(1.0-r1/sqrt(2))))
    
    r1 = bisect(f1,1.0,2.0,tol=1.0e-10,rtol=1.0e-10)
    print('r1 = {0:14.12f}; compare with {1:14.12f}: diff = {2:14.12f}'.format(r1,sqrt(2.0),abs(1.0-r1/sqrt(2))))

    r2 = bisect(f2,2.0,4.0,args=(0.0,9.0))
    print('r2 = {1:11.8f}; compare with {1:11.8f}: diff = {2:11.8f}'.format(r2,2.25,abs(1.0-r2/2.25)))
