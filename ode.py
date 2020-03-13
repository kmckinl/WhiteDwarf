########################################################################
#   Team Owners: Edward Brown
#   Michigan State University
#   
#   routines for advancing the solution of z'(t) = f(z,t) from t to t+h
#   1. Forward Euler
#   2. Second-order Runge-Kutta
#   3. Fourth-order Runge-Kutta
#   
#   All routines have the same calling sequence:
#   znew = <stepper>(f,t,z,h,args=())
#   Arguments:
#      f     := user-defined function; RHS of z'(t) = f(t,z,...)
#      t     := current value of indendent variable
#      z     := current value of dependent variable
#      h     := desired step
#      args  := tuple of optional arguments to pass to f
########################################################################

def fEuler(f,t,z,h,args=()):
   """Takes one forward euler step.
      znew = fEuler(f,t,z,h,args=())
      Arguments:
         f     := user-defined function; RHS of z'(t) = f(t,z,...)
         t     := current value of indendent variable
         z     := current value of dependent variable
         h     := desired step
         args  := tuple of optional arguments to pass to f
      Returns: znew := value of z at t+h
   """
   return z + h*f(t,z,*args)

def rk2(f,t,z,h,args=()):
   """Takes one second-order Runge-Kutta step.
      znew = rk2(f,t,z,h,args=())
      Arguments:
         f     := user-defined function; RHS of z'(t) = f(t,z,...)
         t     := current value of indendent variable
         z     := current value of dependent variable
         h     := desired step
         args  := tuple of optional arguments to pass to f
      Returns: znew := value of z at t+h
   """
   zp = z + 0.5*h*f(t,z,*args)
   return z + h*f(t+0.5*h,zp,*args)

def rk4(f,t,z,h,args=()):
   """Takes one fourth-order Runge-Kutta step.
      znew = rk4(f,t,z,h,args=())
      Arguments:
         f     := user-defined function; RHS of z'(t) = f(t,z,...)
         t     := current value of indendent variable
         z     := current value of dependent variable
         h     := desired step
         args  := tuple of optional arguments to pass to f
      Returns: znew := value of z at t+h
   """
   k1 = f(t,z,*args)
   k2 = f(t+0.5*h,z+0.5*h*k1,*args)
   k3 = f(t+0.5*h,z+0.5*h*k2,*args)
   k4 = f(t+h,z+h*k3,*args)
   return z + h*(k1+2.0*k2+2.0*k3+k4)/6.0


if (__name__ == '__main__'):
    """
    test routine with z = [cos(2*pi*t/P), sin(2*pi*t/P)
    z' = f = [-2.0*pi/P*sin(2.0*pi*t/P),2.0*pi/P*cos(2.0*pi*t/P)]
    
    written as a series of odes, this is 
        [u,v]' = [-2.0*pi/P * v, 2.0*pi/P * u]
    """
    
    from numpy import pi, sin, cos, zeros
    
    def f(t,z,omega):
        """
        RHS of equation z' = f(t,z,omega) with solution z = [cos(omega*t), sin(omega*t)]
        """
        zp = zeros(2)
        zp[0] = -omega * z[1]
        zp[1] = omega * z[0]
        return zp
    
    def soln(t,omega):
        """
        Analytical solution of z' = f(t,z,omega) with solution z = [cos(omega*t), sin(omega*t)]
        """
        phase = omega*t
        zsoln = zeros(2)
        zsoln[0] = cos(phase)
        zsoln[1] = sin(phase)
        return zsoln

    def do_one(stepper):
        # Period and stepsize
        P = 1.0
        omega = 2.0*pi/P
        h = P/200.0
    
        # initial conditions
        z = zeros(2)
        t = 0.0
        z[0] = 1.0
        z[1] = 0.0
    
        # print a table of the integration: zs is the exact solution, z is the 
        # numerical one
        header = 'step','t','z[0]','zs[0]','z[1]','zs[1]'
        cnt = 0
        # print every 5th line.
        print_interval = 5

        print('{0:>5}{1:>7}  {2:>7}{3:>7}  {4:>7}{5:>7}'.format(*header))
        while t < P:
            z = stepper(f,t,z,h,args=(omega,))
            t += h
            zs = soln(t,omega)
            cnt += 1
            if (cnt % print_interval == 0):
                print('{0:5d}{1:7.3f}  {2:7.3f}{3:7.3f}  {4:7.3f}{5:7.3f}'.format(cnt,t,z[0],zs[0],z[1],zs[1]))

    print('\n====================Forward Euler====================')
    do_one(fEuler)
    print('\n================2nd order Runge-Kutta================')
    do_one(rk2)
    print('\n================4th order Runge-Kutta================')
    do_one(rk4)
    