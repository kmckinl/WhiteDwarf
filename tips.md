# More Python tips

## `math` vs `numpy`

Many of you want to import both the `math` and the `numpy` modules.  This is unneccesary and can lead to errors.  For numerical work, forget the `math` module and just use `numpy`. To explain this, look at the following `python` session.

    In [1]: import math as m

    In [2]: import numpy as np

    In [3]: x = np.array([0.1,0.2,0.3])

    In [4]: m.cos(x)
    ---------------------------------------------------------------------------
    TypeError                                 Traceback (most recent call last)
    /Users/edward/Teaching/AST304/<ipython-input-4-d291f21f594a> in <module>()
    ----> 1 m.cos(x)

    TypeError: only length-1 arrays can be converted to Python scalars

    In [5]: np.cos(x)
    Out[5]: array([ 0.99500417,  0.98006658,  0.95533649])

See the difference?  The routines in `math` can only work on scalars, i.e., one number at a time.  The routines in `numpy` are more general: they can work on a whole array at once.

## Avoid `pylab`

The `pylab` module is good for quick calculations at the command line, but you should avoid it when writing codes.  To understand why, I first need to explain what it does.  The `pylab` module simply loads the `numpy` and `matplotlib` modules, but it loads *everything* in those modules into the current namespace. 

Why does this matter?  When you write `import numpy as np`, you are making all of the numpy objects available, but to access them, you have to prepend `np.`.  This means that python distinguishes between a numpy routine `np.func` and your routine `func`. In contrast, when you import `pyplot`, the numpy version of `func` and your version have the same name.  Thus, writing `x=func()` may lead to unexpected results: which `func` is python using?

It's best to simply use the following.

    import numpy as np
    import matplotlib.pyplot as plt

This gives you all of the functionality without the potential name clashes.  If you don't want to write `np.pi` all of the time, you can also do

    from numpy import pi

which makes `pi` accessible without the `np.` prefix.

## Use scientific notation

Write `2.4e4` rather than `2.4 * 10**4`.  Why?  The first version is simply read as a number; the second version is read as an expression, which causes python to first compute `10**4` and then multiply it by `2.4`.  If you have this in a routine that is called often, the second version can potentially increase your runtime&mdash;it is a much more expensive operation.

## Define intermediate variables when appropriate

Suppose you have a formula that uses `T_9`, meaning temperature in units of GigaKelvin.  The code, however, passes `T` in units of Kelvin.  Rather than do the conversion everywhere in the formula, you could write

    def f(T):
        T_9 = T * 1.0e-9    # convert to GigaKelven
        return < some expression in terms of T_9 >
    
This makes the formula much cleaner and easier to read (and therefore easier to debug!).  

More generally, code according to the DRY (don't repeat yourself) principle: if you do the same calculation in several places, encapsulate it into a single variable that is defined in one, and only one, place.
