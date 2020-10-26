import numpy as np

def TrapezoidalMethod(f, a, b, n):
    I = 0
    for i in np.linspace(a,b,num=n,endpoint=False)[1:]:
        I = I + f(i)
    I = ( b - a ) * ( I * 2 + f(a) + f(b) ) / ( 2 * n)
    return I