import numpy as np

def TrapezoidalMethod(f, a, b, n):
    I = 0
    for i in np.linspace(a,b,num=n,endpoint=False)[1:]:
        I = I + f(i)
    I = ( b - a ) * ( I * 2 + f(a) + f(b) ) / ( 2 * n)
    return I

def Simpsons13Rule(f, a, b ):
    return (b - a) / 6 * ( f(a)  + 4 * f((a + b) / 2) + f(b) )

def Simpsons13Rule_n(f, a, b, n):
    h = (b - a)/n
    even = []
    odd = []
    for i, v in enumerate(np.linspace(a, b, num= n , endpoint=False)[1:]):
        if not(i % 2):
            even.append(v)
        elif (i % 2):
            odd.append(v)
    I = (b - a) * (f(a) + 4 * sum(map(f, even)) + 2 * sum(map(f, odd)) + f(b))/( 3 * n )
    return I

def GaussQuadrature_nPoint(f, a, b, n):
    X, C = np.polynomial.legendre.leggauss(n)

    g = lambda t: f(t *(b - a)/2 + (b + a)/2) * (b - a)/2

    I = 0
    for i in range(1 - 1, n + 1 - 1):
        I += C[i] + g(X[i])
    return I