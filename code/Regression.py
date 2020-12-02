import numpy as np
from typing import Callable

import LinEq as eq

def LinearRegression(f, x :np.ndarray):
    y = np.array(list(map(f, x)))
    n = x.shape[0]
    k2 = ( n * x.dot(y) - x.sum() * y.sum()) / ( n * (x**2).sum() - x.sum()**2 )
    k1 = np.average(y) - k2 * np.average(x)
    return lambda x: k1 + k2 * x

def ExponentialRegression(f, x :np.ndarray, l = None):
    y = np.vectorize(f)(x)
    A_l = lambda l: (y.dot(np.exp(l * x)) / np.sum(np.exp(2 * l * x)))
    f_l = lambda l: (x * y).dot(np.exp(l * x)) - A_l(l) * x.dot(np.exp(2 * l * x))

    if l == None:
        print('find root for the returned function and supply it as l')
        return f_l
    else:
        A = A_l(l)
        return lambda x: A * np.exp( l * x )

def PolynomialRegression(f : Callable, x: np.ndarray, order: int) -> np.ndarray:
    y = np.vectorize(f)(x)
    
    T = []
    for i in range(order * 2 - 1):
        T.append(x**i)
    T = np.array(T)

    A = np.zeros((order, order))
    for i, j in enumerate(T, start= - order + 1):
        diag = np.diag(np.fliplr(A),k=-i)
        diag.setflags(write=True)
        diag.fill(np.sum(j))
    C = T[:order].dot(y).T
    poly = eq.Gaussian_Elimination(A,C[...,None])
    return poly