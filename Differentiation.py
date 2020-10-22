import numpy as np
import Roots
import LinEq 

def DividedDifference(f, x_i, x_i1):
    return (f(x_i1) - f(x_i))/(x_i1 - x_i)

def Direct_fit_Polynomial_Method(f , p :list, x :int):
    A = np.array([[a**i for i in range(0,len(p))] for a in p])
    B = np.array([[f(i)] for i in p])
    Vfit_coeff = LinEq.Gaussian_Elimination(A,B).tolist()
    Vfit_d = Roots.Polynomial([i*Vfit_coeff[i] for i in range(1,len(Vfit_coeff))])
    return Vfit_d(x)