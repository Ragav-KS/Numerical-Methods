import  numpy as np
import itertools as it
import portion as por

import LinEq
import  Roots

def DirectFitPolynomial(f, p:list):
    """Interpolation using Direct Fit Method

    Returns the Polynomial Coefficients

    Args:
        f (Funtion): The function that you want to interpolate on

        Should take one float argument and should return a single float output

        p (list): The x values of choice for interpolation

    Returns:
        List: A list of Coefficients of the Direct Fit Polynomial
    """    
    # Coefficient Matrix
    A = np.array([[a**i for i in range(0,len(p))] for a in p])

    # Constant Matrix
    B = np.array([[f(i)] for i in p])

    # Solve using Gaussian elimination
    coeffs = LinEq.Gaussian_Elimination(A,B)
    
    # Return Coefficients
    return coeffs.tolist()


def LagrangePolynomial(f, P_choice :list, get_LoopFunc = False):
    
    def LagrangePolynomial_LoopFunc(f, p):
        def Lagrange_fit(f, p, x):
            n = len(p)
            f_n = 0
            for i in range(n):
                L_i = 1
                for j in range(n):
                    if i==j: continue
                    L_i = ((x - p[j])/(p[i] - p[j]))*L_i
                f_n = f_n + L_i * f(p[i])
            return f_n
        return lambda x: Lagrange_fit(f, p, x)

    if get_LoopFunc == True:
        return LagrangePolynomial_LoopFunc(f, P_choice)
        
    order = len(P_choice)
    f_L_CoeffsArray = np.array([0] * order)
    for x in range(1, order + 1):
        p = P_choice[:x-1] + P_choice[x:]
        f_L_t = np.array([])
        for i in range(order):
            alpha = 0
            for comb in it.combinations(p ,i):
                beta = 1
                for element in comb:
                    beta = beta * element
                alpha = alpha + beta
            f_L_t = np.insert(f_L_t, 0, (alpha*( (-1)**i )) )
        beta = 1
        for i in p:
            beta = beta * (P_choice[x-1] - i)
        f_L_CoeffsArray = f_L_CoeffsArray + f(P_choice[x-1]) * f_L_t / beta
        
    return f_L_CoeffsArray.tolist()

def LinearSplineMethod(f, xi):
    Spl = por.IntervalDict()
    for i in range(len(xi) - 1):
        m = (f(xi[i+1]) - f(xi[i])) / (xi[i+1] - xi[i])
        Spl[por.closed(xi[i],xi[i+1])] = [ f(xi[i]) - m * xi[i] , m ]
    return lambda x: Roots.Polynomial(Spl[x])(x)