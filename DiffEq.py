import numpy as np

def EulerMethod(f, xi, yi, xmax, h):
    result = {}
    
    n = int((xmax - xi)/h) + 1
    
    x = np.linspace(xi, xmax, n)
    y = yi
    
    result[0] = {'x': x[0], 'y': y}
    
    for i in range(1, n):
        y = y + f(x[i-1], y) * h
        result[i] = {'x': x[i],'y': y}
    return result

def HeunsMethod(f_xy, xi, yi, h, xmax, its):
    n = int((xmax - xi)/h) + 1

    x = np.linspace(xi, xmax, n)
    y = yi

    result = {}
    result[0] = {'x': x[0], 'y': y}

    for i in range(1, n):
        y_p = y + f_xy(x[i-1], y) * h
        for _ in range(its):
            y_p = y + h * (f_xy(x[i-1], y) + f_xy(x[i], y_p))/2
        y = y_p
        result[i] = {'x': x[i],'y': y}
    return result

def MidPointMethod(f_xy, xi, yi, h, xmax):

    n = int((xmax - xi)/h) + 1

    x = np.linspace(xi, xmax, n)
    y = yi

    result = {}
    result[0] = {'x': x[0], 'y': y}
    
    for i in range(1, n):
        y_1b2 = y + f_xy(x[i-1], y) * h/2
        y = y + f_xy((x[i] + x[i-1])/2, y_1b2) * h
        result[i] = {'x': x[i],'y': y}
    return result

def Butcher_Tableau(method):
    if method == 'midpoint':
        C = [0, 0.5]
        A = [[0, 0], [0.5, 0]]
        B = [ 0 , 1]
    elif method == '4th Order':
        C = [0, 0.5, 0.5, 1]
        A = [[0, 0, 0, 0], [0.5, 0, 0, 0], [0, 0.5, 0, 0], [0, 0, 1, 0]]
        B = [ 1/6 , 1/3, 1/3, 1/6]
    return {'s':len(C),'C': C,'A': A,'B': B}

def RungeKutta_General(F, xi :float, yi :list, h :float, xmax :float, But_t ):
    itr = int((xmax - xi)/h) + 1
    x = np.linspace(xi, xmax, itr)
    result = {}
    yn = yi.copy()
    var = len(yn)
    for n in range(itr - 1):
        xn = x[n]
        hk = np.array([[np.nan]*But_t['s']] * var)
        for i in range(But_t['s']):
            xt = xn + But_t['C'][i] * h
            yt = yn.copy()
            for m in range(var):
                for j in range(i):
                    yt[m] = yt[m] + But_t['A'][i][j] * hk[m, j]
            for m in range(var):
                hk[m, i] = h * F[m](xt, *yt)

        for i in range(var):
            yn[i] = yn[i] + np.array(But_t['B']).dot(hk[i])
        result[x[n+1]] = [*reversed(yn)]
    return result