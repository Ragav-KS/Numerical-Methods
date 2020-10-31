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

    if method == 'Euler':
        C =  [  0   ]
        A = [[  0   ]]
        B =  [  1   ]
    elif method == 'Midpoint':
        C =  [  0,      0.5 ]
        A = [[  0,      0   ],
             [  0.5,    0   ]]
        B =  [  0 ,     1   ]
    elif method == '3rd Order':
        C =  [  0,      1/2,    1   ]
        A = [[  0,      0,      0   ],
             [  1/2,    0,      0   ],
             [  -1,     2,      0   ]]
        B =  [  1/6,    2/3,    1/6 ]
    elif method == '4th Order':
        C =  [  0,      1/2,    1/2,    1   ]
        A = [[  0,      0,      0,      0   ],
             [  1/2,    0,      0,      0   ],
             [  0,      1/2,    0,      0   ],
             [  0,      0,      1,      0   ]]
        B =  [  1/6,    1/3,    1/3,    1/6 ]

    return {'s':    len(C),
            'C':    np.array(C),
            'A':    np.array(A),
            'B':    np.array(B)}

def RungeKutta_General(F :list, xi :float, yi :list, h :float, xmax :float, Bt ):
    itr = int((xmax - xi)/h) + 1
    x = np.linspace(xi, xmax, itr)

    result = {}
    yn = yi
    var = len(yn)
    result[xi] = [*reversed(yn)]

    for n in range(itr - 1):
        xn = x[n]

        hk = np.array( [ [0.0] * Bt['s'] ] * var )

        # k_i
        for i in range(Bt['s']):
            xt = xn + Bt['C'][i] * h
            yt = yn.copy()

            for m in range(var):
                yt[m] += Bt['A'][i].dot(hk[m])
            for m in range(var):
                hk[m, i] = h * F[m](xt, *yt)

        # y_{n+1}
        for i in range(var):
            yn[i] = yn[i] + np.array(Bt['B']).dot(hk[i])
        result[x[n+1]] = [*reversed(yn)]
    return result