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