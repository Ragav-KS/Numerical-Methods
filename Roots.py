def Polynomial(coeff):
    def Polynomial_proto(x, coeffs):
        coeffs = coeffs.copy()
        # coeffs.reverse()
        o = len(coeffs)
        y = 0
        for i in range(o):
            y += coeffs[i]*x**i
        return y
    return lambda x: Polynomial_proto(x,coeff)

def init_values(func, init_a, maximum = 100, step = 10):
    temp_y = func(init_a)

    if temp_y == 0:
        return (init_a,init_a)
    
    is_neg = (temp_y < 0)
    prev_i = init_a

    for i in range(init_a,maximum,step):
        temp_y = func(i)
        if (is_neg ^ (temp_y < 0)) == True:
            return (prev_i, i)
        prev_i = i

def bisection_method(f_x, a, b, max_i):
    if f_x(a) == 0:
        return a, {}
    if f_x(b) == 0:
        return b, {}
    
    if a > b:
        t = a
        a = b
        b = t
    x_bis = 0
    y = 0

    x_values = {}

    for i in range(max_i):
        x_bis = (a + b)/2
        y = f_x(x_bis)
        x_values[i] = {'a': a, 'b': b, 'x_bis': x_bis, 'y': y}
        if y > 0:
            b = x_bis
        elif y < 0:
            a = x_bis
        elif y == 0:
            return x_bis, x_values

    return x_bis,x_values

def falsePos_method(f_x, a, b, max_i):
    if f_x(a) == 0:
        return (a, {})
    if f_x(b) == 0:
        return (b, {})
    
    if a > b:
        t = a
        a = b
        b = t
    x_false = 0
    y = 0

    result_table = {}

    for i in range(max_i):
        f_a = f_x(a)
        f_b = f_x(b)
        x_false = ((a * f_b) - (b * f_a))/(f_b - f_a)
        y = f_x(x_false)
        result_table[i] = {'a': a, 'b': b, 'x_false': x_false, 'y': y}
        if y > 0:
            b = x_false
        elif y < 0:
            a = x_false
        elif y == 0:
            return (x_false, result_table)

    return (x_false, result_table)

def Newton_Raphson_Method(f, fd,  x = 100, max_i = 100):
    y = f(x)

    result_table = {}
    result_table[0] = {'x': x, 'y': y}
    if y == 0:
        return (x, result_table)

    for i in range(1, max_i):
        x = x - f(x)/fd(x)
        y = f(x)
        result_table[i] = {'x': x, 'y': y}
        if y == 0:
            return (x, result_table)

    return (x, result_table)

def Secant_Method(f, xii, xi, max_i= 100):
    xy = {}

    xy[0] = {'x': xii, 'y': f(xii)}
    xy[1] = {'x': xi, 'y': f(xi)}

    for i in range(1, max_i):
        x = xy[i]['x'] - ( xy[i]['y'] * (xy[i - 1]['x'] - xy[i]['x']))/( xy[i -1]['y'] - xy[i]['y'])
        xy[i + 1] = {'x': x, 'y': f(x)}
        if xy[i + 1]['y'] == 0:
            return xy
    return xy