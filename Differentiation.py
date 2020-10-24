def Diff_Polynomial(CoeffList :list):
    return [i*CoeffList[i] for i in range(1,len(CoeffList))]

def DividedDifference(f, x_i, x_i1):
    return (f(x_i1) - f(x_i))/(x_i1 - x_i)
