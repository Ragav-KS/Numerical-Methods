import numpy as np

def Gaussian_Elimination(A, B):
    # initialization
    M = np.concatenate((A,B), axis=1)
    r, c = M.shape[0], (M.shape[1] - 1)

    # Forward Elimination
    for i in range(r):
        p = i + M[i:,i].argmax()
        M[p],M[i] = tuple(M[[i, p]])
        M[i] = M[i]/M[i, i]
        for j in range(i + 1, r):
            M[j] = M[j] - (M[j][i]/M[i][i])*M[i]

    # Back Substitution
    X = np.array([0.0]*c)
    Ax = M[:, :c]
    Cx= M[:, c]
    for i in range(r-1,-1,-1):
        X[i] = (Cx[i] - Ax[i].dot(X))/Ax[i,i]

    # return
    return X