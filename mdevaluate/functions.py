import numpy as np

def kww(t,A, τ ,β):
    return A*np.exp(-(t/τ)**β)

def kww_1e(A, τ, β):
    return τ*(-np.log(1/(np.e*A)))**(1/β)

