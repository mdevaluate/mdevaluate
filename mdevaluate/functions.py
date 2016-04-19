import numpy as np


def kww(t,A, τ ,β):
    return A*np.exp(-(t/τ)**β)


def kww_1e(A, τ, β):
    return τ*(-np.log(1/(np.e*A)))**(1/β)


def cole_davidson(w, A, b, t0):
    P = np.arctan(w*t0)
    return A * np.cos(P)**b * np.sin(b*P)


def cole_cole(w, A, b, t0):
    return A*(w*t0)**b * np.sin(np.pi*b/2) / (1 + 2*(w*t0)**b * np.cos(np.pi*b/2) + (w*t0)**(2*b))


def hav_neg(w, A, e, d, t0):
    P = np.arctan((w*t0)**d * np.sin(np.pi*d/2) / (1 + (w*t0)**d * np.cos(np.pi*d/2)))
    return A * np.sin(e * P) / (1 + 2*(w*t0)**d * np.cos(np.pi*d/2) + (w*t0)**(2*d))**(e/2)
